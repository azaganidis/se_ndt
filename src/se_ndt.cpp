#include <ndt_map/pointcloud_utils.h>
#include <Eigen/QR>
#include <Eigen/StdVector>
#include <thread>
#include <se_ndt/se_ndt.hpp>
#include <numeric>
#include <other.hpp>
#include <profiler.hpp>
#include <vector>
#include <csignal>
#define COND_SCORE -30000
#define HIST_BINS1 9
#define HIST_BINS2 9
#define RANGE1 25 
#define RANGE2 50 
#define  HIST_ACCURATE_DISTANCE 10//Distance where the histogram is expected to have the same value
#define HIST_FREQ_METER 5 //Recomended 10
#define SIMILARITY_THRESHOLD 0.03

using namespace std;
using namespace perception_oru;


NDTMatch_SE::NDTMatch_SE(initializer_list<float> b,initializer_list<int> c,initializer_list<float> d,int nIn,int max_iter):resolutions(b),resolutions_order(c),size(d),NumInputs(nIn)
{
    matcher.NumInputs=NumInputs;
    matcher.ITR_MAX =max_iter;
    matcher.step_control=true;
    matcher.n_neighbours=2;
    T.setIdentity();
    Td.setIdentity();
    max_size = *max_element(std::begin(size), std::end(size));

    map_=allocateMap(resolutions,size,NumInputs);
    mapLocal=allocateMap(resolutions,size,NumInputs);
    mapLocal_prev=allocateMap(resolutions,size,NumInputs);
}
NDTMatch_SE::~NDTMatch_SE()
{
    for (std::map<int,perception_oru::NDTHistogram*>::iterator hist_it= key_hists.begin();
            hist_it != key_hists.end(); ++hist_it) 
        delete hist_it->second;
    destroyMap(map_,resolutions.size(),NumInputs);
    destroyMap(mapLocal,resolutions.size(),NumInputs);
    destroyMap(mapLocal_prev,resolutions.size(),NumInputs);
}

#ifdef GL_VISUALIZE
std::thread* NDTMatch_SE::visualize()
{
    std::thread *t1=new std::thread(&NDTMatch_SE::visualize_thread, this);
    //t1.detach();
    return t1;
}
void NDTMatch_SE::visualize_thread()
{
    if(viewer==NULL)
    {
        viewer = new NDTViz(true);
        viewer->win3D->start_main_loop_own_thread();
    }
    float occupancy=0;
    std::vector<int > sem(NumInputs);
    std::iota(sem.begin(),sem.end(),0);
    std::vector<int > res=std::vector<int>{resolutions.size()-1};
    while(viewer->win3D->isOpen())
    {
         usleep(10000);
        if (viewer->win3D->keyHit())
        {
            int key = viewer->win3D->getPushedKey();
            switch (key){
                case '+':occupancy++;break;
                case '-':occupancy--;break;
                case '*':occupancy*=2;break;
                case '/':occupancy/=2;break;
            }
            Eigen::Vector3d cpose = pose_graph.get();
            viewer->setCameraPointingToPoint(cpose(0),cpose(1),cpose(2));
            viewer->plotNDTSAccordingToOccupancy(occupancy, map_,res, sem);
            viewer->addPoseGraph(pose_graph.forGL);
            viewer->win3D->keyHitReset();
            if(key=='q')
            {
                viewer->win3D->win_close();
                delete viewer->win3D;
                delete viewer;
                return;
            }
        }
    }
}
#endif

int NDTMatch_SE::find_start(std::map<int,perception_oru::NDTHistogram*>::iterator pi, float max_size)
{
    Eigen::Vector3d cntr = pose_graph.get(pi->first);
    for(auto si=pi;si!=key_hists.begin();--si)
    {
        Eigen::Vector3d df = cntr-pose_graph.get(si->first);
        if((df.cwiseAbs().array()>max_size).any())
            return si->first;
    }
    return key_hists.begin()->first;
}

Eigen::Matrix<double, 7,7> NDTMatch_SE::getPoseInformation(Eigen::Affine3d T, perception_oru::NDTMap*** m1, perception_oru::NDTMap*** m2, bool inverse=false)
{
	Eigen::MatrixXd Covariance(6,6);
	matcher.covariance(m1,m2,T,resolutions,Covariance, true);
    Eigen::Vector3d ea = T.rotation().matrix().eulerAngles(0,1,2);
    Eigen::MatrixXd J=getJacobian(ea);
    Eigen::MatrixXd JI=Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd>(J).pseudoInverse();
    Eigen::Matrix<double, 7,7> QCov=JI.transpose()*Covariance*JI; 
	return QCov;
}
NDTMap*** NDTMatch_SE::allocateAndLoad(int start_id,int id,double *central_pose){
    perception_oru::NDTMap ***mapT;
	mapT=new perception_oru::NDTMap ** [resolutions.size()];
    for(unsigned int i=0;i<resolutions.size();i++)
    {
      mapT[i]= new perception_oru::NDTMap * [NumInputs];
      #pragma omp parallel num_threads(N_THREADS)
      {
        #pragma omp for
        for(unsigned int j=0;j<NumInputs;j++)
        {
            perception_oru::LazyGrid* grid = new perception_oru::LazyGrid(resolutions[i]);
            grid->semantic_index=j;
            mapT[i][j]=new perception_oru::NDTMap(grid,true);
            mapT[i][j]->guessSize(0,0,0,size[i],size[i],size[i]);
            mapT[i][j]->loadSaved(start_id, id, central_pose);
        }
      }
    }
    return mapT;
}
bool NDTMatch_SE::matchToSaved(NDTMap*** M, Eigen::Affine3d &Td_, Eigen::Matrix<double,7,7> &Ccl,int id)
{
    Eigen::Vector3d pose_ref=pose_graph.get(id);
    int start_index=find_start(key_hists.find(id), max_size);
    ofstream pose_out("LC.txt",std::ofstream::out|std::ofstream::app);
    NDTMap*** R = allocateAndLoad(start_index,id, pose_ref.data());
    
    matcher.n_neighbours=4;
    for(auto i:resolutions_order)
    {
        matcher.current_resolution=resolutions.at(i);
        matcher.match(R[i],M[i],Td_,true);
    }
    matcher.n_neighbours=2;
    bool result=matcher.score_best<COND_SCORE;
    if(result)
    {
        Ccl=getPoseInformation(Td_, R,M,true);
        pose_out<<"Success ";
    }
    else
        pose_out<<"Failure ";
    pose_out<<matcher.score_best<<" F:"<<id<<" T:"<<num_clouds<<" "<<
        Td_.translation().transpose()<<" "<<
        Eigen::Quaterniond(Td_.rotation()).coeffs().transpose()<<std::endl;
    destroyMap(R, resolutions.size(), NumInputs);
    return result;
}
#include <chrono>

void print_pose(ostream& os, Eigen::Affine3d& T){
    os<<T.translation().transpose()<<" ";
    os<<Eigen::Quaterniond(T.rotation()).coeffs().transpose()<<std::endl;
}
void NDTMatch_SE::slam(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud)
{
    ofstream pose_out("LC1.txt",std::ofstream::out|std::ofstream::app);
    //Profiler pr;
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud=getSegmentsFast(cloud,NumInputs);
	for(unsigned int i=0;i<resolutions.size();i++)
		loadMap(mapLocal[i],laserCloud);
    Eigen::Matrix<double,7,7> C;
	if(num_clouds>0)
	{
        Eigen::Affine3d tmpT=T*Td;
        for(auto i:resolutions_order)
        {
            matcher.current_resolution=resolutions.at(i);
            matcher.match(map_[i],mapLocal[i],tmpT,true);
        }
        double score_local = matcher.score_best;
        Td=T.inverse()*tmpT;
        T=tmpT;
        pose_graph.addPose(pose_graph.getT()*Td,num_clouds);
        pose_graph.addConstraint(Td,num_clouds-1, num_clouds,
                getPoseInformation(T,map_, mapLocal,true));
        pose_out<<num_clouds-1<<" "<<num_clouds<<" ";print_pose(pose_out,Td);
	}
	else{
        pose_graph.addPose(Eigen::Affine3d::Identity(),num_clouds);
    }

    std::vector<thread> tc;
    {
        Eigen::Affine3d T_=T;
        for(size_t i=0;i<laserCloud.size();i++)
            tc.push_back( std::thread([i,&T_, laserCloud](){
                        perception_oru::transformPointCloudInPlace(T_, *laserCloud[i]); }));
        for(auto& t:tc)t.join();
    }
    #pragma omp parallel num_threads(N_THREADS)
    {
        #pragma omp for
        for(size_t i=0;i<laserCloud.size();i++)
        {
            for(int rez=0;rez<resolutions.size();rez++)
            {
                map_[rez][i]->addPointCloud(T.translation(), *laserCloud[i],0.06, 100, 0.03, 255);
                map_[rez][i]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE,1e9,255,T.translation(),0.01);
            }
        }
    }

    perception_oru::NDTHistogram *hist = new perception_oru::NDTHistogram (mapLocal[1], 1, HIST_BINS1, HIST_BINS2, NumInputs,RANGE1, RANGE2);
    if(num_clouds==0){
        key_hists[num_clouds]=hist;
        num_clouds++;
        return;
    }
    float search_distance=HIST_ACCURATE_DISTANCE+(num_clouds -last_loop_close_id)*0.10;

    std::vector<int> poseIdxSearch;
    poseIdxSearch=pose_graph.get_in_range<perception_oru::NDTHistogram*>(key_hists,search_distance,100);
    if(poseIdxSearch.size()>0)
    {
        std::map<int, Eigen::Affine3d,std::less<int>,
            Eigen::aligned_allocator<std::pair<const int, Eigen::Affine3d> > > hist_rotations;
        #pragma omp parallel num_threads(N_THREADS)
        {
            #pragma omp for
            for(std::vector<int>::iterator it=poseIdxSearch.begin();
                    it!=poseIdxSearch.end();++it)
            {
                Eigen::Affine3d T_=Eigen::Affine3d::Identity();
                hist->bestFitToHistogram(*key_hists[*it],T_,false);
                double histogram_score = hist->getSimilarity(*key_hists[*it],T_);
                if(histogram_score<=SIMILARITY_THRESHOLD)
                    hist_rotations[*it]=T_;
            }
        }
        for(const auto& it:hist_rotations)
        {
            int id = it.first;
            Eigen::Matrix<double, 7,7> Ccl;
            Eigen::Affine3d TL=Eigen::Affine3d::Identity();
            bool matched= matchToSaved(map_,TL, Ccl,id);
            if(matched)
            {
                pose_out<<id<<" "<<num_clouds<<" ";print_pose(pose_out,TL);
                /*
                transformMap(map_,TL,resolutions.size(),NumInputs);
                Eigen::Affine3d TC= pose_graph.getT(id).inverse()*TL*T;
                pose_out<<"TC: ";print_pose(pose_out,TC);
                pose_graph.addConstraint(TC,id, num_clouds,Ccl);
                pose_graph.constantPose(id);
                pose_graph.solve();
                pose_graph.variablePose(id);
                T=pose_graph.getT();
                last_loop_close_id=num_clouds;
                */
            }
        }
        /*
        {
            int id = it.first;
            Eigen::Vector3d pose_ref = pose_graph.get(id);
            Eigen::Matrix<double, 7,7> Ccl;
            Eigen::Affine3d TC,TL, TP=pose_graph.getT(id);
            TL=pose_graph.getT();
            bool matched= matchToSaved(mapLocal_prev,TL, Ccl,id);
            if(matched)
            {
                TC = TP.inverse()*TL;
                pose_out<<"Tconstraint: ";print_pose(pose_out,TC);
                pose_graph.addConstraint(TC,id, num_clouds,Ccl);
                //pose_graph.constantPose(id);
                pose_graph.solve();
                //pose_graph.variablePose(id);
                last_loop_close_id=num_clouds;
                pose_out<<"Tprev: ";print_pose(pose_out,T);
                TP=pose_graph.getT();
                pose_out<<"Tnew: ";print_pose(pose_out,TP);
                TL=T.inverse()*pose_graph.getT();
                pose_out<<"Td: ";print_pose(pose_out,TL);
                transformMap(map_, TL, resolutions.size(), NumInputs);
                T=pose_graph.getT();
            }
        }
        */
    }
    if(pose_graph.distance(key_hists.rbegin()->first)>HIST_FREQ_METER)
        key_hists[num_clouds]=hist;
    else
        delete hist;
    num_clouds++;
    pose_out.close();
}
void NDTMatch_SE::slamSimple(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud)
{
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud=getSegmentsFast(cloud,NumInputs);
	for(unsigned int i=0;i<resolutions.size();i++)
		loadMap(mapLocal[i],laserCloud);
    Eigen::Matrix<double,7,7> C;
	if(num_clouds>0)
	{
        for(auto i:resolutions_order)
        {
            matcher.current_resolution=resolutions.at(i);
            matcher.match(mapLocal_prev[i],mapLocal[i],Td,true);
        }
        double score_local = matcher.score_best;
        pose_graph.addPose(pose_graph.getT()*Td,num_clouds);
        pose_graph.addConstraint(Td,num_clouds-1, num_clouds,
                getPoseInformation(Td,mapLocal_prev, mapLocal,true));
	}
	else{
        pose_graph.addPose(Eigen::Affine3d::Identity(),num_clouds);
    }
    perception_oru::NDTMap ***mT;
    mT=mapLocal_prev;
    mapLocal_prev=mapLocal;
    mapLocal=mT;

    perception_oru::NDTHistogram *hist = new perception_oru::NDTHistogram (mapLocal_prev[1], 1, HIST_BINS1, HIST_BINS2, NumInputs,RANGE1, RANGE2);
    if(key_hists.size()==0){
        save_map(num_clouds,mapLocal_prev);
        key_hists[num_clouds]=hist;
        num_clouds++;
        return;
    }
    ///LOOP CLOSE SEARCH
    std::vector<int> poseIdxSearch;
    float search_distance=HIST_ACCURATE_DISTANCE+(num_clouds -last_loop_close_id)*0.10;

    poseIdxSearch=pose_graph.get_in_range<perception_oru::NDTHistogram*>(key_hists,search_distance,100);
    if(poseIdxSearch.size()>0)
    {
        std::map<int, Eigen::Affine3d,std::less<int>,
            Eigen::aligned_allocator<std::pair<const int, Eigen::Affine3d> > > hist_rotations;
        std::mutex result_mutex;
        #pragma omp parallel num_threads(N_THREADS)
        {
            #pragma omp for
            for(std::vector<int>::iterator it=poseIdxSearch.begin();
                    it!=poseIdxSearch.end();++it)
            {
                Eigen::Affine3d T_=Eigen::Affine3d::Identity();
                hist->bestFitToHistogram(*key_hists[*it],T_,false);
                double histogram_score = hist->getSimilarity(*key_hists[*it],T_);
                if(histogram_score<=SIMILARITY_THRESHOLD)
                {
                    hist_rotations[*it]=T_;
                }
            }
        }
        for(const auto& it:hist_rotations)
        {
            int id = it.first;
            Eigen::Matrix<double, 7,7> Ccl;
            Eigen::Affine3d TL;
            TL=pose_graph.getT(id).inverse()*pose_graph.getT();
            load_map(id,map_);
            for(auto i:resolutions_order)
            {
                matcher.current_resolution=resolutions.at(i);
                matcher.match(map_[i],mapLocal_prev[i],TL,true);
            }
            if(matcher.score_best<COND_SCORE)
            {
                Ccl=getPoseInformation(Td,map_, mapLocal_prev,true);
                pose_graph.addConstraint(TL,id, num_clouds,Ccl);
                pose_graph.solve();
                last_loop_close_id=num_clouds;
            }
            destroyMap(map_,resolutions.size(),NumInputs);
            map_=allocateMap(resolutions,size,NumInputs);
        }
    }
    if(pose_graph.distance(key_hists.rbegin()->first)>HIST_FREQ_METER)
    {
        save_map(num_clouds,mapLocal_prev);
        key_hists[num_clouds]=hist;
    }
    else
        delete hist;
    num_clouds++;
}
