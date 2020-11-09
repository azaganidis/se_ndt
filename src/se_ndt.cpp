#include <ndt_map/pointcloud_utils.h>//P
#include <Eigen/QR>
#include <Eigen/StdVector>
#include <thread>
#include <se_ndt/se_ndt.hpp>
#include <numeric>
#include <other.hpp>
#include <profiler.hpp>
#include <vector>
#include <csignal>
#define COND_SCORE -10000
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
	firstRun=true;
    matcher.NumInputs=NumInputs;
    matcher.ITR_MAX =max_iter;
    matcher.step_control=true;
    matcher.n_neighbours=2;
    T.setIdentity();
    Td.setIdentity();
    max_size = *max_element(std::begin(size), std::end(size));

    map=allocateMap(resolutions,size,NumInputs);
    mapLocal=allocateMap(resolutions,size,NumInputs);
    mapLocal_prev=allocateMap(resolutions,size,NumInputs);
}
NDTMatch_SE::~NDTMatch_SE()
{
    for (std::map<int,perception_oru::NDTHistogram*>::iterator hist_it= key_hists.begin();
            hist_it != key_hists.end(); ++hist_it) 
        delete hist_it->second;
    destroyMap(map,resolutions.size(),NumInputs);
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
            viewer->plotNDTSAccordingToOccupancy(occupancy, map,res, sem);
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
void loadSaved(perception_oru::NDTMap*** map,int nRes, int nIn, int startID, int stopID,double *pose)
{
    for(unsigned int i=0;i<nRes;i++)
        for(unsigned int j=0;j<nIn;j++)
            map[i][j]->loadSaved(startID,stopID, poes);
}
bool NDTMatch_SE::matchToSaved(Eigen::Affine3d &Td_, Eigen::Vector3d &pose_ref, int start_index, Eigen::Matrix<double,7,7> &Ccl,int stop_index, int current)
{
    ofstream pose_out("LC.txt",std::ofstream::out|std::ofstream::app);
    bool result=true;
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
            //mapT[i][j]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE);
            mapT[i][j]->loadSaved(start_index, stop_index, pose_ref.data());
        }
            }
    }
    for(auto i:resolutions_order)
    {
        matcher.current_resolution=resolutions.at(i);
        matcher.match(mapT[i],map[i],Td_,true);
    }
    pose_out<<"MAP_REG Score: "<<matcher.score_best<<"REF: "<<stop_index<< " CURR: "<<current<<"\n\t";
    if(matcher.score_best<COND_SCORE)
    {
        Ccl=getPoseInformation(Td_, mapT,map,true);
        /*
        for(unsigned int i=0;i<resolutions.size();i++)
            #pragma omp parallel num_threads(N_THREADS)
            {
                #pragma omp for
                for(unsigned int j=0;j<NumInputs;j++)
                {
                    map[i][j]->transformNDTMap(Td_);
                }
            }
        */
        pose_out<<"MAP_REG_SUCCESS "<<Td_.translation().transpose()<<std::endl;
    }
    else
    {
        pose_out<<"MAP_REG_FAIL "<<Td_.translation().transpose()<<std::endl;
        Td_.setIdentity();
        result=false;
    }
    for(unsigned int i=0;i<resolutions.size();i++)
    {
        #pragma omp parallel num_threads(N_THREADS)
        {
            #pragma omp for
            for(auto j=0;j<NumInputs;j++)
                delete mapT[i][j];
        }
        delete[] mapT[i];
    }
    delete[] mapT;
    return result;
}
#include <chrono>

Eigen::Affine3d NDTMatch_SE::slam(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud)
{
    //Profiler pr;
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud=getSegmentsFast(cloud,NumInputs);
	for(unsigned int i=0;i<resolutions.size();i++)
		loadMap(mapLocal[i],laserCloud);
    Eigen::Matrix<double,7,7> C;
	if(!firstRun)
	{
        for(auto i:resolutions_order)
        {
            matcher.current_resolution=resolutions.at(i);
            matcher.match(mapLocal_prev[i],mapLocal[i],Td,true);
        }
        double score_local = matcher.score_best;
        T=T*Td;
        pose_graph.addPose(T,num_clouds);
        pose_graph.addConstraint(Td,num_clouds-1, num_clouds,
                getPoseInformation(Td,mapLocal_prev, mapLocal,true));
	}
	else{
        pose_graph.addPose(T,num_clouds);
        firstRun=false;
    }
    perception_oru::NDTMap ***mT;
    mT=mapLocal_prev;
    mapLocal_prev=mapLocal;
    mapLocal=mT;

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
                map[rez][i]->addPointCloud(T.translation(), *laserCloud[i],0.06, 100, 0.03, 255);
                map[rez][i]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE,1e9,255,T.translation(),0.01);
            }
        }
    }

    perception_oru::NDTHistogram *hist = new perception_oru::NDTHistogram (mapLocal_prev[1], 1, HIST_BINS1, HIST_BINS2, NumInputs,RANGE1, RANGE2);
    if(key_hists.size()==0){
        key_hists[num_clouds]=hist;
        num_clouds++;
        return T;
    }
    ofstream pose_out("LC.txt",std::ofstream::out|std::ofstream::app);
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
                    const std::lock_guard<std::mutex> lock(result_mutex);
                    pose_out<<"HIST_SCORE_G: "<<histogram_score<<std::endl;
                }
                else
                {
                    const std::lock_guard<std::mutex> lock(result_mutex);
                    pose_out<<"HIST_SCORE_L: "<<histogram_score<<std::endl;
                }
            }
        }
        for(const auto& it:hist_rotations)
        {
            int id = it.first;
            int loadID=find_start(key_hists.find(id), max_size);
            Eigen::Vector3d pose_ref = pose_graph.get(id);
            Eigen::Matrix<double, 7,7> Ccl;
            Eigen::Affine3d TC,TL, TP=pose_graph.getT(id);
            TL.setIdentity();
            //TL.matrix().block<3,3>(0,0)=it.second.matrix().block<3,3>(0,0);
            //TL.matrix().block<3,1>(0,3)=pose_ref-T.translation();
            //tLC=T_prev.inverse();
            bool matched= matchToSaved(TL,pose_ref, loadID,Ccl,id,num_clouds);
            if(matched)
            {
                TC = TP.inverse()*T*TL;
                pose_graph.addConstraint(TC,id, num_clouds,Ccl);
                pose_graph.solve();
                last_loop_close_id=num_clouds;
            }
        }
    }
    if(pose_graph.distance(key_hists.rbegin()->first)>HIST_FREQ_METER)
    {
        key_hists[num_clouds]=hist;
    }
    else
        delete hist;
    num_clouds++;
    //pr.elapsed(0);
    pose_out.close();
	return pose_graph.getT();
}

void NDTMatch_SE::slamSimple(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud)
{
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud=getSegmentsFast(cloud,NumInputs);
	for(unsigned int i=0;i<resolutions.size();i++)
		loadMap(mapLocal[i],laserCloud);
    Eigen::Matrix<double,7,7> C;
	if(!firstRun)
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
        firstRun=false;
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
            load_map(id,map);
            for(auto i:resolutions_order)
            {
                matcher.current_resolution=resolutions.at(i);
                matcher.match(map[i],mapLocal_prev[i],TL,true);
            }
            if(matcher.score_best<COND_SCORE)
            {
                Ccl=getPoseInformation(Td,map, mapLocal_prev,true);
                pose_graph.addConstraint(TL,id, num_clouds,Ccl);
                pose_graph.solve();
                last_loop_close_id=num_clouds;
            }
            destroyMap(map,resolutions.size(),NumInputs);
            map=allocateMap(resolutions,size,NumInputs);
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
