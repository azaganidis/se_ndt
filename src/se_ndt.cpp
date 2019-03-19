#include <Eigen/QR>
#include <se_ndt/se_ndt.hpp>
#include <numeric>
#define COND_SCORE -10000
#define HIST_BINS1 60
#define HIST_BINS2 60
#define RANGE1 25 
#define MGFEAT 20
#define RANGE2 50 
#define  HIST_ACCURATE_DISTANCE 10//Distance where the histogram is expected to have the same value
#define HIST_FREQ_METER 5 //Recomended 10
#ifdef GROUND_TRUTH_POSE 
       #define NO_LOOP_CLOSE
#endif
#define SIMILARITY_THRESHOLD 0.2
//#include <chrono>
//#include <iomanip>

using namespace std;
int edge_counter=0;
void writeG2O_Vertex(Eigen::Affine3d T, int i)
{
    Eigen::Quaterniond qR(T.rotation().matrix());
    std::ofstream g2("se_loop.g2o", std::ofstream::out|std::ofstream::app);
    g2<<"VERTEX_SE3:QUAT "<<i<<" "<<T(0,3)<<" "<<T(1,3)<<" "<<T(2,3);
    g2<<" "<<qR.coeffs().x();
    g2<<" "<<qR.coeffs().y();
    g2<<" "<<qR.coeffs().z();
    g2<<" "<<qR.coeffs().w();
        g2<<std::endl;
    g2.close();
}
void writeG2O_Edge(Eigen::Affine3d T, int i, int j, Eigen::Matrix<double,7,7> C)
{
    Eigen::Quaterniond qR(T.rotation().matrix());
    std::ofstream g2("se_loop.g2o", std::ofstream::out|std::ofstream::app);
    g2<<"EDGE_SE3:QUAT "<<i<<" "<<j<<" "<<T(0,3)<<" "<<T(1,3)<<" "<<T(2,3);
    g2<<" "<<qR.coeffs().x();
    g2<<" "<<qR.coeffs().y();
    g2<<" "<<qR.coeffs().z();
    g2<<" "<<qR.coeffs().w();
    for(int i=0;i<6;i++)
        for(int j=i;j<6;j++)
            g2<<" "<<C(i,j);
    g2<<std::endl;
    g2.close();
}
void NDTMatch_SE::loadMap(perception_oru::NDTMap **map,std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> input_clouds,float sensor_range)
{
    #pragma omp parallel num_threads(8)
    {
        #pragma omp for
        for(size_t i=0;i<input_clouds.size();i++)
        {
            map[i]->loadPointCloud(*input_clouds[i],sensor_range);
            map[i]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE);
        }
    }
}
Eigen::Matrix<double,6,6> getHes(Eigen::Matrix<double,6,6> Hessian,Eigen::Matrix<double,6,1> score_gradient)
{
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,6,6> > Sol (Hessian);
        Eigen::Matrix<double,6,1> evals = Sol.eigenvalues().real();
        double minCoeff = evals.minCoeff();
        double maxCoeff = evals.maxCoeff();
        if(minCoeff < 0)  //|| evals.minCoeff()) // < 10e-5*evals.maxCoeff()) 
		{
			Eigen::Matrix<double,6,6> evecs = Sol.eigenvectors().real();
			double regularizer = score_gradient.norm();
			regularizer = regularizer + minCoeff > 0 ? regularizer : 0.001*maxCoeff - minCoeff;
			//double regularizer = 0.001*maxCoeff - minCoeff;
			Eigen::Matrix<double,6,1> reg;
			//ugly
			reg<<regularizer,regularizer,regularizer,regularizer,regularizer,regularizer;
			evals += reg;
			Eigen::Matrix<double,6,6> Lam;
			Lam = evals.asDiagonal();
			Hessian = evecs*Lam*(evecs.transpose());
		}
		return Hessian;
}
perception_oru::NDTMap** NDTMatch_SE::loadSavedMap(int index)
{
    perception_oru::NDTMap*** loaded_map =new perception_oru::NDTMap** [resolutions.size()]; 
    for(unsigned int i=0;i<resolutions.size();i++)
    {
        float cur_res=(resolutions[i%resolutions.size()]);
        loaded_map[i]= new perception_oru::NDTMap * [NumInputs];
        for(size_t j=0;j<matcher.NumInputs;j++)
        {

            perception_oru::LazyGrid *grid = new perception_oru::LazyGrid(cur_res);
            grid->semantic_index=j;
            perception_oru::NDTMap *mapTMP=new perception_oru::NDTMap(grid);
            loaded_map[i][j]=mapTMP;
            loaded_map[i][j]->guessSize(0,0,0,size[0],size[1],size[2]);
            loaded_map[i][j]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE);
        }
    }
}
NDTMatch_SE::NDTMatch_SE(initializer_list<float> b,initializer_list<int> c,initializer_list<float> d,int nIn,int max_iter):resolutions(b),resolutions_order(c),size(d),NumInputs(nIn), poses(new pcl::PointCloud<pcl::PointXYZL>), toRVIZ(nIn)
{
	firstRun=true;
    matcher.NumInputs=NumInputs;
    matcher.ITR_MAX =max_iter;
    matcher.step_control=true;
    matcher.n_neighbours=2;
    T.setIdentity();
    Td.setIdentity();
    max_size = *max_element(std::begin(size), std::end(size));

	map=new perception_oru::NDTMap ** [resolutions.size()];
	mapLocal=new perception_oru::NDTMap ** [resolutions.size()];
	mapLocal_prev=new perception_oru::NDTMap ** [resolutions.size()];
	for(unsigned int i=0;i<resolutions.size();i++)
	{
        map[i]= new perception_oru::NDTMap * [NumInputs];
        mapLocal[i]= new perception_oru::NDTMap * [NumInputs];
        mapLocal_prev[i]= new perception_oru::NDTMap * [NumInputs];
        float cur_res=resolutions.at(i);
        for(size_t j=0;j<NumInputs;j++)
        {
            perception_oru::LazyGrid *grid = new perception_oru::LazyGrid(cur_res);
            grid->semantic_index=j;
            map[i][j]=new perception_oru::NDTMap(grid);
            map[i][j]->guessSize(0,0,0,size[0],size[1],size[2]);
            map[i][j]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE);

            perception_oru::LazyGrid *grid2 = new perception_oru::LazyGrid(cur_res);
            grid2->semantic_index=j;
            mapLocal[i][j]=new perception_oru::NDTMap(grid2);
            mapLocal[i][j]->guessSize(0,0,0,size[0],size[1],size[2]);
            mapLocal[i][j]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE);

            perception_oru::LazyGrid *grid3 = new perception_oru::LazyGrid(cur_res);
            grid3->semantic_index=j;
            mapLocal_prev[i][j]=new perception_oru::NDTMap(grid3);
            mapLocal_prev[i][j]->guessSize(0,0,0,size[0],size[1],size[2]);
            mapLocal_prev[i][j]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE);
        }
	}
}
NDTMatch_SE::~NDTMatch_SE()
{
    for(unsigned int i=0;i<resolutions.size();i++)
    {
        for(auto j=0;j<NumInputs;j++)
        {
            delete map[i][j];
            delete mapLocal[i][j];
            delete mapLocal_prev[i][j];
        }
        delete[] map[i];
        delete[] mapLocal[i];
        delete[] mapLocal_prev[i];
    }
    delete[] map;
    delete[] mapLocal;
    delete[] mapLocal_prev;
}

#ifdef GL_VISUALIZE
void NDTMatch_SE::visualize()
{
    thread t1(&NDTMatch_SE::visualize_thread, this);
    t1.detach();
}
void NDTMatch_SE::visualize_thread()
{
    if(viewer==NULL)
    {
        viewer = new NDTViz(true);
        viewer->win3D->start_main_loop_own_thread();
    }
    float occupancy=128;
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
                case 'q':delete viewer;break;
            }
            viewer->plotNDTSAccordingToOccupancy(occupancy, map,std::vector<int>{1}, std::vector<int>{0,1,2,3,4,5,6,7});
            viewer->win3D->keyHitReset();
        }
    }
}
#endif
vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> NDTMatch_SE::getSegmentsFast(pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudIn)
{
	int cloudSize = laserCloudIn->points.size();
	vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud;
	for(int i=0;i<NumInputs;i++)
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT(new pcl::PointCloud<pcl::PointXYZ>);
		laserCloud.push_back(cloudT);
	}
	for(int i=0;i<cloudSize;i++)
	{
            auto pnt=laserCloudIn->points[i];
            int atr = round(pnt.intensity);
            if(NumInputs==1)
                atr=0;
            pcl::PointXYZ point;
            point.x=pnt.x;
            point.y=pnt.y;
            point.z=pnt.z;
            laserCloud[atr]->points.push_back(point);
    }
	return laserCloud;
}

Eigen::Affine3d NDTMatch_SE::matchFaster(Eigen::Affine3d Tinit, pcl::PointCloud<pcl::PointXYZI>::Ptr cloud)
{
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud=getSegmentsFast(cloud);
	for(unsigned int i=0;i<resolutions.size();i++)
		loadMap(mapLocal[i],laserCloud);
	if(!firstRun)
	{
		for(auto i:resolutions_order)
		{
			matcher.current_resolution=resolutions.at(i);
			matcher.match(map[i],mapLocal[i],Tinit,true);
		}
		//std::cout<<getHes(matcher.HessianF,matcher.score_gradientF).inverse()<<std::endl;
	}
	else firstRun=false;
	perception_oru::NDTMap ***mapT;
	mapT=map;
	map=mapLocal;
	mapLocal=mapT;
	return Tinit;
}
Eigen::Affine3d NDTMatch_SE::matchFaster(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud)
{
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud=getSegmentsFast(cloud);
    Eigen::Affine3d Tinit;
    Tinit.setIdentity();

	for(unsigned int i=0;i<resolutions.size();i++)
		loadMap(mapLocal[i],laserCloud);
	if(!firstRun)
	{
		for(auto i:resolutions_order)
		{
			matcher.current_resolution=resolutions.at(i);
			matcher.match(map[i],mapLocal[i],Tinit,true);
		}
		//std::cout<<getHes(matcher.HessianF,matcher.score_gradientF).inverse()<<std::endl;
	}
	else firstRun=false;
	perception_oru::NDTMap ***mapT;
	mapT=map;
	map=mapLocal;
	mapLocal=mapT;
	return Tinit;
}
Eigen::Affine3d NDTMatch_SE::matchFaster_OM(Eigen::Affine3d Tinit, pcl::PointCloud<pcl::PointXYZI>::Ptr cloud)
{
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud=getSegmentsFast(cloud);
	for(unsigned int i=0;i<resolutions.size();i++)
		loadMap(mapLocal[i],laserCloud);
	if(!firstRun)
	{
		for(auto i:resolutions_order)
		{
			matcher.current_resolution=resolutions.at(i);
			matcher.match(map[i],mapLocal[i],Tinit,true);
		}
	}
	else firstRun=false;
    std::vector<thread> tc;
    for(size_t i=0;i<laserCloud.size();i++)
        tc.push_back( std::thread([i,&Tinit, laserCloud](){
                    perception_oru::transformPointCloudInPlace(Tinit, *laserCloud[i]);
                }));
    for(auto& t:tc)t.join();
    for(size_t i=0;i<laserCloud.size();i++)
    {
		for(int rez=0;rez<resolutions.size();rez++)
		{
            map[rez][i]->addPointCloud(Tinit.translation(), *laserCloud[i],0.06, 100, 0.03, 255);
            //map[rez][i]->loadPointCloud(*laserCloud[i],sensor_range);
            map[rez][i]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE,1e9,255,Tinit.translation(),0.01);
        }
    }
	return Tinit;
}
/*void NDTMatch_SE::loadMap(int start_index, int stop_index, pcl::PointXYZL &pose)
{
    perception_oru::NDTMap ***mapT;
	mapT=new perception_oru::NDTMap ** [resolutions.size()];
    for(unsigned int i=0;i<resolutions.size();i++)
    {
        mapT[i]= new perception_oru::NDTMap * [NumInputs];
        #pragma omp parallel num_threads(8)
        {
            #pragma omp for
            for(unsigned int j=0;j<NumInputs;j++)
            {
                perception_oru::LazyGrid* grid = new perception_oru::LazyGrid(resolutions[i]);
                grid->semantic_index=j;
                mapT[i][j]=new perception_oru::NDTMap(grid);
                mapT[i][j]->guessSize(0,0,0,size[0],size[1],size[2]);
                //mapT[i][j]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE);
                mapT[i][j]->loadSaved(start_index, stop_index, pose);
                std::vector<perception_oru::NDTCell*> tmV= mapT[i][j]->getAllCells();//deleted in rviz. If no rviz, comment out.
//                std::cerr<<i<<" "<<j<<" "<<tmV.size()<<"\n";
                toRVIZ.insert(toRVIZ.end(), tmV.begin(), tmV.end());
            }
        }
    }
}*/
/*
int NDTMatch_SE::find_start(pcl::PointCloud<pcl::PointXYZL>::iterator pi, float max_size)
{
    int eI = std::distance(poses->begin(), pi);
    int sI;
    for(sI=eI; sI>0; sI--)
    {
        if( abs(poses->at(sI).x-pi->x)>max_size ||
            abs(poses->at(sI).y-pi->y)>max_size ||
            abs(poses->at(sI).z-pi->z)>max_size)
            break;
    }
    return poses->at(sI).label;
}*/
int NDTMatch_SE::find_start(pcl::PointCloud<pcl::PointXYZL>::iterator pi, float max_size)
{
    pcl::PointCloud<pcl::PointXYZL>::iterator si;
    for(si=pi; si!=poses->begin(); --si)
    {
        if( abs(si->x-pi->x)>max_size || abs(si->y-pi->y)>max_size || abs(si->z-pi->z)>max_size)
            return si->label;
    }
    return 0;
}
bool NDTMatch_SE::matchToSaved(Eigen::Affine3d &Td_, pcl::PointXYZL &pose_end, pcl::PointXYZL &pose_current, int start_index, int iP, Eigen::Matrix<double,7,7> &Ccl)
{
    bool result=true;
    Profiler pr;
    pr.start();
    int stop_index=pose_end.label;
    perception_oru::NDTMap ***mapT;
	mapT=new perception_oru::NDTMap ** [resolutions.size()];
    for(unsigned int i=0;i<resolutions.size();i++)
    {
        mapT[i]= new perception_oru::NDTMap * [NumInputs];
            #pragma omp parallel num_threads(8)
            {
                #pragma omp for
        for(unsigned int j=0;j<NumInputs;j++)
        {
            perception_oru::LazyGrid* grid = new perception_oru::LazyGrid(resolutions[i]);
            grid->semantic_index=j;
            mapT[i][j]=new perception_oru::NDTMap(grid,true);
            mapT[i][j]->guessSize(0,0,0,size[0],size[1],size[2]);
            //mapT[i][j]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE);
            mapT[i][j]->loadSaved(start_index, stop_index, pose_end);
        }
            }
    }
    //pr.check(0);
//    std::cerr<<"LOADED"<<std::endl;
    bool good_registration=true;
    for(auto i:resolutions_order)
    {
        matcher.current_resolution=resolutions.at(i);
        if(!matcher.match(mapT[i],map[i],Td_,true))
            good_registration=true;
        if(Td_.translation().norm()>(pose_current.label-stop_index)*0.5)
        {
            good_registration=false;
        }
    }
//    std::cerr<<"MATCHED. NDTs local: "<<sum_dist_local<<" NDTs loaded: "<<sum_dist_load<<" INDXS "<<start_index<< " "<<stop_index<<std::endl;
    //pr.check(1);
    std::cerr<<"Loop closed. Score: "<<matcher.score_best<< " START: "<<start_index<<" STOP: "<<stop_index<<"\n\t";
    //std::cerr<<"Gfeat: "<<(gfeat[pose_current.label]-gfeat[stop_index]).norm()<<"\n\t";
    if(matcher.score_best<COND_SCORE&&good_registration)
    {
        Ccl=getPoseInformation(Td, mapT,map,true);
//        CovSum=covS[iP]+Ccl;
        for(unsigned int i=0;i<resolutions.size();i++)
            #pragma omp parallel num_threads(8)
            {
                #pragma omp for
                for(unsigned int j=0;j<NumInputs;j++)
                {
                    map[i][j]->transformNDTMap(Td_);
                }
            }
    }
    else
    {
        if(!good_registration)
            std::cerr<<"NOT GOOD. \n\t";
        std::cerr<<"Transform: "<<Td_.translation().transpose()<<"\n\t";
        Td_.setIdentity();
        result=false;
    }
    for(unsigned int i=0;i<resolutions.size();i++)
    {
        #pragma omp parallel num_threads(8)
        {
            #pragma omp for
            for(auto j=0;j<NumInputs;j++)
                delete mapT[i][j];
        }
        delete[] mapT[i];
    }
    delete[] mapT;
    //pr.check(2);
    return result;
}
#include <chrono>
bool fixRotation(Eigen::Affine3d &T){
/*
    auto t=T.matrix();
    Eigen::Vector3d x =t.row(0).head(3);
    Eigen::Vector3d y =t.row(1).head(3);
    Eigen::Vector3d z =t.row(2).head(3);
    double error=x.dot(y);
    Eigen::Vector3d x_ort=x-error/2*y;
    Eigen::Vector3d y_ort=y-error/2*x;
    Eigen::Vector3d z_ort=x_ort.cross(y_ort);
    T.matrix().row(0).head(3)= 1/2*(3-x_ort.dot(x_ort))*x_ort;
    T.matrix().row(1).head(3)= 1/2*(3-y_ort.dot(y_ort))*y_ort;
    T.matrix().row(2).head(3)= 1/2*(3-z_ort.dot(z_ort))*z_ort;
    T.matrix().row(3)<<0,0,0,1;
    */
    Eigen::Matrix3d m=T.rotation().matrix();
    auto qr = m.colPivHouseholderQr();
    //T.matrix().block<3,3>(0,0)=qr.householderQ().setLength(qr.nonzeroPivots()) ;
    //T.matrix().row(3)<<0,0,0,1;
    return qr.isInvertible();
}
Eigen::Affine3d NDTMatch_SE::mapUpdate(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud, bool Does_nothing)
{
    Profiler pr;
    pr.start();
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud=getSegmentsFast(cloud);
	for(unsigned int i=0;i<resolutions.size();i++)
		loadMap(mapLocal[i],laserCloud);
    Eigen::Affine3d T_prev= T;
    Eigen::Matrix<double,7,7> C;
	if(!firstRun)
	{
#ifndef GROUND_TRUTH_POSE 
        for(auto i:resolutions_order)
        {
            matcher.current_resolution=resolutions.at(i);
            matcher.match(mapLocal_prev[i],mapLocal[i],Td,true);
        }
        double score_local = matcher.score_best;
        T=T_prev*Td;
        writeG2O_Edge(Td,num_clouds-1, num_clouds, getPoseInformation(Td,mapLocal_prev, mapLocal,true));

        //MAP REG
        /*
        matcher.current_resolution=resolutions.at(1);
        matcher.match(map[1],mapLocal[1],T,true);
        double score_map = matcher.score_best;

        if(!fixRotation(T) || score_map>score_local || !T.matrix().allFinite()|| !T.inverse().matrix().allFinite())
        {
            T=T_prev*Td;
            std::cerr<<"Match to map not successfull. Using result from local matching only.\n";
        }
        C=getPoseInformation(T, map, mapLocal,true);
        writeG2O_Edge(T_prev.inverse()*T,num_clouds-1, num_clouds, C);
        */
#endif
	}
	else firstRun=false;

    writeG2O_Vertex(T, num_clouds);

    T_prev=T;
    perception_oru::NDTMap ***mapT;
    mapT=mapLocal_prev;
    mapLocal_prev=mapLocal;
    mapLocal=mapT;

    std::vector<thread> tc;
    for(size_t i=0;i<laserCloud.size();i++)
        tc.push_back( std::thread([i,&T_prev, laserCloud](){
                    perception_oru::transformPointCloudInPlace(T_prev, *laserCloud[i]);
                }));
    for(auto& t:tc)t.join();
//    std::cerr<<"START_LOADING"<<std::endl;
    #pragma omp parallel num_threads(8)
    {
        #pragma omp for
        for(size_t i=0;i<laserCloud.size();i++)
        {
            for(int rez=0;rez<resolutions.size();rez++)
            {
                map[rez][i]->addPointCloud(T.translation(), *laserCloud[i],0.06, 100, 0.03, 255);
                //map[rez][i]->loadPointCloud(*laserCloud[i],sensor_range);
                map[rez][i]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE,1e9,255,T.translation(),0.01);

            if(rez==1){//Dont print them for now
                std::vector<perception_oru::NDTCell*> tmV= map[rez][i]->getAllCells();
                toRVIZ[i].insert(toRVIZ[i].end(), tmV.begin(), tmV.end());
            }
                //mapLocal_prev[rez][i]->loadPointCloud(*laserCloud[i],sensor_range);
                //mapLocal_prev[rez][i]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE,1e9,255,T.translation(),0.01);
            }
        }
    }
//    std::cerr<<"END_LOADING"<<std::endl;
    Eigen::Vector3d p_=T.translation();
    pcl::PointXYZL pose_current;
    pose_current.x=p_[0];
    pose_current.y=p_[1];
    pose_current.z=p_[2];
    pose_current.label=num_clouds;

    perception_oru::NDTHistogram *hist = new perception_oru::NDTHistogram (map[1], 1, HIST_BINS1, HIST_BINS2, NumInputs,RANGE1, RANGE2);
    double entropy = hist->calculateEntropy();
//    std::ofstream ofLog;
    //ofLog.open("LOG.txt", std::ofstream::out|std::ofstream::app);
    //ofLog<<entropy<<std::endl;

    if(hists.size()==0){
        hists.push_back(*hist);
        poses->push_back(pose_current);
//        covS.push_back(CovSum);
        Ti.push_back(T);
        num_clouds++;
        return T;
    }
    bool point_dist = (pcl::geometry::distance(poses->back(), pose_current)<HIST_FREQ_METER);
    bool entropy_critirion=entropy<hists.back().ENTROPY;//Minimum Entropy, fartherst from uniform distribution.
    bool insert_=!point_dist;
    if(point_dist&&entropy_critirion)
    {
        poses->points.pop_back();
        hists.pop_back();
        //covS.pop_back();
        Ti.pop_back();
        insert_=true;
    }
    if(insert_)
    {
        hists.push_back(*hist);
        //covS.push_back(CovSum);
        poses->push_back(pose_current);
        Ti.push_back(T);
    }
    ///LOOP CLOSE SEARCH
    std::vector<int> poseIdxSearch;
    std::vector<float> poseDistSearch;
    pose_kdtree.setInputCloud(poses);
    //std::cerr<<poses->size()<<std::endl;
    float search_distance=HIST_ACCURATE_DISTANCE+(num_clouds -last_loop_close_id)*0.05;
//    float search_distance=HIST_ACCURATE_DISTANCE+2*radius;
#ifndef NO_LOOP_CLOSE
    std::cerr<<search_distance<<std::endl;
    if(pose_kdtree.radiusSearch(pose_current,search_distance, poseIdxSearch, poseDistSearch))
    {
        std::vector<double> p_vH(poseIdxSearch.size(),10000);
//        std::vector<double> p_vF(poseIdxSearch.size(),10000);
        #pragma omp parallel num_threads(8)
        {
            #pragma omp for
            for(int i=0;i<poseIdxSearch.size();i++)
            {
                int pi=poseIdxSearch[i];
//                Eigen::Matrix<double,7,7> CD=CovSum-covS[pi];
//                double radiusP=sqrt(CD(0,0)+CD(1,1)+CD(2,2));
                if(pcl::geometry::distance(pose_current, poses->at(pi))>(search_distance))
//                if(pcl::geometry::distance(pose_current, poses->at(pi))>(2*radiusP)+HIST_ACCURATE_DISTANCE)
                    continue;
                if(num_clouds - poses->at(pi).label <200)
                    continue;
                double histogram_score = hists.at(pi).getSimilarity(*hist);
                //double feature_vector_score = (gfeat[pose_current.label]-gfeat[poses->at(pi).label]).norm()/1024;
                double similarity =histogram_score;
                if(num_clouds - last_loop_close_id<100 && similarity>last_loop_close_sim*0.8)
                    continue;
                //p_vF[i]=feature_vector_score;
                p_vH[i]=histogram_score;
            }
        }
        auto p_v=p_vH;
        std::vector<double>::iterator min_it=std::min_element(std::begin(p_v), std::end(p_v));
        int iS=(std::distance(std::begin(p_v), min_it));
        //ofLog<<*min_it<<" "<<p_vH[iS]<<endl;//0.014 0.3
        if(p_vH[iS]<=SIMILARITY_THRESHOLD&&*min_it!=0)
        //if(p_vF[iS]<0.02&&p_vH[iS]<0.35&&*min_it!=0)
        {
            int index_in_poses = poseIdxSearch[iS];
            auto pose_it=poses->begin();
            auto hist_it=hists.begin();
            std::advance(pose_it, index_in_poses);
            std::advance(hist_it, index_in_poses);
            std::cerr<<p_vH[iS]<<endl;
            int loadStartIndex=find_start(pose_it, max_size);
            Eigen::Matrix<double, 7,7> Ccl;
            T_prev=Ti[index_in_poses].inverse()*T;
            bool matched= matchToSaved(T_prev,*pose_it, pose_current, loadStartIndex, index_in_poses, Ccl);
            if(matched)
            {
                T=T_prev*T;
                T_prev=Ti[index_in_poses].inverse()*T;
                writeG2O_Edge(T_prev, poses->at(index_in_poses).label, num_clouds, Ccl);
                last_loop_close_sim=*min_it;
                last_loop_close_id=num_clouds;
                pose_current.x=T.translation()[0];
                pose_current.y=T.translation()[1];
                pose_current.z=T.translation()[2];
                //ofLog<<*min_it<<" "<<pose_current<<std::endl;
                ofstream pose_out("LC",std::ofstream::out|std::ofstream::app);
                pose_out<<pose_current<<" "<<*pose_it<<std::endl;
                pose_out.close();
            }
            std::cerr<<*min_it<<" CLOUD: "<<pose_current<<" TO: "<<*pose_it<<index_in_poses<<" T "<<T_prev.translation().transpose()<<std::endl;
            if(!matched)//Remove the pose from the graph, since it already gave a false.
            {
                poses->erase(pose_it);
                hists.erase(hist_it);
            }
        }
    }
#endif
    num_clouds++;
//    std::cerr<<"NEXT CLOUD"<<std::endl;
    pr.elapsed(0);
	return T;
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

Eigen::Matrix<double,7,6> getJacobian(Eigen::VectorXd v)
{
	Eigen::Matrix<double,7,6> m;
	double x=v(0)/2;
	double y=v(1)/2;
	double z=v(2)/2;
	double cx=cos(x);
	double cy=cos(y);
	double cz=cos(z);
	double sx=sin(x);
	double sy=sin(y);
	double sz=sin(z);
	m(0,3)=0;
	m(1,3)=0;
	m(2,3)=0;
	m(3,3)=-sx*sy*cz+cx*cy*sz;
	m(4,3)=-sx*cy*sz-cx*sy*cz;
	m(5,3)=+cx*cy*cz+sx*sy*sz;
	m(6,3)=-sx*cy*cz+cx*sy*sz;

	m(0,4)=0;
	m(1,4)=0;
	m(2,4)=0;
	m(3,4)=+cx*cy*cz-sx*sy*sz;
	m(4,4)=-cx*sy*sz-sx*cy*cz;
	m(5,4)=-sx*sy*cz-cx*cy*sz;
	m(6,4)=-cx*sy*cz+sx*cy*sz;

	m(0,5)=0;
	m(1,5)=0;
	m(2,5)=0;
	m(3,5)=-cx*sy*sz+sx*cy*cz;
	m(4,5)=+cx*cy*cz+sx*sy*sz;
	m(5,5)=-sx*cy*sz-cx*sy*cz;
	m(6,5)=-cx*cy*sz+sx*sy*cz;

	m=m/2;

	m(0,0)=1;//1
	m(1,0)=0;
	m(2,0)=0;
	m(3,0)=0;
	m(4,0)=0;
	m(5,0)=0;
	m(6,0)=0;

	m(0,1)=0;
	m(1,1)=1;//1	
	m(2,1)=0;
	m(3,1)=0;
	m(4,1)=0;
	m(5,1)=0;
	m(6,1)=0;

	m(0,2)=0;
	m(1,2)=0;
	m(2,2)=1;//1
	m(3,2)=0;
	m(4,2)=0;
	m(5,2)=0;
	m(6,2)=0;

	return m;
}
