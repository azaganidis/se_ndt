#include <se_ndt/se_ndt.hpp>
#include <numeric>
#define COND_SCORE -25000
#define COND_SIMILARITY 0.2
#define HIST_BINS1 50
#define HIST_BINS2 50
//#include <chrono>
//#include <iomanip>

using namespace std;

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
NDTMatch_SE::NDTMatch_SE(initializer_list<float> b,initializer_list<int> c,initializer_list<float> d,initializer_list<int> inputs,int max_iter):resolutions(b),resolutions_order(c),size(d),inputs(inputs), poses(new pcl::PointCloud<pcl::PointXYZL>)
{
	NumInputs=inputs.size();
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
    viewer->plotNDTSAccordingToOccupancy(occupancy, map,std::vector<int>{0}, std::vector<int>{0});
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
            viewer->plotNDTSAccordingToOccupancy(occupancy, map,std::vector<int>{0}, std::vector<int>{0});
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
            if(!inputs.at(atr))
                continue;
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
void NDTMatch_SE::loadMap(int start_index, int stop_index, pcl::PointXYZL &pose)
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
                std::vector<perception_oru::NDTCell*> tmV= mapT[i][j]->getAllCells();
//                std::cerr<<i<<" "<<j<<" "<<tmV.size()<<"\n";
                toRVIZ.insert(toRVIZ.end(), tmV.begin(), tmV.end());
            }
        }
    }
}
bool NDTMatch_SE::matchToSaved(Eigen::Affine3d &Td_, int pose_index, pcl::PointXYZL pose_current)
{
    bool result=true;
    Profiler pr;
    pr.start();
    pcl::PointXYZL pose=poses->points.at(pose_index);
    int start_pi=pose_index;
    while(start_pi>=0)
    {
        if( abs(poses->points.at(start_pi).x-pose.x)<max_size &&
            abs(poses->points.at(start_pi).y-pose.y)<max_size &&
            abs(poses->points.at(start_pi).z-pose.z)<max_size)
            start_pi--;
        else
            break;
    }
    int start_index=start_pi>=0?poses->points.at(start_pi).label:0;
    int stop_index=poses->points.at(pose_index).label;
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
            std::vector<perception_oru::NDTCell*> tmV= mapT[i][j]->getAllCells();
            toRVIZ.insert(toRVIZ.end(), tmV.begin(), tmV.end());
        }
            }
    }
    //pr.check(0);
//    std::cerr<<"LOADED"<<std::endl;
    Td_.setIdentity();
    Td_.translation()[0]=(pose.x - pose_current.x);
    Td_.translation()[1]=(pose.y - pose_current.y);
    Td_.translation()[2]=(pose.z - pose_current.z);
    bool good_registration=true;
    for(auto i:resolutions_order)
    {
        matcher.current_resolution=resolutions.at(i);
        if(!matcher.match(mapT[i],map[i],Td_,true))
            good_registration=false;
    }
//    std::cerr<<"MATCHED. NDTs local: "<<sum_dist_local<<" NDTs loaded: "<<sum_dist_load<<" INDXS "<<start_index<< " "<<stop_index<<std::endl;
    //pr.check(1);
    std::cerr<<"Loop closed. Score: "<<matcher.score_best<< " START: "<<start_index<<" STOP: "<<stop_index<<"\n\t";
    if(matcher.score_best<COND_SCORE&&good_registration)
    {
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
        int sc=0;
        for(unsigned int i=0;i<resolutions.size();i++)
            for(unsigned int j=0;j<NumInputs;j++)
            {
                std::vector<perception_oru::NDTCell*> tmV= mapT[i][j]->getAllCells();
                sc+=tmV.size();
                //toRVIZ.insert(toRVIZ.end(), tmV.begin(), tmV.end());
            }
        if(!good_registration)
            std::cerr<<"NOT GOOD. ";
        std::cerr<<"NCELLS "<<sc<<"\n\t";
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
Eigen::Affine3d NDTMatch_SE::mapUpdate(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud, bool constant_velocity)
{
    //Profiler pr;
    //pr.start();
    bool ignore_map=true;
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud=getSegmentsFast(cloud);
	for(unsigned int i=0;i<resolutions.size();i++)
		loadMap(mapLocal[i],laserCloud);
    Eigen::Affine3d T_prev= T;
    if(constant_velocity)
        T=T_prev*Td;
	if(!firstRun)
	{
		for(auto i:resolutions_order)
		{
			matcher.current_resolution=resolutions.at(i);
            if(ignore_map)
                matcher.match(mapLocal_prev[i],mapLocal[i],Td,true);
            else
                matcher.match(map[i],mapLocal[i],T,true);
		}
	}
	else firstRun=false;
    if(ignore_map)
    {
        perception_oru::NDTMap ***mapT;
        mapT=mapLocal_prev;
        mapLocal_prev=mapLocal;
        mapLocal=mapT;
    }
    if(constant_velocity&&!ignore_map)
        Td=T_prev.inverse()*T;
    //Td.translation()<<8,0,0;
    if(ignore_map)
        T=T_prev*Td;
    T_prev=T;
    //std::cout<<matcher.score_best<<std::endl;
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

    perception_oru::NDTHistogram *hist = new perception_oru::NDTHistogram (map[resolutions.size()-1], 1, HIST_BINS1, HIST_BINS2, NumInputs,25, 50);
    double entropy = hist->calculateEntropy();
    std::ofstream ofLog;
    ofLog.open("LOG.txt", std::ofstream::out|std::ofstream::app);

    if(hists.size()>0){
        bool point_dist = (pcl::geometry::distance(poses->back(), pose_current)<10);
        bool entropy_more=hists.back().ENTROPY<entropy;
        bool insert_=!point_dist;
        if(point_dist&&entropy_more)
        {
            //std::cerr<<"POP "<<entropy<<" "<<hists.back().ENTROPY<<std::endl;
            poses->points.pop_back();
            hists.pop_back();
            insert_=true;
        }
        if(insert_)
        {
            hists.push_back(*hist);
            poses->push_back(pose_current);
            //std::cerr<<"INSERT "<<entropy<<" "<<hists.back().ENTROPY<<std::endl;
        }
        ///LOOP CLOSE SEARCH
        std::vector<int> poseIdxSearch;
        std::vector<float> poseDistSearch;
        pose_kdtree.setInputCloud(poses);
        //std::cerr<<poses->size()<<std::endl;
        if(pose_kdtree.radiusSearch(pose_current, max_size, poseIdxSearch, poseDistSearch))
        {
            std::vector<double> p_v(poseIdxSearch.size(),1);
            #pragma omp parallel num_threads(8)
            {
                #pragma omp for
                for(int i=0;i<poseIdxSearch.size();i++)
                {
                    int pi=poseIdxSearch[i];
                    if(num_clouds - poses->at(pi).label <200)
                        continue;
                    double similarity=hists.at(pi).getSimilarity(*hist);
                    if(num_clouds - last_loop_close_id<100 && similarity>last_loop_close_sim*0.8)
                        continue;
                    //std::cerr<<"DIFF "<<num_clouds-last_loop_close_id<<std::endl;
                    p_v[i]=similarity;
                }
            }
            std::vector<double>::iterator min_it=std::min_element(std::begin(p_v), std::end(p_v));
            ofLog<<*min_it<<std::endl;
            if(*min_it<COND_SIMILARITY&&*min_it!=0)
            {
                int index_in_search=(std::distance(std::begin(p_v), min_it));
                int index_in_poses = poseIdxSearch[index_in_search];
                bool matched= matchToSaved(T_prev, index_in_poses, pose_current);
                if(matched)
                {
                    T=T_prev*T;
                    last_loop_close_sim=*min_it;
                    last_loop_close_id=num_clouds;
                    pose_current.x=T.translation()[0];
                    pose_current.y=T.translation()[1];
                    pose_current.z=T.translation()[2];
                    //ofLog<<*min_it<<" "<<pose_current<<std::endl;
                }
                std::cerr<<*min_it<<" CLOUD: "<<pose_current<<" TO: "<<poses->at(index_in_poses)<<" T "<<T_prev.translation().transpose()<<std::endl;
            }
        }
    }
    else
    {
            hists.push_back(*hist);
            poses->push_back(pose_current);
    }
    num_clouds++;
    //pr.elapsed(1);
//    std::cerr<<"NEXT CLOUD"<<std::endl;
	return T;
}

Eigen::Matrix<double, 6,6> NDTMatch_SE::getPoseCovariance(Eigen::Affine3d T)
{
	Eigen::MatrixXd Covariance(6,6);
	Eigen::Matrix<double,6,6> Covariance_;
	matcher.covariance(map,mapLocal,T,resolutions,Covariance);
	Covariance_=Covariance;
	return Covariance_;
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
