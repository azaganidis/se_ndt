#include <se_ndt/ndt_fuser_hmt_se.h>
using namespace std;

namespace lslgeneric {
NDTFuserHMT_SE::NDTFuserHMT_SE(Eigen::Affine3d a,initializer_list<float> b,initializer_list<int> c,initializer_list<float> d,initializer_list<int> e,initializer_list<float> ig,float removeP,int max_iter):Tnow(a),resolutions(b),resolutions_order(c),size(d),tails(e),ignore(ig),removeProbability(removeP)
{
	//unimplemented
	bool preLoad=false;
	firstRun=true;
	beHMT=false;
	visualize=false;
	///////////////
	sensor_pose.setIdentity();
	sensor_range = 100.0;
	prefix="";
	    ctr =0;
	checkConsistency = false;
	    translation_fuse_delta = 0.0;
	    rotation_fuse_delta = 0.0;
	    max_translation_norm = 1.0;
	    max_rotation_norm = M_PI;
	localMapSize<<sensor_range,sensor_range,*(size.begin()+2);
    bool step_control=true;
	    hmt_map_dir="map";
	fuseIncomplete=false;
	vector<int> tails_t(tails);
	NumInputs=count_tails(tails_t)+70*std::count(tails_t.begin(),tails_t.end(),117);

		matcher.NumInputs=NumInputs;
		matcher.ITR_MAX =max_iter;
		matcher.step_control=step_control;
	    matcher2D.ITR_MAX = max_iter;
	    matcher2D.step_control=step_control;
		fAddTimes =NULL;
		fRegTimes=NULL;
/*
	    char fname[1000];
	    snprintf(fname,999,"%s_addTime.txt",prefix.c_str());
	    fAddTimes = fopen(fname,"w");
		*/
	    //std::cout<<"MAP: resolution: "<<resolutions.at(0)<<" size "<<*(size.begin()+0)<<" "<<*(size.begin()+1)<<" "<<*(size.begin()+2)<<" sr "<<sensor_range<<std::endl;

#ifndef NO_NDT_VIZ
        if(visualize_){
          viewer = new NDTViz(visualize);
          viewer->win3D->start_main_loop_own_thread(); // Very very ugly to start it here... FIX ME.
        }
#endif
	if(beHMT)
	cout<<"HMT not implemented"<<endl;
		/*
		if(preLoad) {
			char fname[1000];
			snprintf(fname,999,"%s/%s_map%s.jff",hmt_map_dir.c_str(),prefix.c_str(),boost::to_string(i).c_str());
			std::cerr<<"Loading "<<fname<<std::endl;
			map[i]->loadFromJFF(fname);
		}
		*/
		//else{
	map=new lslgeneric::NDTMap ** [resolutions.size()];
	mapLocal=new lslgeneric::NDTMap ** [resolutions.size()];
	for(auto i=0;i<resolutions.size();i++)
	{
		map[i]=initMap(NumInputs,{resolutions.at(i)},size);
		mapLocal[i]=initMap(NumInputs,{resolutions.at(i)},size);
	}
		//}
	Tlast_fuse = Tnow;
	Todom = Tnow;
	if(visualize) 
	{
#ifndef NO_NDT_VIZ
	  //      # error compiling with visualization
		viewer->plotNDTSAccordingToOccupancy(-1,map[0]); 
		//viewer->plotLocalNDTMap(cloud,resolution);
		viewer->addTrajectoryPoint(Tnow.translation()(0),Tnow.translation()(1),Tnow.translation()(2)+0.5,1,0,0);
		viewer->addTrajectoryPoint(Todom.translation()(0),Todom.translation()(1),Todom.translation()(2)+0.5,0,1,0);
		viewer->displayTrajectory();
		//viewer->setCameraPointing(Tnow.translation()(0),Tnow.translation()(1),Tnow.translation()(2)+3);
		viewer->repaint();	
#endif
	}
}


    /**
     *
     *
     */
Eigen::Affine3d NDTFuserHMT_SE::update(Eigen::Affine3d Tmotion, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,initializer_list<vector<double> > attributes)
{
	Todom = Todom * Tmotion; //we track this only for display purposes!
	double t0=0,t1=0,t2=0,t3=0,t4=0,t5=0,t6=0;
	///Set the cloud to sensor frame with respect to base
	lslgeneric::transformPointCloudInPlace(sensor_pose, *cloud);
	t0 = getDoubleTime();
	///Create local map
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud=getSegments(cloud,attributes,tails,ignore,removeProbability);
	for(int i=0;i<resolutions.size();i++)
		loadMap(mapLocal[i],laserCloud,NumInputs);
	//pass through ndlocal and set all cells with vertically pointing normals to non-gaussian :-O
	/*SpatialIndex *index = ndlocal.getMyIndex();
	  typename SpatialIndexctorItr it = index->begin();
	  while (it != index->end())
	  {
	  NDTCell *cell = dynamic_cast<NDTCell*> (*it);
	  if(cell!=NULL)
	  {
	  if(cell->hasGaussian_)
	  {
	  if(cell->getClass() == NDTCell::HORIZONTAL) {
	  cell->hasGaussian_ = false;
	  }
	  }
	  }
	  it++;
	  }*/

	t1 = getDoubleTime();
	Eigen::Affine3d Tinit = Tnow * Tmotion;
	t2 = getDoubleTime();
	bool matcher_res=true;
	if(!firstRun)
		for(auto i:resolutions_order)
		{
			matcher.current_resolution=resolutions.at(i);
			matcher_res&=matcher.match(map[i],mapLocal[i],Tinit,true);
		}
	//std::cout<<getHes(matcher.HessianF,matcher.score_gradientF).inverse()<<std::endl;
	if(matcher_res || fuseIncomplete)
	{
	t3 = getDoubleTime();
	Eigen::Affine3d diff = (Tnow * Tmotion).inverse() * Tinit;

	if((diff.translation().norm() > max_translation_norm || 
			diff.rotation().eulerAngles(0,1,2).norm() > max_rotation_norm) && checkConsistency)
	{
		fprintf(stderr,"****  NDTFuserHMT_L -- ALMOST DEFINATELY A REGISTRATION FAILURE *****\n");
		Tnow = Tnow * Tmotion;
		//save offending map:
		//map->writeToJFF("map.jff");
		//ndlocal.writeToJFF("local.jff");
	}
	else
	{
		Tnow = Tinit;
		//std::cout<<Tinit.translation().norm()<<std::endl;
		//Tnow = Tnow * Tmotion;
		for(int j=0;j<NumInputs;j++)
			lslgeneric::transformPointCloudInPlace(Tnow, *laserCloud[j]);
		Eigen::Affine3d spose = Tnow*sensor_pose;
		Eigen::Affine3d diff_fuse = Tlast_fuse.inverse()*Tnow;
		if(diff_fuse.translation().norm() > translation_fuse_delta ||
			diff_fuse.rotation().eulerAngles(0,1,2).norm() > rotation_fuse_delta ||firstRun)
		{
			firstRun=false;
		//std::cout<<"F: "<<spose.translation().transpose()<<" "<<spose.rotation().eulerAngles(0,1,2).transpose()<<std::endl;
		t4 = getDoubleTime();

		for(unsigned int i=0;i<resolutions.size();i++)
			for(int j=0;j<NumInputs;j++)
			{
				map[i][j]->addPointCloudMeanUpdate(spose.translation(),*laserCloud[j],localMapSize, 1e5, 1250, *(size.begin()+2)/2, 0.06);
//		map[i][j]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE, 1e5, 1250, spose.translation(), 0.06);
				map[i][j]->computeNDTCells();
			}
		t5 = getDoubleTime();
		//map->addPointCloudMeanUpdate(spose.translation(),cloud,localMapSize, 1e5, 25, 2*map_size_z, 0.06);
		//map->addPointCloud(spose.translation(),cloud, 0.06, 25);
		//map->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE, 1e5, 255, spose.translation(), 0.1);
		//t4 = getDoubleTime();
		//std::cout<<"match: "<<t3-t2<<" addPointCloud: "<<t5-t4<<" ndlocal "<<t1-t0<<" total: "<<t5-t0<<std::endl;
		Tlast_fuse = Tnow;
		if(visualize) //&&ctr%20==0) 
		{
#ifndef NO_NDT_VIZ
			if(ctr%2==0) {
			viewer->plotNDTSAccordingToOccupancy(-1,map[0]); 
			//viewer->plotLocalNDTMap(cloud,resolution); 
			}
			viewer->addTrajectoryPoint(Tnow.translation()(0),Tnow.translation()(1),Tnow.translation()(2)+0.2,0,1,0);
			viewer->addTrajectoryPoint(Todom.translation()(0),Todom.translation()(1),Todom.translation()(2)+0.2,0.5,0,0.5);
			viewer->displayTrajectory();
			//viewer->setCameraPointing(Tnow.translation()(0),Tnow.translation()(1),Tnow.translation()(2)+3);
			viewer->repaint();
			viewer->win3D->process_events();
#endif
		}
		ctr++;
		}
	}
	}
	else
	{
		t3 = getDoubleTime();
		Tnow = Tnow * Tmotion;
	}

	t6 = getDoubleTime();
	if(fAddTimes!=NULL) 
	{
	    fprintf(fAddTimes,"%lf %lf %lf\n",t3-t2,t5-t4,t6-t0);
	    fflush(fAddTimes);
	}

	return Tnow;
}
Eigen::Affine3d NDTFuserHMT_SE::match(Eigen::Affine3d Tmotion, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,initializer_list<vector<double> > attributes)
{
	Todom = Todom * Tmotion; //we track this only for display purposes!
	lslgeneric::transformPointCloudInPlace(sensor_pose, *cloud);
	laserCloud_c.clear();
	laserCloud_c=getSegments(cloud,attributes,tails,ignore,removeProbability);
	for(int i=0;i<resolutions.size();i++)
		loadMap(mapLocal[i],laserCloud_c,NumInputs);
	Eigen::Affine3d Tinit = Tnow * Tmotion;
	bool matcher_res=true;
	if(!firstRun)
		for(auto i:resolutions_order)
		{
			matcher.current_resolution=resolutions.at(i);
			matcher_res&=matcher.match(map[i],mapLocal[i],Tinit,true);
		}
	if(matcher_res || fuseIncomplete)
	{
		Eigen::Affine3d diff = (Tnow * Tmotion).inverse() * Tinit;

		if((diff.translation().norm() > max_translation_norm || 
				diff.rotation().eulerAngles(0,1,2).norm() > max_rotation_norm) && checkConsistency)
		{
			fprintf(stderr,"****  NDTFuserHMT_L -- ALMOST DEFINATELY A REGISTRATION FAILURE *****\n");
			Tnow = Tnow * Tmotion;
		}
		else
		{
			Tnow = Tinit;
			for(int j=0;j<NumInputs;j++)
				lslgeneric::transformPointCloudInPlace(Tnow, *laserCloud_c[j]);
			Eigen::Affine3d diff_fuse = Tlast_fuse.inverse()*Tnow;
			if(diff_fuse.translation().norm() > translation_fuse_delta ||
				diff_fuse.rotation().eulerAngles(0,1,2).norm() > rotation_fuse_delta ||firstRun)
			{
				firstRun=false;
				canUpdate=true;
			}
		}
	}
	else
		Tnow = Tnow * Tmotion;
	return Tnow;
}
Eigen::Matrix<double, 6,6> NDTFuserHMT_SE::getPoseCovariance(Eigen::Affine3d T)
{
	Eigen::MatrixXd Covariance(6,6);
	Eigen::Matrix<double,6,6> Covariance_;
	matcher.covariance(map,mapLocal,T,resolutions,Covariance);
	Covariance_=Covariance;
	return Covariance_;
}
bool NDTFuserHMT_SE::updateMap()
{
	if(!canUpdate)return false;
	Eigen::Affine3d spose = Tnow*sensor_pose;
	for(unsigned int i=0;i<resolutions.size();i++)
		for(int j=0;j<NumInputs;j++)
		{
			map[i][j]->addPointCloudMeanUpdate(spose.translation(),*laserCloud_c[j],localMapSize, 1e5, 1250, *(size.begin()+2)/2, 0.06);
//		map[i][j]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE, 1e5, 1250, spose.translation(), 0.06);
			map[i][j]->computeNDTCells();
		}
	Tlast_fuse = Tnow;
	ctr++;
	canUpdate=false;
	return true;
}
}

