#ifndef NDT_VIZ_HH
#define NDT_VIZ_HH

#include <ndt_map/ndt_map.h>
#include <ndt_visualisation/NDTVizGlut.hh>
#include <pthread.h>

//#warning "ALWAYS PLACE THE ndt_viz.h BEFORE THE ROS HEADERS!!!"


class NDTViz {

    public:
	NDTVizGlut *win3D;
	NDTVizGlutPointCloudColor gl_points;
	NDTVizGlutPointCloudColor gl_particles;
	NDTVizGlutPointCloudColor gl_pointcloud;
	NDTVizGlutSetOfLines gl_laserlines;
	NDTVizGlutEllipsoids gl_ellipsoids;
        NDTVizGlutSetOfLines gl_grid;
        NDTVizGlutSetOfLines gl_polygonboundary;
        NDTVizGlutSetOfLines gl_pathline;

	NDTViz(bool allocate_new_window=true)

	{
	    if(allocate_new_window) 
	    {
	      win3D = new NDTVizGlut();
	      int argc=0;
	      char** argv = NULL;
	      win3D->win_run(&argc,argv);
	    }
	    else 
	    {
		win3D = NULL;
	    }
	    
	    gl_points.setPointSize(4.5);
	    gl_particles.setPointSize(2.5);

	}
	virtual ~NDTViz () {
	    if(win3D!=NULL) {
		win3D->win_close();
		delete win3D;
	    }
	}

	void startEventLoop() {
	  	
	}

	void repaint(){
	    win3D->repaint();
	}

	void clear(){
	  win3D->clearScene();
	}

	void clearTrajectoryPoints(){
	  gl_points.clear();
	}

	void addTrajectoryPoint(float x, float y, float z, float R=1.0, float G = 1.0, float B = 1.0){
	  gl_points.push_back(x, y, z, R ,G,B);
	}
	void displayTrajectory(){
	  win3D->addObject(&gl_points);
	}

	void clearParticles(){ gl_particles.clear();}
	void addParticle(float x, float y, float z, float R=1.0, float G = 1.0, float B = 1.0){
	    gl_particles.push_back(x, y, z, R ,G,B);
	}
	void displayParticles(){
	  win3D->addObject(&gl_particles);
	}
	void setCameraPointingToPoint(double x, double y, double z) {
	  win3D->setCameraPointingToPoint(x,y,z);
	}
	void setCameraPointing(double x, double y, double z) {
	  win3D->setCameraPointingToPoint(x,y,z);
	}

	
  void addGrid(size_t nbCells, double cellSize, double lineWidth = 1., double R = 1.0, double G = 0.4, double B = 0., double A = 0.8) {
          double length = nbCells*cellSize;
          double offset = length * 0.5;
          for (size_t i = 0; i <= nbCells; i++) {
            gl_grid.appendLine(i*cellSize-offset, -offset, 0, i*cellSize-offset, nbCells*cellSize-offset, 0);
            
            for (size_t i = 0; i <= nbCells; i++) {
              gl_grid.appendLine(-offset, i*cellSize-offset, 0, nbCells*cellSize-offset, i*cellSize-offset, 0);
            }
          }
          gl_grid.setColor4(R,G,B,A);
          gl_grid.setThickness(lineWidth);
          win3D->addObject(&gl_grid);
        }
  

	/**
	  * Add the laser scan to the scen 
	  */
	void addScan(Eigen::Vector3d orig, pcl::PointCloud<pcl::PointXYZ> &cloud, double R=1.0,double G=1.0,double B=1.0){
	  
	  gl_laserlines.clear();
	  for(unsigned int i=0;i<cloud.points.size();i+=2){
	    gl_laserlines.appendLine(orig(0),orig(1),orig(2), cloud.points[i].x, cloud.points[i].y, cloud.points[i].z);
	  }
	  gl_laserlines.setColor(R,G,B);
	  win3D->addObject(&gl_laserlines);
	}

        void addPolygon(std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d> > &pts, double R=1.0, double G=1.0, double B=1.0)  {
          gl_polygonboundary.clear();
          if (pts.size() < 2)
            return;
          for (unsigned int i = 0; i < pts.size(); i++) {
            unsigned int j = i+1;
            if (i >= pts.size() - 1) {
              j = 0; // final one - connect to first
            }
            gl_polygonboundary.appendLine(pts[i](0),pts[i](1),0.,
                                          pts[j](0),pts[j](1),0.);
          }
          gl_polygonboundary.setColor(R,G,B);
          gl_polygonboundary.setThickness(10.);
          win3D->addObject(&gl_polygonboundary);
        }

        void addPathLine(std::vector<Eigen::Vector2d, Eigen::aligned_allocator<Eigen::Vector2d> > &pts, double R=1.0, double G=1.0, double B=1.0)  {
          gl_pathline.clear();
          if (pts.size() < 2)
            return;
          for (unsigned int i = 0; i < pts.size()-1; i++) {
            unsigned int j = i+1;
            gl_pathline.appendLine(pts[i](0),pts[i](1),0.,
                                   pts[j](0),pts[j](1),0.);
          }
          gl_pathline.setColor(R,G,B);
          gl_pathline.setThickness(5.);
          win3D->addObject(&gl_pathline);
        }



	void addPointCloud(pcl::PointCloud<pcl::PointXYZ> &cloud, double R=1.0,double G=1.0,double B=1.0){

	  gl_pointcloud.setPointSize(2.5);
	  for(unsigned int i=0;i<cloud.points.size();i+=2){
	    gl_pointcloud.push_back(cloud.points[i].x, cloud.points[i].y, cloud.points[i].z,R,G,B);
	  }
	  win3D->addObject(&gl_pointcloud);
	}

	void plotNDTMap(perception_oru::NDTMap *map, double R=1.0,double G=1.0,double B=1.0, bool heightCoding=false, bool setCameraPos = true ){
	    if(win3D == NULL) return;
	    std::vector<perception_oru::NDTCell*> global_ndts;
	    global_ndts = map->getAllCells();
	    
	    win3D->clearScene();
	    gl_ellipsoids.clear();
	    
	    unsigned int accepted_ndts=1;
	    double x = 0,y=0,s=0;
	    //	    fprintf(stderr,"-NDT:%u-",global_ndts.size());
	    for(unsigned int i=0;i<global_ndts.size();i++){
		Eigen::Vector3d m = global_ndts[i]->getMean();
		if(!global_ndts[i]->hasGaussian_) continue;
		x+=m[0];
		y+=m[1];
		s+=1;

		accepted_ndts++;
		NDTVizGlutEllipsoid objEllip;
		Eigen::Matrix3d cov = global_ndts[i]->getCov();
		objEllip.setCovMatrix(cov);
		objEllip.setLocation(m[0], m[1], m[2]);
		double hval=1.0;
		if(heightCoding){
		    hval = fabs(m[2]+1.5)/6.0;
		}
		if(global_ndts[i]->getOccupancy() > 0){
		    objEllip.setColor(R/hval,G/hval,B/hval,1.0);
		}else{
		    objEllip.setColor(1.0,0,0,1.0);
		}

		objEllip.enableDrawSolid3D(true);
		gl_ellipsoids.push_back(objEllip);
	    }
	    win3D->addObject(&gl_ellipsoids);
	    if(setCameraPos) win3D->setCameraPointingToPoint(x/s,y/s,3.0);
	    win3D->repaint();
	    //fprintf(stderr,"(%lf %lf) s=%lf\n",x/s,y/s,s);
	    for(unsigned int i=0;i<global_ndts.size();i++) delete global_ndts[i];

	}

	void plotNDTSAccordingToCost(float occupancy, double MAX_COST, perception_oru::NDTMap *map){

	    if(win3D == NULL) return;
	    std::vector<perception_oru::NDTCell*> global_ndts;
	    global_ndts = map->getAllCells();
	    //	    fprintf(stderr," NUM NDT: %u ", global_ndts.size());

	    win3D->clearScene();
	    gl_ellipsoids.clear();
	    unsigned int accepted_ndts=1;
	    double x = 0,y=0,s=0;
	    for(unsigned int i=0;i<global_ndts.size();i++){
		Eigen::Vector3d m = global_ndts[i]->getMean();
		if(!global_ndts[i]->hasGaussian_) continue;
		x+=m[0];
		y+=m[1];
		s+=1;
		if(global_ndts[i]->getOccupancy()>occupancy){
		    accepted_ndts++;
		    NDTVizGlutEllipsoid objEllip;
		    Eigen::Matrix3d cov = global_ndts[i]->getCov();
		    objEllip.setCovMatrix(cov);
		    objEllip.setLocation(m[0], m[1], m[2]);
		    double c = global_ndts[i]->cost;
		    if(c > MAX_COST) c=MAX_COST;
		    c = 1 - c/MAX_COST; 
		    objEllip.setColor(c,c,c,0.6);

		    objEllip.enableDrawSolid3D(true);
		    gl_ellipsoids.push_back( objEllip );
		}else if(global_ndts[i]->getOccupancy()<-0){
		}
	    }
	    win3D->addObject(&gl_ellipsoids);
	    //win3D->setCameraPointingToPoint(x/s,y/s,3.0);
	    win3D->repaint();
	    fprintf(stderr,"(%lf %lf) s=%lf\n",x/s,y/s,s);
	    for(unsigned int i=0;i<global_ndts.size();i++) delete global_ndts[i];

	}
	/* /\** plots ndts according to the cell class, with an occupancy cutoff *\/ */
	void plotNDTSAccordingToClass(float occupancy, perception_oru::NDTMap *map){
	  
	  std::cout << "PLOT CLASS" << std::endl;
	    if(win3D == NULL) return;
	    std::vector<perception_oru::NDTCell*> global_ndts;
	    global_ndts = map->getAllCells();
	    //	    fprintf(stderr," NUM NDT: %u ", global_ndts.size());

	    Eigen::Vector3d cFlat, cInclined, cRough, cVert, cUnknown, color;
	    cFlat<<0,0.9,0;
	    cInclined<<0.1,0.2,0.5;
	    cRough<<0.9,0,0;
	    cVert<<0,0,0.9;
	    cUnknown<<0,0,0;

	    win3D->clearScene();
	    gl_ellipsoids.clear();

	    unsigned int accepted_ndts=1;
	    double x = 0,y=0,s=0;
	    for(unsigned int i=0;i<global_ndts.size();i++){
		Eigen::Vector3d m = global_ndts[i]->getMean();
		if(!global_ndts[i]->hasGaussian_) continue;
		x+=m[0];
		y+=m[1];
		s+=1;
		if(global_ndts[i]->getOccupancy()>occupancy){
		    accepted_ndts++;
		    NDTVizGlutEllipsoid objEllip;
		    Eigen::Matrix3d cov = global_ndts[i]->getCov();
		    objEllip.setCovMatrix(cov);
		    objEllip.setLocation(m[0], m[1], m[2]);
		    switch(global_ndts[i]->getClass()) {
			case perception_oru::NDTCell::HORIZONTAL :
			    color = cFlat;
			    break;
			case perception_oru::NDTCell::VERTICAL :
			    color = cVert;
			    break;
			case perception_oru::NDTCell::INCLINED :
			    color = cInclined;
			    break;
			case perception_oru::NDTCell::ROUGH :
			    color = cRough;
			    break;
			default:
			    color = cUnknown;
		    }
		    objEllip.setColor(color(0),color(1),color(2),0.6);

		    objEllip.enableDrawSolid3D(true);
		    gl_ellipsoids.push_back(objEllip);

		}else if(global_ndts[i]->getOccupancy()<-0){

		}
	    }
	    //win3D->setCameraPointingToPoint(x/s,y/s,3.0);
	    win3D->addObject(&gl_ellipsoids);
	    win3D->repaint();
	    fprintf(stderr,"(%lf %lf) s=%lf\n",x/s,y/s,s);
	    for(unsigned int i=0;i<global_ndts.size();i++) delete global_ndts[i];

	}

	void plotNDTSAccordingToOccupancy(float occupancy, perception_oru::NDTMap ***maps,std::vector<int> res, std::vector<int> sem){
	  if(win3D == NULL) return;
	    win3D->clearScene();
	    gl_ellipsoids.clear();
        for(auto _indexX:res)
		{
            for(auto _indexY:sem)
			{
				std::vector<perception_oru::NDTCell*> tempMap=(maps[_indexX][_indexY]->getAllCells());
				for(unsigned int i=0;i<tempMap.size();i++)
				{
					Eigen::Vector3d m = tempMap[i]->getMean();
					if(!tempMap[i]->hasGaussian_) continue;
					if(tempMap[i]->getOccupancy()>occupancy)
					{
						NDTVizGlutEllipsoid objEllip;
						Eigen::Matrix3d cov = tempMap[i]->getCov();
						objEllip.setCovMatrix(cov);
						objEllip.setLocation(m[0], m[1], m[2]);
						
			//		    objEllip.setColor((m[2]+2.0)/10.0,0,0,0.6);

						int R,G,B;
						R=(_indexY+1)&1;
						G=(_indexY+1)&2;
						B=(_indexY+1)&4;

						objEllip.setColor(R,G,B,0.6);

						objEllip.enableDrawSolid3D(true);
						gl_ellipsoids.push_back(objEllip);
					}
					else if(tempMap[i]->getOccupancy()<-0){std::cout<<"beep"<<std::endl; }
				}
				win3D->addObject(&gl_ellipsoids);
				win3D->repaint();
				for(unsigned int i=0;i<tempMap.size();i++) delete tempMap[i];
			}
		}

	}

	void plotLocalNDTMap(pcl::PointCloud<pcl::PointXYZ> &cloud, double resolution, double R=0, double G=1, double B=0, bool heightCoding=true){
	    if(win3D == NULL) return;

	    perception_oru::NDTMap ndlocal(new perception_oru::LazyGrid(resolution));
	    ndlocal.addPointCloudSimple(cloud);
	    ndlocal.computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE);

	    std::vector<perception_oru::NDTCell*> ndts;
	    ndts = ndlocal.getAllCells();
	    std::cout<<"LOCAL: "<<ndts.size()<<std::endl;

	    for(unsigned int i=0;i<ndts.size();i++){
		Eigen::Vector3d m = ndts[i]->getMean();
		//if(m[2]>3.0) continue;

		NDTVizGlutEllipsoid objEllip;
		Eigen::Matrix3d cov = ndts[i]->getCov();
		objEllip.setCovMatrix(cov);
		objEllip.setLocation(m[0], m[1], m[2]);
		double hval=1.0;
		if(heightCoding){
		    hval = fabs(m[2]+1.5)/20.0;
		}
		objEllip.setColor(R/hval,G/hval,B/hval,0.4);
		objEllip.enableDrawSolid3D(true);
		gl_ellipsoids.push_back(objEllip);
	    }
	    win3D->addObject(&gl_ellipsoids);
	    win3D->repaint();
	    for(unsigned int i=0;i<ndts.size();i++){
		delete ndts[i];
	    }

	}

	void plotLocalConflictNDTMap(perception_oru::NDTMap *map, pcl::PointCloud<pcl::PointXYZ> &cloud,
		double resolution, double R=1, double G=0, double B=0, bool heightCoding=false, double maxz=0){
	    if(win3D == NULL) return;

	    perception_oru::NDTMap ndlocal(new perception_oru::LazyGrid(resolution));

	    ndlocal.addPointCloudSimple(cloud);
	    ndlocal.computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE);

	    std::vector<perception_oru::NDTCell*> ndts;
	    ndts = ndlocal.getAllCells();
	    gl_ellipsoids.clear();
	    win3D->clearScene();

	    std::cout<<"LOCAL: "<<ndts.size()<<std::endl;

	    for(unsigned int i=0;i<ndts.size();i++){


		Eigen::Vector3d m = ndts[i]->getMean();
		if(m[2]>maxz) continue;
		pcl::PointXYZ p;
		p.x = m[0]; p.y=m[1]; p.z = m[2];
		perception_oru::NDTCell *map_cell=NULL;
		map->getCellAtPoint(p, map_cell);
		if(map_cell == NULL) continue;

		if(map_cell->getOccupancy()>0.5) continue;

		//	if(m[2]>3.0) continue;
		NDTVizGlutEllipsoid objEllip;
		Eigen::Matrix3d cov = ndts[i]->getCov();
		objEllip.setCovMatrix(cov);
		objEllip.setLocation(m[0], m[1], m[2]);
		double hval=1.0;
		if(heightCoding){
		    hval = fabs(m[2]+1.5)/20.0;
		}
		objEllip.setColor(R/hval,G/hval,B/hval,0.6);
		objEllip.enableDrawSolid3D(true);
		gl_ellipsoids.push_back(objEllip);
	    }
	    win3D->addObject(&gl_ellipsoids);
	    win3D->repaint();
	    for(unsigned int i=0;i<ndts.size();i++){
		delete ndts[i];
	    }

	}

	void plotNDTTraversabilityMap(perception_oru::NDTMap *map){
	    if(win3D == NULL) return;
	    std::vector<perception_oru::NDTCell*> global_ndts;
	    global_ndts = map->getAllCells();

	    gl_ellipsoids.clear();
	    win3D->clearScene();

	    for(unsigned int i=0;i<global_ndts.size();i++){
		Eigen::Vector3d m = global_ndts[i]->getMean();
		if(!global_ndts[i]->hasGaussian_) continue;

		NDTVizGlutEllipsoid objEllip;
		Eigen::Matrix3d cov = global_ndts[i]->getCov();
		objEllip.setCovMatrix(cov);
		objEllip.setLocation(m[0], m[1], m[2]);

		//perception_oru::CellClass CC;
		//perception_oru::NDTCell<PointT>::CellClass CC;

		//CC = global_ndts[i]->getClass();
		// {HORIZONTAL=0, VERTICAL, INCLINED, ROUGH, UNKNOWN};
		if(global_ndts[i]->getClass() == perception_oru::NDTCell::HORIZONTAL){
		    objEllip.setColor(0,1.0,0,1.0);
		}else if(global_ndts[i]->getClass() == perception_oru::NDTCell::VERTICAL){
		    objEllip.setColor(1.0,0,0,1.0);
		}else if(global_ndts[i]->getClass() == perception_oru::NDTCell::INCLINED){
		    objEllip.setColor(1.0,1.0,0,1.0);
		}else if(global_ndts[i]->getClass() == perception_oru::NDTCell::ROUGH){
		    objEllip.setColor(0,0,1.0,1.0);
		}else{
		    objEllip.setColor(1.0,1.0,1.0,1.0);
		}

		objEllip.enableDrawSolid3D(true);
		gl_ellipsoids.push_back(objEllip);
	    }
	    win3D->addObject(&gl_ellipsoids);
	    win3D->repaint();
	    for(unsigned int i=0;i<global_ndts.size();i++) delete global_ndts[i];

	}
	/* ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// */
	/* ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// */
	/* ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// */
	/* ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////// */
	/* /\** */
	/*  *  Computes and visualizes an NDT depth image, based on maximum likelihood estimation */
	/*  *  */
	/*  *\/ */
	/* void ndtCoolCam(perception_oru::NDTMap *map ,const Eigen::Affine3d &spos, double maxDist=70.0,  */
	/* 	unsigned int Nx=800, unsigned int Ny=600, double fx=800, double fy=600){ */
	/*     Eigen::Matrix3d K; */
	/*     K << fx,0,(double)Nx/2.0, */
	/*       0, fy, (double)Ny/2.0, */
	/*       0,0,1; */

	/*     Eigen::Matrix3d invK; */
	/*     double det=0; */
	/*     bool exists=false; */
	/*     K.computeInverseAndDetWithCheck(invK,det,exists); */

	/*     if(!exists){ */
	/* 	fprintf(stderr,"ndtCoolCam::ERROR: NO INVERSE!!!!\n"); */
	/* 	return; */
	/*     } */

	/*     Eigen::Matrix3d m1,m2; ///Orientation matrix, for rotating the camera to world frame */
	/*     m1 = Eigen::AngleAxisd(-M_PI/2.0, Eigen::Vector3d::UnitX()) */
	/* 	* Eigen::AngleAxisd(0, Eigen::Vector3d::UnitY()) */
	/* 	* Eigen::AngleAxisd(0, Eigen::Vector3d::UnitZ()); */

	/*     m2 = Eigen::AngleAxisd(0, Eigen::Vector3d::UnitX()) */
	/* 	* Eigen::AngleAxisd(0, Eigen::Vector3d::UnitY()) */
	/* 	* Eigen::AngleAxisd(-M_PI/2.0, Eigen::Vector3d::UnitZ()); */

	/*     mrpt::opengl::CTexturedPlanePtr gl_img 		=  mrpt::opengl::CTexturedPlane::Create(0.5,-0.5,-0.5,0.5); */
	/*     mrpt::opengl::COpenGLViewportPtr viewInt; */
	/*     //mrpt::utils::CImage  img(Nx,Ny,CH_GRAY); */
	/*     mrpt::utils::CImage  img(Nx,Ny,CH_RGB); */


	/*     //mrpt::opengl::COpenGLScenePtr &scene = win3D->get3DSceneAndLock();	 */

	/*     //mrpt::opengl::CSetOfLinesPtr obj = mrpt::opengl::CSetOfLines::Create();		 */
	/*     for(unsigned int i=0; i<Nx*Ny;i++){ */
	/* 	Eigen::Vector3d pix(i%Nx,(int)i/Nx,1); */
	/* 	Eigen::Vector3d dirv = invK*pix; ///< direction vector based on camera model */
	/* 	dirv.normalize(); */

	/* 	Eigen::Vector3d dirv_w = spos.rotation()*(m2*(m1*dirv)); ///Direction vector in world frame */
	/* 	//double max_depth = 20; */
	/* 	double depth = map->getDepth(spos.translation(), dirv_w, maxDist); */

	/* 	if(depth>maxDist) img.setPixel(i%Nx,i/Nx,0); */
	/* 	else{ */
	/* 	    float x = (depth/maxDist); */
	/* 	    float r = ((x >= 3.0/8.0) & (x < 5.0/8.0))*(4. * x - 3./2.)+((x >= 5./8.) & (x < 7./8.))+(x >= 7./8.) * (-4. * x + 9./2.); */
	/* 	    float g = ((x >= 1./8.) & (x < 3./8.))*(4. * x - 1./2.)+((x >= 3./8.) & (x < 5./8.))+((x >= 5./8.) & (x < 7./8.))*(-4. * x + 7./2.); */
	/* 	    float b = (x < 1./8.)*(4. * x + 1./2.)+((x >= 1./8.) & (x < 3./8.))+((x >= 3./8.) & (x < 5./8.))*(-4. * x + 5./2.); */
	/* 	    size_t color=((unsigned char)(r*255))*255*255 + ((unsigned char)(g*255))*255 + ((unsigned char)(b*255)); */

	/* 	    img.setPixel(i%Nx,i/Nx,color); */
	/* 	    //fprintf(stderr,"(%.2lf) ",depth); */
	/* 	} */
	/* 	/\* */
	/* 	   Eigen::Vector3d ray_endpos=spos.translation() + dirv_w * max_depth; */
	/* 	   obj->appendLine(spos.translation()(0),spos.translation()(1),spos.translation()(2), */
	/* 	   ray_endpos(0), ray_endpos(1), ray_endpos(2));*\/ */

	/*     } */
	/*     //scene->insert(obj); */

	/*     mrpt::opengl::COpenGLScenePtr &scene = win3D->get3DSceneAndLock();	 */
	/*     scene->clear(); */

	/*     viewInt = scene->createViewport("view2d"); */
	/*     scene->insert( gl_img, "view2d"); */

	/*     viewInt->setViewportPosition(0, 0, -1.0, -1.0); */

	/*     viewInt->setTransparent(false); */

	/*     viewInt->getCamera().setOrthogonal(true); */
	/*     viewInt->getCamera().setAzimuthDegrees(90); */
	/*     viewInt->getCamera().setElevationDegrees(90); */
	/*     viewInt->getCamera().setZoomDistance(1.0); */

	/*     gl_img->assignImage_fast(img); */

	/*     win3D->unlockAccess3DScene(); */
	/*     win3D->repaint(); */
	/* }	 */
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



};

#endif
