#include <set>
#include <cstdlib>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include "sstream"
#include "string"
#include "iostream"
#include "boost/serialization/serialization.hpp"
#include "boost/serialization/base_object.hpp"
#include "boost/shared_ptr.hpp"
#include "boost/serialization/vector.hpp"

#include <Eigen/Eigen>
#include <pcl/point_cloud.h>
//#include <pcl/features/feature.h>
#include <pcl/common/common.h>

#include <string>
#include <climits>

//#include <ndt_map/oc_tree.h>
#include <ndt_map/ndt_map.h>

#include <cstring>
#include <cstdio>

namespace perception_oru
{
inline bool isnan(const double* const pt)
{
     return std::isnan(pt[0]) ||std::isnan(pt[1]) ||std::isnan(pt[2]);
}
/**
* loadPointCloud - You can call this if you are only interested in dealing with one scan
* without need for fusing several ones or representing empty space and occupancy
*
* Otherwise you should always call addPointCloud (or if you don't want occupancy then addPointCloudSimple)
*
* \param pc the PointCloud that is to be loaded
* \note every subsequent call will destroy the previous map!
*/
void NDTMap::loadPointCloud(const std::vector<const double* const> &pc, double range_limit)
{
    if(index_ != NULL)
    {
        //std::cout<<"CLONE INDEX\n";
        SpatialIndex *si = index_->clone();
        //cout<<"allocating index\n";
        if(!isFirstLoad_)
        {
            //std::cout<<"deleting old index\n";
            delete index_;
        }
        isFirstLoad_ = false;
        index_ = si;
    }
    else
    {
        //NULL index in constructor, abort!
        //ERR("constructor must specify a non-NULL spatial index\n");
        return;
    }
    
    LazyGrid *lz = dynamic_cast<LazyGrid*>(index_);
    if(lz == NULL)
    {
	fprintf(stderr,"Unfortunately This works only with Lazygrid!\n");
	exit(1);
    }

    if(index_ == NULL)
    {
        //ERR("Problem creating index, unimplemented method\n");
        return;
    }

    std::vector<const double* const>::const_iterator it = pc.begin();
    if(guess_size_) 
    {
	double maxDist = 0;//, distCeil = 200;

	Eigen::Vector3d centroid(0,0,0);
	int npts = 0;
	while(it!=pc.end())
	{
	    if(isnan(*it))
        {
            it++;
            continue;
        }
        Eigen::Map<const Eigen::Vector3d> d(*it);
	    if(range_limit>0)
            if(d.norm()>range_limit){
                it++;continue;}
	    centroid += d;
	    it++;
	    npts++;
	}

	centroid /= (double)npts;
	double maxz=-1000, minz=10000;
	//Eigen::Vector4f centroid(0,0,0,0);
	//pcl::compute3DCentroid(pc,centroid);

	//compute distance to furthest point
	it = pc.begin();
	while(it!=pc.end())
	{
	    if(isnan(*it))
        {
            it++;
            continue;
        }
        Eigen::Map<const Eigen::Vector3d> d(*it);
	    if(range_limit>0)
            if(d.norm()>range_limit){
                it++;
                continue;
            }
        Eigen::Vector3d dist_v=centroid-d;
	    double dist = dist_v.norm();
	    maxDist = (dist > maxDist) ? dist : maxDist;
	    maxz = (dist_v(2) > maxz) ? dist_v(2) : maxz;
	    minz = (dist_v(2) < minz) ? dist_v(2) : minz;
	    it++;
	}
	// cout<<"Points = " <<pc.points.size()<<" maxDist = "<<maxDist<<endl;
	NDTCell *ptCell = new NDTCell();
	index_->setCellType(ptCell);
	delete ptCell;
	index_->setCenter(centroid.data());

	if(map_sizex >0 && map_sizey >0 && map_sizez >0)
	{
	    index_->setSize(map_sizex,map_sizey,map_sizez);
	}
	else
	{
	    index_->setSize(4*maxDist,4*maxDist,3*(maxz-minz));
	}
    }
    else
    {
	//set sizes
	NDTCell *ptCell = new NDTCell();
	index_->setCellType(ptCell);
	delete ptCell;
	index_->setCenter(centerx,centery,centerz);
	if(map_sizex >0 && map_sizey >0 && map_sizez >0)
	{
	    index_->setSize(map_sizex,map_sizey,map_sizez);
	}
    }

    //    ROS_INFO("centroid is %f,%f,%f", centroid(0),centroid(1),centroid(2));
    //    ROS_INFO("maxDist is %lf", maxDist);

    it = pc.begin();
    while(it!=pc.end())
    {
        Eigen::Map<const Eigen::Vector3d> d(*it);
        if(isnan(*it))
        {
            it++;
            continue;
        }
        if(range_limit>0)
        {
            if(d.norm()>range_limit)
            {
                it++;
                continue;
            }
        }
        index_->addPoint(*it);
	NDTCell *ptCell;
	lz->getNDTCellAt(*it,ptCell);
#ifdef REFACTORED
	if(ptCell!=NULL) {
	//    std::cerr<<"adding to update set\n";
	    update_set.insert(ptCell);
	}
#endif
        it++;
    }

    isFirstLoad_ = false;
}

/**
* loadPointCloudCentroid - A special load function to enable the matching of centroids (create alligned maps)
* \param pc the PointCloud that is to be loaded
* \note every subsequent call will destroy the previous map!
*/
void NDTMap::loadPointCloudCentroid(const std::vector<const double* const> &pc, const Eigen::Vector3d &origin,
	const Eigen::Vector3d &old_centroid,const Eigen::Vector3d &map_size, double range_limit){
    
    if(index_ != NULL){
        SpatialIndex *si = index_->clone();
        if(!isFirstLoad_) delete index_;
        
        isFirstLoad_ = false;
        index_ = si;
    }
    else{
	return;
    }

    if(index_ == NULL) return; ///Should never happen!

    NDTCell *ptCell = new NDTCell();
    index_->setCellType(ptCell);
    delete ptCell;

    LazyGrid *lz = dynamic_cast<LazyGrid*>(index_);
    if(lz == NULL)
    {
	fprintf(stderr,"Unfortunately This works only with Lazygrid!\n");
	exit(1);
    }

    Eigen::Vector3d diff = origin - old_centroid;  //-origin;
    double cx=0, cy=0, cz=0;
    lz->getCellSize(cx, cy, cz);

    ///How many cell to each direction is the new origin from old one
    Eigen::Vector3d centroid(cx,cy,cz);
    centroid=old_centroid+floor(diff/centroid)*centroid;

    index_->setCenter(centroid.data());
    //index_->setCenter(origin(0),origin(1),origin(2));
    index_->setSize(map_size(0),map_size(1),map_size(2));
    //lz->initializeAll();

   //temporarily commented fprintf(stderr,"centroid is %lf,%lf,%lf (origin: %lf %lf %lf) (map_size %lf %lf %lf) N=%d", centroid(0),centroid(1),centroid(2), origin(0),origin(1),origin(2), map_size(0), map_size(1), map_size(2),(int)pc.size());
    // ROS_INFO("centroid is %f,%f,%f", centroid(0),centroid(1),centroid(2));
    // ROS_INFO("maxDist is %lf", maxDist);

    std::vector<const double* const>::const_iterator it = pc.begin();
    while(it!=pc.end())
    {
        Eigen::Map<Eigen::Vector3d> d;
        if(isnan(*it))
        {
            it++;
            continue;
        }
        
            //        std::cout << "centoroid add point [" << it->x << "," << it->y << "," <<it->z <<std::endl;

        if(range_limit>0)
        {
            if((d-origin).norm()>range_limit)
            {
            it++;
            continue;
                }
        }
        
        //fprintf(stderr,"HEP!");
        index_->addPoint(*it);
        NDTCell *ptCell=NULL;
        lz->getNDTCellAt(*it,ptCell);
#ifdef REFACTORED
        if(ptCell!=NULL) {
              update_set.insert(ptCell);
              //          std::cout << "insert" << std::endl;
            }
            // else {
            //   std::cout << "invalid cell..." << *it << std::endl;
            //   int indX, indY, indZ;
            //   lz->getIndexForPoint(*it, indX, indY, indZ);
            //   std::cout << "ind : " << indX << "," << indY << "," << indZ << std::endl;
            // }
            {
              // int indX, indY, indZ;
              // lz->getIndexForPoint(*it, indX, indY, indZ);
              // //          std::cout << "ind : " << indX << "," << indY << "," << indZ << std::endl;
              // double dx,dy,dz;
              // lz->getCenter(dx,dy,dz);
              // //          std::cout << "center : " << dx << "," << dy << "," << dz << std::endl;
              // int indX, indY, indZ;
              // lz->getGridSize(indX, indY, indZ);
              // //          std::cout << "gridsize : " << indX << "," << indY << "," << indZ << std::endl;
            }
#endif
        it++;
    }

    isFirstLoad_ = false;
}



/**
 * Just adds the points, without raytracing and such
 */
void NDTMap::addPointCloudSimple(const std::vector<const double* const> &pc,double maxz)
{
    if(isFirstLoad_)
    {
        loadPointCloud( pc);
        return;
    }

    std::vector<const double* const>::const_iterator it = pc.begin();
    it = pc.begin();
    LazyGrid *lz = dynamic_cast<LazyGrid*>(index_);
    if(lz == NULL)
    {
	fprintf(stderr,"Unfortunately This works only with Lazygrid!\n");
	exit(1);
    }

    while(it!=pc.end())
    {
        if(isnan(*it))
        {
            it++;
            continue;
        }
        if(it->z>maxz)
        {
            it++;
            continue;
        }
        index_->addPoint(*it);
	NDTCell *ptCell;
	lz->getNDTCellAt(*it,ptCell);
#ifdef REFACTORED
	if(ptCell!=NULL) update_set.insert(ptCell);
#endif
        it++;
    }

}

/**
* Add a distribution to the map
*/
void NDTMap::addDistributionToCell(const Eigen::Matrix3d &ucov, const Eigen::Vector3d &umean, unsigned int numpointsindistribution, 
                                   float r, float g,float b,  unsigned int maxnumpoints, float max_occupancy)
{
    pcl::PointXYZ pt;
    pt.x = umean[0];
    pt.y = umean[1];
    pt.z = umean[2];
    LazyGrid *lz = dynamic_cast<LazyGrid*>(index_);
    if(lz==NULL)
    {
        fprintf(stderr,"NOT LAZY GRID!!!\n");
        exit(1);
    }
    NDTCell *ptCell = NULL; 
    //    lz->getNDTCellAt(pt,ptCell);
    lz->getCellAtAllocate(pt,ptCell);
    
    if(ptCell != NULL)
    {
	//std::cout<<"BEFORE\n";
	//std::cout<<"had G "<<ptCell->hasGaussian_<<" occ "<<ptCell->getOccupancy()<<std::endl;
//	std::cout<<umean.transpose()<<" "<<numpointsindistribution <<std::endl;
//	std::cout<<ucov<<std::endl;
//	std::cout<<ptCell->getMean().transpose()<<std::endl;
//	std::cout<<ptCell->getCov()<<std::endl;
	ptCell->updateSampleVariance(ucov, umean, numpointsindistribution, true, max_occupancy,maxnumpoints);
	ptCell->setRGB(r,g,b);
//	std::cout<<"AFTER\n";
//	std::cout<<ptCell->getMean().transpose()<<std::endl;
//	std::cout<<ptCell->getCov()<<std::endl;
    }
    else {
      //      std::cerr << "addDistributionToCell: failed to get a cell to add the distribution to" << std::endl;
    }
}

///Get the cell for which the point fall into (not the closest cell)
bool NDTMap::getCellAtPoint(const pcl::PointXYZ &refPoint, NDTCell *&cell)
{
    LazyGrid *lz = dynamic_cast<LazyGrid*>(index_);
    if(lz==NULL)
    {
        fprintf(stderr,"NOT LAZY GRID!!!\n");
        exit(1);
    }
    lz->getNDTCellAt(refPoint,cell);
    return (cell != NULL);
}

bool NDTMap::getCellAtPoint(const pcl::PointXYZ &refPoint, NDTCell *&cell) const 
{
    LazyGrid *lz = dynamic_cast<LazyGrid*>(index_);
    if(lz==NULL)
    {
        fprintf(stderr,"NOT LAZY GRID!!!\n");
        exit(1);
    }
    lz->getNDTCellAt(refPoint,cell);
    return (cell != NULL);
}

///Get the cell for which the point fall into (not the closest cell)
bool NDTMap::getCellAtAllocate(const pcl::PointXYZ &refPoint, NDTCell *&cell)
{
    LazyGrid *lz = dynamic_cast<LazyGrid*>(index_);
    if(lz==NULL)
    {
        fprintf(stderr,"NOT LAZY GRID!!!\n");
        exit(1);
    }
    lz->getCellAtAllocate(refPoint,cell);
    return (cell != NULL);
}

bool NDTMap::getCellAtAllocate(const pcl::PointXYZ &refPoint, NDTCell *&cell) const 
{
    LazyGrid *lz = dynamic_cast<LazyGrid*>(index_);
    if(lz==NULL)
    {
        fprintf(stderr,"NOT LAZY GRID!!!\n");
        exit(1);
    }
    lz->getCellAtAllocate(refPoint,cell);
    return (cell != NULL);
}

/**
 * Adds a new cloud: NDT-OM update step
 */
void NDTMap::addPointCloud(const Eigen::Vector3d &origin, const std::vector<const double* const> &pc, double classifierTh, double maxz, 
	double sensor_noise, double occupancy_limit)
{
// 	std::cout << "Good function: addPointCloud" << std::endl;
// 	exit(0);
    if(isFirstLoad_)
    {
			loadPointCloud( pc);
			return;
    }
    if(index_ == NULL)
    {
			//ERR("Problem creating index, unimplemented method\n");
			return;
    }
    std::vector<const double* const>::const_iterator it = pc.begin();

    LazyGrid *lz = dynamic_cast<LazyGrid*>(index_);
    if(lz==NULL)
    {
			fprintf(stderr,"NOT LAZY GRID!!!\n");
			exit(1);
    }
    //auto cA=lz->size();
    lz->setSensorPose(origin.data());
    pcl::PointXYZ pt;
    NDTCell* ptCell = NULL; 

#ifdef REFACTORED
    std::vector< NDTCell*> cells; 
    bool updatePositive = true;
    double max_range = 200.;

    while(it!=pc.end())
    {
            Eigen::Map<const Eigen::Vector3d> pt(*it);
			if(isnan(*it))
			{
					it++;
					continue;
			}

            Eigen::Vector3d diff=pt-origin;
			double l = diff.norm();

			if(l>max_range)
			{
          fprintf(stderr,"addPointCloud::Very long distance (%lf) :( \n",l);
					it++;
					continue;
			}

			cells.clear();
// 			std::cout << "Tracing the line" << std::endl;
			if(!lz->traceLine(origin,*it,diff,maxz,cells)) {
					it++;
					continue;
			}
// 			else{
// 				std::cout << "That is not true so it's good" << std::endl;
// 			}

			for(unsigned int i=0; i<cells.size(); i++)
			{
					ptCell = cells[i];
					if(ptCell != NULL)
					{
				double l2target = 0;
				if(ptCell->hasGaussian_)
                {
                    Eigen::Map<const Eigen::Vector3d> pend(*it);
                    Eigen::Vector3d out, vpt;
                    double lik = ptCell->computeMaximumLikelihoodAlongLine(origin.data(), *it, out);
                    l2target = (out-pend).norm();

                    double dist = (origin-out).norm();
                    if(dist > l) continue; ///< don't accept points further than the measurement

                    l2target = (out-pend).norm(); ///<distance to endpoint

                    double sigma_dist = 0.5 * (dist/30.0); ///test for distance based sensor noise
                    double snoise = sigma_dist + sensor_noise;
                    double thr =exp(-0.5*(l2target*l2target)/(snoise*snoise)); ///This is the probability of max lik point being endpoint
                    lik *= (1.0-thr);
                    if(lik<0.3) continue;
                    lik = 0.1*lik+0.5; ///Evidence value for empty - alpha * p(x);
                    double logoddlik = log( (1.0-lik)/(lik) );
                    //ptCell->updateEmpty(logoddlik,l2target);
                    //fprintf(stderr,"[E=%.2lf] ", logoddlik);
                    ptCell->updateOccupancy(logoddlik, occupancy_limit);
                    if(ptCell->getOccupancy()<=0) ptCell->hasGaussian_ = false; 
                }
				else
				{
					 // ptCell->updateEmpty(-0.2,l2target); ///The cell does not have gaussian, so we mark that we saw it empty...
					 ptCell->updateOccupancy(-0.2, occupancy_limit);
					 if(ptCell->getOccupancy()<=0) ptCell->hasGaussian_ = false; 
				}
					}
					//update_set.insert(ptCell);
			}
	if(updatePositive) {
	    ptCell = dynamic_cast<NDTCell*>(index_->addPoint(*it));
	    if(ptCell!=NULL) {
				update_set.insert(ptCell);
	    }
	}
	it++;
    }
    isFirstLoad_ = false;
    //auto cB=lz->size();
    //std::cerr<<"SIZEDIFF  "<<cA<<"\t"<<cB<<"\t"<<cB-cA<<std::endl;

#else
    double centerX,centerY,centerZ;
    lz->getCenter(centerX, centerY, centerZ);
    double cellSizeX,cellSizeY,cellSizeZ;
    lz->getCellSize(cellSizeX, cellSizeY, cellSizeZ);
    int sizeX,sizeY,sizeZ;
    lz->getGridSize(sizeX, sizeY, sizeZ);

    NDTCell ****dataArray = lz->getDataArrayPtr();

    double min1 = std::min(cellSizeX,cellSizeY);
    double min2 = std::min(cellSizeZ,cellSizeY);

    double resolution = std::min(min1,min2); ///Select the smallest resolution

    while(it!=pc.end())
    {
        if(isnan(*it))
        {
            it++;
            continue;
        }
        Eigen::Vector3d diff;
        diff << it->x-origin(0), it->y-origin(1), it->z-origin(2);

        double l = diff.norm();
        int N = l / (resolution);

        if(l>200)
        {
            fprintf(stderr,"addPointCloud::Very long distance (%lf) :( \n",l);
            it++;
            continue;
        }
        if(resolution<0.01)
        {
            fprintf(stderr,"Resolution very very small (%lf) :( \n",resolution);
            it++;
            continue;
        }


        if(N <= 0)
        {
            fprintf(stderr,"N=%d (r=%lf l=%lf) :( ",N,resolution,l);
            it++;
            continue;
        }

        diff = diff/(float)N;

        bool updatePositive = true;
        if(it->z>maxz)
        {
            it++;
            continue;
        }

        int idxo=0, idyo=0,idzo=0;
        for(int i=0; i<N-1; i++)
        {
            pt.x = origin(0) + ((float)(i+1)) *diff(0);
            pt.y = origin(1) + ((float)(i+1)) *diff(1);
            pt.z = origin(2) + ((float)(i+1)) *diff(2);
            int idx,idy,idz;

            idx = (int)(((pt.x - centerX)/cellSizeX+0.5) + sizeX/2);
            idy = (int)(((pt.y - centerY)/cellSizeY+0.5) + sizeY/2);
            idz = (int)(((pt.z - centerZ)/cellSizeZ+0.5) + sizeZ/2);


            ///We only want to check every cell once, so
            ///increase the index if we are still in the same cell
            if(idx == idxo && idy==idyo && idz ==idzo)
            {
                continue;
            }
            else
            {
                idxo = idx;
                idyo=idy,idzo=idz;
            }
            /// Check the validity of the index
            if(idx < sizeX && idy < sizeY && idz < sizeZ && idx >=0 && idy >=0 && idz >=0)
            {
                //fprintf(stderr,"(in)");
                ptCell = (dataArray[idx][idy][idz]);
            }
            else
            {
                //fprintf(stderr,"(out)");
                continue;
            }

            if(ptCell != NULL)
            {
                double l2target = 0;
                if(ptCell->hasGaussian_)
                {
                    Eigen::Vector3d out, pend,vpt;
                    pend << it->x,it->y,it->z;
                    double lik = ptCell->computeMaximumLikelihoodAlongLine(origin.data(), pt, out);
                    l2target = (out-pend).norm();

                    double dist = (origin-out).norm();
                    if(dist > l) continue; ///< don't accept points further than the measurement

                    l2target = (out-pend).norm(); ///<distance to endpoint

                    double sigma_dist = 0.5 * (dist/30.0); ///test for distance based sensor noise
                    double snoise = sigma_dist + sensor_noise;
                    double thr =exp(-0.5*(l2target*l2target)/(snoise*snoise)); ///This is the probability of max lik point being endpoint
                    lik *= (1.0-thr);
                    lik = 0.2*lik+0.5; ///Evidence value for empty - alpha * p(x);
                    double logoddlik = log( (1.0-lik)/(lik) );
                    ptCell->updateEmpty(logoddlik,l2target);
                }
                else
                {
                    ptCell->updateEmpty(-0.2,l2target); ///The cell does not have gaussian, so we mark that we saw it empty...
                }
            }
            else
            {
                index_->addPoint(pt); ///Add fake point to initialize!
            }
        }

        if(updatePositive) index_->addPoint(*it);
        it++;
    }
    isFirstLoad_ = false;
#endif
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/**
* Add new pointcloud to map - Updates the occupancy using the mean values of 
* a local map generated from an observation
* 
* Performs raytracing, updates conflicts and adds points to cells
* computeNDTCells must be called after calling this
*
* @param &origin is the position of the sensor, from where the scan has been taken from.
* @param &pc is the pointcloud to be added
* @param &localmapsize The dimensions of the local map used for computing the gaussians
* @param maxnumpoints Defines the forgetting factor (default 100000) the smaller the value the faster the adaptation
* @param occupancy_limit Clamping threshold for log-odds value
* @param maxz threshold for the maximum z-coordinate value for the measurement point_cloud
* @param sensor_noise The expected standard deviation of the sensor noise
*/
void  NDTMap::addPointCloudMeanUpdate(const Eigen::Vector3d &origin, 
	const std::vector<const double* const> &pc, 
	const Eigen::Vector3d &localmapsize,
	unsigned int maxnumpoints, float occupancy_limit,double maxz, double sensor_noise){

    if(isFirstLoad_){
	//fprintf(stderr,"First load... ");
	loadPointCloud(pc);
	computeNDTCells();
	//fprintf(stderr,"DONE!");
	return;
    }
    if(index_ == NULL){
	return;
    }
    LazyGrid *lz = dynamic_cast<LazyGrid*>(index_);
    if(lz==NULL){
	fprintf(stderr,"NOT LAZY GRID!!!\n");
	exit(1);
    }
    ///fprintf(stderr,"UPDATING\n");
    double centerX,centerY,centerZ;
    lz->getCenter(centerX, centerY, centerZ);
    double cellSizeX,cellSizeY,cellSizeZ;
    lz->getCellSize(cellSizeX, cellSizeY, cellSizeZ);
    int sizeX,sizeY,sizeZ;
    lz->getGridSize(sizeX, sizeY, sizeZ);

    Eigen::Vector3d old_centroid;
    old_centroid(0) = centerX;
    old_centroid(1) = centerY;
    old_centroid(2) = centerZ;

    double min1 = std::min(cellSizeX,cellSizeY);
    double min2 = std::min(cellSizeZ,cellSizeY);
    double resolution = std::min(min1,min2); ///Select the smallest resolution
		//fprintf(stderr,"1:");
    ///Lets first create a local ndmap (this really works only if we have lazy grid as well)
    perception_oru::NDTMap ndlocal(new perception_oru::LazyGrid(resolution));
    ndlocal.loadPointCloudCentroid(pc,origin, old_centroid, localmapsize, 70.0); ///FIXME:: fixed max length
    
    ///Use Student-T
    //ndlocal.computeNDTCells(CELL_UPDATE_MODE_STUDENT_T);
    
    /// Use Sample variance
    ndlocal.computeNDTCells();
		
    ///TODO: This returns now copies --- should be pointers?
    std::vector<perception_oru::NDTCell*> ndts;
    ndts = ndlocal.getAllCells();
		//fprintf(stderr,"2(%u):",ndts.size());
    NDTCell *ptCell=NULL;

    Eigen::Matrix<double,3,1> pt;
    //fprintf(stderr,"NUM = %d\n",ndts.size());
    int num_high = 0;
    for(unsigned int it=0;it<ndts.size();it++)
    {
	if(ndts[it] == NULL){
	    fprintf(stderr,"NDTMap::addPointCloudMeanUpdate::GOT NULL FROM MAP -- SHOULD NEVER HAPPEN!!!\n");
	    continue;
	}

	if(!ndts[it]->hasGaussian_){
	    fprintf(stderr,"NDTMap::addPointCloudMeanUpdate::NO GAUSSIAN!!!! -- SHOULD NEVER HAPPEN!!!\n");
	    continue;
	}

	int numpoints = ndts[it]->getN();

	if(numpoints<=0){
	    fprintf(stderr,"addPointCloudMeanUpdate::Number of points in distribution<=0!!");
	    continue;
	}
	Eigen::Vector3d diff;
	Eigen::Vector3d m = ndts[it]->getMean();
	diff = m-origin;

	double l = diff.norm();
	int NN = l / (resolution);

	if(l>200){
	    fprintf(stderr,"addPointCloudMeanUpdate::Very long distance (%lf) :( \n",l);
	    continue;
	}
	if(resolution<0.01){
	    fprintf(stderr,"addPointCloudMeanUpdate::Resolution very very small (%lf) :( \n",resolution);
	    continue;
	}
	if(NN < 0){ 
	    fprintf(stderr,"addPointCloudMeanUpdate::N=%d (r=%lf l=%lf) :( ",NN,resolution,l);
	    continue;
	}

	bool updatePositive = true;
	if(m(2)>maxz){
	    num_high++;
	    updatePositive = false; ///Lets update negative even though the measurement was too high
	}
	//std::cout<<"c: "<<ndts[it]->getCov()<<std::endl;

	diff = diff/(float)NN;

	//fprintf(stderr,"3(%d):",NN);
	int idxo=0, idyo=0,idzo=0;
	for(int i=0; i<NN-2; i++)  
	{
        pt=origin+ ((float)(i+1)) *diff;
	    int idx,idy,idz;

	    idx = (int)(((pt(0) - centerX)/cellSizeX+0.5) + sizeX/2.0);
	    idy = (int)(((pt(1) - centerY)/cellSizeY+0.5) + sizeY/2.0);
	    idz = (int)(((pt(2) - centerZ)/cellSizeZ+0.5) + sizeZ/2.0);


	    ///We only want to check every cell once, so
	    ///increase the index if we are still in the same cell
	    if(idx == idxo && idy==idyo && idz ==idzo){
		continue;
	    }else{
		idxo = idx;idyo=idy;idzo=idz;
	    }
	    ptCell = NULL;
	    /// Check the validity of the index
	    lz->getNDTCellAt(idx,idy,idz,ptCell);
	    /*
	       if(ptCell != NULL){
	       if(!ptCell->hasGaussian_){
	       ptCell->updateOccupancy(-0.85*numpoints,occupancy_limit);
	       }else{
	       ptCell->updateOccupancy(-0.2*numpoints,occupancy_limit);
	       ptCell = lz->getClosestNDTCell(pt);
	       }
	       }*/

	    if(ptCell != NULL){
		double l2target = 0;
		if(ptCell->hasGaussian_)
		{
		    Eigen::Vector3d out, vpt;

		    double lik = ptCell->computeMaximumLikelihoodAlongLine(origin.data(), pt.data(), out);
		    l2target = (out-m).norm();

		    double dist = (origin-out).norm();
		    if(dist > l) continue; ///< don't accept points further than the measurement

		    l2target = (out-m).norm(); ///<distance to endpoint

		    double sigma_dist = 0.5 * (dist/30.0); ///test for distance based sensor noise
		    //double sigma_dist = 0;
		    double snoise = sigma_dist + sensor_noise;
		    double thr =exp(-0.5*(l2target*l2target)/(snoise*snoise)); ///This is the probability of max lik point being endpoint
		    lik *= (1.0-thr);
		    //lik = 0.4*lik+0.5+0.05; ///Evidence value for empty - alpha * p(x);
		    lik = 0.2*lik+0.5; ///Evidence value for empty - alpha * p(x);
		    double logoddlik = log( (1.0-lik)/(lik) );

		    //fprintf(stderr,"l=%lf ",numpoints*logoddlik);
		    ptCell->updateOccupancy(numpoints*logoddlik,occupancy_limit);
		    if(ptCell->getOccupancy()<0) ptCell->hasGaussian_ = false;
		}
		else
		{
		    ptCell->updateOccupancy(-0.85*numpoints,occupancy_limit); ///The cell does not have gaussian, so we mark that we saw it empty...
		    //ptCell->updateEmpty(-0.2*numpoints,l2target); ///The cell does not have gaussian, so we mark that we saw it empty...
		}
	    }else{
		ptCell = dynamic_cast<NDTCell*>(index_->addPoint(pt.data())); ///Add fake point to initialize!
	    }
	}
	if(updatePositive){
	    Eigen::Matrix3d ucov = ndts[it]->getCov();
	    float r,g,b;
	    ndts[it]->getRGB(r,g,b);
	    addDistributionToCell(ucov, m, numpoints, r, g,b,maxnumpoints,occupancy_limit); ///<FIXME:: local implementation can be faster?
	}
    }
    //fprintf(stderr,"| Gaussians=%d Invalid=%d |", ndts.size(), num_high);
    for(unsigned int i=0;i<ndts.size();i++) delete ndts[i];
    isFirstLoad_ = false;
}


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/**
 * Adds one measurement to the map using NDT-OM update step
 * @return true if an inconsistency was detected
 */
bool NDTMap::addMeasurement(const Eigen::Vector3d &origin, Eigen::Vector3d endpoint, double classifierTh, double maxz, double sensor_noise)
{

    if(index_ == NULL)
    {
	return false;
    }

    bool retval = false;

    LazyGrid *lz = dynamic_cast<LazyGrid*>(index_);
    double centerX,centerY,centerZ;
    lz->getCenter(centerX, centerY, centerZ);
    double cellSizeX,cellSizeY,cellSizeZ;
    lz->getCellSize(cellSizeX, cellSizeY, cellSizeZ);
    int sizeX,sizeY,sizeZ;
    lz->getGridSize(sizeX, sizeY, sizeZ);

    //Cell<PointT> ****dataArray = lz->getDataArrayPtr();

    double min1 = std::min(cellSizeX,cellSizeY);
    double min2 = std::min(cellSizeZ,cellSizeY);

    double resolution = std::min(min1,min2); ///Select the smallest resolution

    if(lz==NULL)
    {
        fprintf(stderr,"NOT LAZY GRID!!!\n");
        exit(1);
    }
    NDTCell *ptCell=NULL;

    Eigen::Matrix<double,3,1> pt;

    Eigen::Vector3d diff=endpoint-origin;


    double l = diff.norm();
    if(l>200)
    {
        fprintf(stderr,"addMeasurement::Very long distance (%lf) :( \n",l);
        return false;
    }
    if(resolution<0.01)
    {
        fprintf(stderr,"Resolution very very small (%lf) :( \n",resolution);
        return false;
    }

    int NN = l / (resolution);
    if(NN <= 0) return false;
    diff = diff/(float)NN;
    //fprintf(stderr," (N=%d) ",NN);
    bool updatePositive = true;
    if(endpoint(2)>maxz)
    {
        return false;
    }

    int idxo=0, idyo=0,idzo=0;
    for(int i=0; i<NN-2; i++)
    {
        pt = origin + ((float)(i+1)) *diff;
        int idx,idy,idz;

        idx = (int)(((pt(0) - centerX)/cellSizeX+0.5) + sizeX/2);
        idy = (int)(((pt(1) - centerY)/cellSizeY+0.5) + sizeY/2);
        idz = (int)(((pt(2) - centerZ)/cellSizeZ+0.5) + sizeZ/2);


        ///We only want to check every cell once, so
        ///increase the index if we are still in the same cell
        if(idx == idxo && idy==idyo && idz ==idzo)
        {
            continue;
        }
        else
        {
            idxo = idx;
            idyo=idy,idzo=idz;
        }
        /// Check the validity of the index
        //if(idx < sizeX && idy < sizeY && idz < sizeZ && idx >=0 && idy >=0 && idz >=0)
        //{
        //    ptCell = dynamic_cast<NDTCell<PointT> *>  (lz->getCellAt(idx,idy,idz)); //dataArray[idx][idy][idz]);
	lz->getNDTCellAt(idx,idy,idz,ptCell);
        //}
        //else
        //{
        //    continue;
        //}

        if(ptCell != NULL)
        {
            double l2target = 0;

            if(ptCell->hasGaussian_)
            {
                Eigen::Vector3d out,vpt;

                double lik = ptCell->computeMaximumLikelihoodAlongLine(origin.data(), pt.data(), out);
                double dist = (origin-out).norm();
                if(dist > l) continue; ///< don't accept points further than the measurement

                l2target = (out-endpoint).norm(); ///<distance to endpoint
                //double thr =exp(-0.5*(l2target*l2target)/(sensor_noise*sensor_noise)); ///This is the probability of max lik point being endpoint
                //ptCell->updateEmpty(lik*(1-thr),l2target);

                double sigma_dist = 0.5 * (dist/30.0); ///test for distance based sensor noise
                double snoise = sigma_dist + sensor_noise;
                double thr =exp(-0.5*(l2target*l2target)/(snoise*snoise)); ///This is the probability of max lik point being endpoint
                lik *= (1.0-thr);
                lik = 0.2*lik+0.5; ///Evidence value for empty - alpha * p(x);
                double logoddlik = log( (1.0-lik)/(lik) );
                ptCell->updateEmpty(logoddlik,l2target);


                /*
                if(lik>thr){
                	retval = true;
                	ptCell->updateEmpty(lik,l2target);
                }*/
            }
            else
            {
                ptCell->updateEmpty(-0.1,l2target); ///The cell does not have gaussian, so we mark that we saw it empty...
            }
        }
        else
        {
            index_->addPoint(pt.data()); ///Add fake point to initialize!
        }
    }

    if(updatePositive) index_->addPoint(endpoint);

    isFirstLoad_ = false;

    return retval;
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Get estimated depth 
/// This goes through now the cells always to max depth and is not the most efficient solution
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

double NDTMap::getDepth(Eigen::Vector3d origin, Eigen::Vector3d dir, double maxDepth){
	Eigen::Vector3d ray_endpos=origin + dir * maxDepth;
	std::vector< NDTCell*> cells; 
	
	Eigen::Vector3d diff = ray_endpos - origin;
	
	LazyGrid *lz = dynamic_cast<LazyGrid*>(index_);
	if(lz==NULL){
		fprintf(stderr,"NOT LAZY GRID!!!\n");
		exit(1);
	}
	
	if(!lz->traceLine(origin, ray_endpos,diff,1000.0, cells)){
		return maxDepth+1.0;
	}
	//fprintf(stderr,"Got trace with %d Cells (%lf)\n",cells.size(),(ray_endpos-origin).norm());
	
	Eigen::Vector3d out;
	bool hasML = false;
	
	for(unsigned int i=0;i<cells.size();i++){
		 if(cells[i]->hasGaussian_){
			double lik = cells[i]->computeMaximumLikelihoodAlongLine(origin.data(), ray_endpos.data(), out);
			if(lik>0.1){
				//fprintf(stderr,"Got ML %lf (%lf)\n",lik,(out-origin).norm());
				hasML = true;
				break;
			}
		} 
	}
	
	if(hasML) return (out - origin).norm();
	
	return (maxDepth+1.0);
}


double NDTMap::getDepthSmooth(Eigen::Vector3d origin,
                                      Eigen::Vector3d dir,
                                      double maxDepth,
                                      int n_neigh,
                                      double weight,
                                      double threshold,
                                      Eigen::Vector3d *hit)
{
    Eigen::Vector3d ray_endpos=origin + dir * maxDepth;
    std::vector< NDTCell*> cells, surr; 

    Eigen::Vector3d diff = ray_endpos - origin;

    LazyGrid *lz = dynamic_cast<LazyGrid*>(index_);
    if(lz==NULL){
	fprintf(stderr,"NOT LAZY GRID!!!\n");
	exit(1);
    }

    if(!lz->traceLine(origin, ray_endpos,diff,1000.0, cells)){
	return maxDepth+1.0;
    }

    Eigen::Vector3d out;
    bool hasML = false;

    for(unsigned int i = 0; i < cells.size(); i++)
    {
	if(cells[i]->hasGaussian_)
	{
	    surr = lz->getClosestNDTCells(cells[i]->getCenterAddr(), n_neigh, true);
	    double like = cells[i]->computeMaximumLikelihoodAlongLine(origin.data(), ray_endpos.data(), out);
	    for (unsigned int k = 1u; k < surr.size(); ++k)
	    {
		like += weight * surr[k]->getLikelihood(out.data());
	    }			
	    if(like > threshold)
	    {
		hasML = true;
		break;
	    }
	} 
    }

    if (hasML)
    {
	if (hit != NULL)
	{
	    *hit = out;
	}
	return (out - origin).norm();
    }

    return (maxDepth+1.0);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void NDTMap::loadPointCloud(const std::vector<const double* const> &pc, const std::vector<std::vector<size_t> > &indices)
{

  //    loadPointCloud(pc);
    // Specific function related to CellVector
    CellVector *cl = dynamic_cast<CellVector*>(index_);
    if (cl != NULL)
    {
        for (size_t i = 0; i < indices.size(); i++)
        {
            cl->addCellPoints(pc, indices[i]);
        }

    }
    else
    {
        //ERR("loading point clouds using indices are currently supported in CellVector index_.");
    }
}

void NDTMap::computeMaximumLikelihoodPointCloudWithRangePairs(const Eigen::Vector3d &origin, 
                                                              const std::vector<const double* const> &pc,
                                                              const Eigen::Vector3d &virtualOrigin,
                                                              std::vector<Eigen::Vector3d> &pc_out,
                                                              std::vector<std::pair<double,double> > &ranges,
                                                              double max_range = 130.) const
{
    pc_out.clear();
    if(isFirstLoad_ || index_ == NULL)
    {
        return;
    }
    std::vector<const double* const>::const_iterator it = pc.begin();

    LazyGrid *lz = dynamic_cast<LazyGrid*>(index_);
    if(lz==NULL)
    {
        fprintf(stderr,"NOT LAZY GRID!!!\n");
        exit(1);
    }
    NDTCell* ptCell = NULL; 

    // Step through the points, compute the directions and for now do the full line tracing...
    std::vector< NDTCell*> cells; 
    bool updatePositive = true;

    while(it!=pc.end())
    {
        Eigen::Map<const Eigen::Vector3d> pt(*it);
        if(isnan(*it))
        {
            it++;
            continue;
        }
        Eigen::Vector3d diff=pt-origin;
        double l = diff.norm();
        
        cells.clear();

        // Trace line uses the diff to check that we never go beyond the diff value...
        // Hacky for now but force the diff to be a bit larger -> 3 meters...
        diff *= ((l+3.)/l);
        if(!lz->traceLine(virtualOrigin,pt,diff,1000.0,cells)) {
            it++;
            continue;
        }
        
        for(unsigned int i=0; i<cells.size(); i++)
        {
            ptCell = cells[i];
            if(ptCell != NULL)
            {
                double l2target = 0;
                if(ptCell->hasGaussian_)
                {
                    Eigen::Vector3d out;
                    double lik = ptCell->computeMaximumLikelihoodAlongLine(origin.data(), *it, out);
                    //std::cerr << "lik : " << lik << std::endl;
                    if (lik > 0.) {
                        ranges.push_back(std::pair<double,double>(l, (out-origin).norm()));
                        pc_out.push_back(out);
                    }
                }
            }
        }
        it++;
    }
}

/*
void NDTMap::computeConflictingPoints(const Eigen::Vector3d &origin,
                              const pcl::PointCloud<pcl::PointXYZ> &pc,
                              pcl::PointCloud<pcl::PointXYZ> &pc_out,
                              pcl::PointCloud<pcl::PointXYZ> &pc2_out,
                              double likelihoodFactor) const
{
    pc_out.clear(); pc2_out.clear();
    if(isFirstLoad_ || index_ == NULL)
    {
        return;
    }
    pcl::PointCloud<pcl::PointXYZ>::const_iterator it = pc.points.begin();

    LazyGrid *lz = dynamic_cast<LazyGrid*>(index_);
    if(lz==NULL)
    {
        fprintf(stderr,"NOT LAZY GRID!!!\n");
        exit(1);
    }
    pcl::PointXYZ po;
    po.x = origin(0); po.y = origin(1); po.z = origin(2);
    NDTCell* ptCell = NULL; 

    while(it!=pc.points.end())
    {
        if(std::isnan(it->x) ||std::isnan(it->y) ||std::isnan(it->z))
        {
            it++;
            continue;
        }
        
        Eigen::Vector3d diff;
        diff << it->x-origin(0), it->y-origin(1), it->z-origin(2);
        double l = diff.norm();
        
        // Get the cell
        ptCell = lz->getClosestNDTCell(*it);
        if (ptCell == NULL) {
            pc_out.push_back(*it);
            it++;
            continue;
        }
        // if (ptCell->getOccupancyRescaled() < 0.8) {
        //     pc_out.push_back(*it);
        //     it++;
        //     continue;
        // }

        if (!ptCell->hasGaussian_) {
            pc_out.push_back(*it);
            it++;
            continue;
        }

        double lik = ptCell->getLikelihood(*it);
        //        std::cerr << "lik : " << lik << " " << std::flush;
        lik = (lik < 0) ? 0 : lik;
        Eigen::Vector3d out, pt;
        pt << it->x, it->y, it->z;
        double lik_trace = ptCell->computeMaximumLikelihoodAlongLine(po, *it, out);
        //        std::cerr << " lik : " << lik << " lik_trace : " << lik_trace << " lik/lik_trace : " << lik/lik_trace << std::endl;

        // if (lik / lik_trace < likelihoodFactor) {
        //     pc_out.push_back(*it);
        // }
        //        if (lik < likelihoodFactor) {
        //        std::cerr << " " << norm << std::flush;
        double norm = (pt - out).norm();
        //        if (norm > 0.8)
        // if (fabs(pt(2) - out(2)) > 0.2)
        // {
        //     pc_out.push_back(*it);
        //     pc2_out.push_back(pcl::PointXYZ(out(0),out(1),out(2)));
        // }
        pc_out.push_back(pcl::PointXYZ(out(0),out(1),out(2)));
        it++;
    }
}
*/


#if 0
void 
NDTMap::computeMaximumLikelihoodPointRangesForPoseSet(const std::vector<Eigen::Affine3d> &poses, 
                                                       const pcl::PointCloud<pcl::PointXYZ> &pc,
                                                       const Eigen::Vector3d &virtualOrigin,
                                                       Eigen::MatrixXd &predictedRanges, 
                                                      Eigen::VectorXd &rawRanges) const {
    if(isFirstLoad_ || index_ == NULL)
    {
        return;
    }
    LazyGrid *lz = dynamic_cast<LazyGrid*>(index_);
    if(lz==NULL)
    {
        fprintf(stderr,"NOT LAZY GRID!!!\n");
        exit(1);
    }

    pcl::PointCloud<pcl::PointXYZ>::const_iterator it = pc.points.begin();

    pcl::PointXYZ po;
    po.x = origin(0); po.y = origin(1); po.z = origin(2);
    NDTCell* ptCell = NULL; 

    // Step through the points, compute the directions and for now do the full line tracing...
    std::vector< NDTCell*> cells; 
    bool updatePositive = true;

    while(it!=pc.points.end())
    {
        if(std::isnan(it->x) ||std::isnan(it->y) ||std::isnan(it->z))
        {
            it++;
            continue;
        }
        
        Eigen::Vector3d diff;
        diff << it->x-origin(0), it->y-origin(1), it->z-origin(2);
        double l = diff.norm();
        
        cells.clear();

        // Trace line uses the diff to check that we never go beyond the diff value...
        // Hacky for now but force the diff to be a bit larger -> 3 meters...
        diff *= ((l+3.)/l);
        if(!lz->traceLine(virtualOrigin,*it,diff,1000.0,cells)) {
            it++;
            continue;
        }
        
        for(unsigned int i=0; i<cells.size(); i++)
        {
            ptCell = cells[i];
            if(ptCell != NULL)
            {
                double l2target = 0;
                if(ptCell->hasGaussian_)
                {
                    Eigen::Vector3d out;
                    double lik = ptCell->computeMaximumLikelihoodAlongLine(po, *it, out);
                    //std::cerr << "lik : " << lik << std::endl;
                    if (lik > 0.) {
                        ranges.push_back(std::pair<double,double>(l, (out-origin).norm()));
                        pc_out.push_back(pcl::PointXYZ(out[0], out[1], out[2]));
                    }
                }
            }
        }
        it++;
    }
}
#endif

void NDTMap::computeMaximumLikelihoodRanges(const Eigen::Vector3d &origin,
                                            const Eigen::VectorXd &rawRanges,
                                            const std::vector<Eigen::Vector3d> &dirs,
                                            Eigen::VectorXd &ranges) const {
    

    if(isFirstLoad_ || index_ == NULL)
    {
        return;
    }

    LazyGrid *lz = dynamic_cast<LazyGrid*>(index_);
    if(lz==NULL)
    {
        fprintf(stderr,"NOT LAZY GRID!!!\n");
        exit(1);
    }
    NDTCell* ptCell = NULL; 

    // Step through the points, compute the directions and for now do the full line tracing...
    std::vector< NDTCell*> cells; 
    bool updatePositive = true;

    double resolution = this->getSmallestCellSizeInMeters();
    ranges.resize(dirs.size());
    ranges.setZero();
    
    for (unsigned int i = 0; i < dirs.size(); i++)
    {
        cells.clear();

        // TODO update traceLine to use the rawRanges as a bias in selecting the cells.
        // Use the difference to be 2*the resolution.
        // Eigen::Vector3d p1 = origin+(rawRanges[i]-resolution)*dirs[i];
        // Eigen::Vector3d p2 = origin+(rawRanges[i]+resolution)*dirs[i];
        // if(!lz->traceLine(p1, p2, dirs[i]*2*resolution,1000.0,cells)) {
        //     continue;
        // }
        Eigen::Vector3d p1 = origin;
        Eigen::Vector3d p2 = origin+dirs[i]*100.;
        // Eigen::Vector3d p2 = origin+(rawRanges[i]+resolution)*dirs[i];
        if(!lz->traceLine(p1, p2, dirs[i]*100.,1000.0,cells)) {
            continue;
        }
        
        double max_lik = 0.;
        double range = 0.;
        for(unsigned int j=0; j<cells.size(); j++)
        {
            ptCell = cells[j];
            if(ptCell != NULL)
            {
                if(ptCell->hasGaussian_)
                {
                    Eigen::Vector3d out;
                    double lik = ptCell->computeMaximumLikelihoodAlongLine(p1, p2, out);
                    
                    if (lik > max_lik) {
                        max_lik = lik;
                        range = (out-origin).norm();
                        ranges[i] = range;
                    }
                }
            }
        }

        //        std::cout << "ranges[" << i << "] : " << ranges[i] << "   -   " << rawRanges[i] << std::endl;

    }    
}


/** Helper function, computes the  NDTCells
*/
void NDTMap::computeNDTCells(int cellupdatemode, unsigned int maxnumpoints, float occupancy_limit, Eigen::Vector3d origin, double sensor_noise)
{
    CellVector *cv = dynamic_cast<CellVector*>(index_);

    typename std::set<NDTCell*>::iterator it = update_set.begin();
    unsigned int num_updates=update_set.size();
#pragma omp parallel num_threads(N_THREADS)
{
    #pragma omp for
    for(int i=0;i<num_updates;i++)
    {
        auto it = update_set.begin();
        for(int j=0;j<i;j++)
            it++;
        NDTCell *cell = *it;
        if(cell!=NULL)
        {
            cell->computeGaussian(cellupdatemode,maxnumpoints, occupancy_limit, origin,sensor_noise);
            ///Process the conflict points
            if(cell->points_.size()>0)
            {
                cell->points_.clear();
            }
            if (cv!=NULL)
            {
                // Set the mean to the cell's centre.
                Eigen::Vector3d mean = cell->getMean();
                cell->setCenter(mean.data());
            }
        }
    }
}
    update_set.clear();

    CellVector *cl = dynamic_cast<CellVector*>(index_);
    if(cl!=NULL)
    {
        cl->initKDTree();
    }
}



/** 
* Computes the normaldistribution parameters and leaves the points a
*/
void NDTMap::computeNDTCellsSimple(bool keepPoints)
{
    CellVector *cv = dynamic_cast<CellVector*>(index_);
    typename SpatialIndex::CellVectorItr it = index_->begin();

    while (it != index_->end())
    {
        NDTCell *cell = (*it);
        if(cell!=NULL)
        {
            cell->computeGaussianSimple();
          
            if (cv!=NULL)
            {
                // Set the mean to the cell's centre.
                Eigen::Vector3d mean = cell->getMean();
                cell->setCenter(mean.data());
            }
        }
        else
        {

        }
        it++;
    }

    CellVector *cl = dynamic_cast<CellVector*>(index_);
    if(cl!=NULL)
    {
        cl->initKDTree();
    }
}

/** output methods for saving the map in the jff format
 */
int NDTMap::writeToJFF(const char* filename)
{

    if(filename == NULL)
    {
        //ERR("problem outputing to jff\n");
        return -1;
    }

    FILE * jffout = fopen(filename, "w+b");

    fwrite(_JFFVERSION_, sizeof(char), strlen(_JFFVERSION_), jffout);

    switch(this->getMyIndexInt())
    {
    case 1:
        writeCellVectorJFF(jffout);
        break;
    case 2:
        //writeOctTreeJFF(jffout);
        break;
    case 3:
        writeLazyGridJFF(jffout);
        break;
    default:
        //ERR("unknown index type\n");
        return -1;
    }

    fclose(jffout);

    return 0;
}


int NDTMap::writeCellVectorJFF(FILE * jffout)
{
    int indexType[1] = {1};
    fwrite(indexType, sizeof(int), 1, jffout);

    // TODO: add CellVector specific stuff

    typename SpatialIndex::CellVectorItr it = index_->begin();
    while (it != index_->end())
    {
        NDTCell *cell = (*it);
        if(cell!=NULL)
        {
            if(cell->hasGaussian_)
            {
                // TODO: add index specific content smartly
                if(cell->writeToJFF(jffout) < 0)
                    return -1;
            }
        }
        else
        {
            // do nothing
        }
        it++;
    }

    return 0;

}


int NDTMap::writeOctTreeJFF(FILE * jffout)
{
    int indexType[1] = {2};
    fwrite(indexType, sizeof(int), 1, jffout);

    // TODO: add OctTree specific stuff

    typename SpatialIndex::CellVectorItr it = index_->begin();
    while (it != index_->end())
    {
        NDTCell *cell = (*it);
        if(cell!=NULL)
        {
            if(cell->hasGaussian_)
            {
                // TODO: add index specific content smartly
                if(cell->writeToJFF(jffout) < 0)
                    return -1;
            }
        }
        else
        {
            // do nothing
        }
        it++;
    }

    return 0;

}

int NDTMap::writeLazyGridJFF(FILE * jffout)
{
    int indexType[1] = {3};
    fwrite(indexType, sizeof(int), 1, jffout);

    // add LazyGrid specific stuff
    double sizeXmeters, sizeYmeters, sizeZmeters;
    double cellSizeX, cellSizeY, cellSizeZ;
    double centerX, centerY, centerZ;
    LazyGrid *ind = dynamic_cast<LazyGrid*>(index_);

    ind->getGridSizeInMeters(sizeXmeters, sizeYmeters, sizeZmeters);
    ind->getCellSize(cellSizeX, cellSizeY, cellSizeZ);
    ind->getCenter(centerX, centerY, centerZ);

    double lazyGridData[9] = { sizeXmeters, sizeYmeters, sizeZmeters,
                               cellSizeX,   cellSizeY,   cellSizeZ,
                               centerX,     centerY,     centerZ
                             };

    fwrite(lazyGridData, sizeof(double), 9, jffout);

    fwrite(ind->getProtoType(), sizeof(NDTCell), 1, jffout);

    // loop through all active cells
    typename SpatialIndex::CellVectorItr it = index_->begin();
    while (it != index_->end())
    {
	if((*it)->writeToJFF(jffout) < 0)	return -1;
        it++;
    }

    return 0;

}

/** method to load NDT maps from .jff files
USAGE:	create NDTMap with desired index and PointType (index type is
checked, but Point type is NOT checked) via e.g.

lslgeneric::NDTMap<pcl::PointXYZ> nd1(
new lslgeneric::LazyGrid<pcl::PointXYZ>(0.4)); --> (*)

and then call

nd1.loadFromJFF("map0027.jff");

 *) use this constructor so index is not initialized and attributes
 can be set manually
 */
int NDTMap::loadFromJFF(const char* filename)
{

    FILE * jffin;

    if(filename == NULL)
    {
        JFFERR("problem outputing to jff");
    }

    jffin = fopen(filename,"r+b");
    if (jffin == NULL) {
      JFFERR("file not found");
      std::cerr << "Failed to open : " << filename << std::endl;
      return -5;
    }

    char versionBuf[16];
    if(fread(&versionBuf, sizeof(char), strlen(_JFFVERSION_), jffin) <= 0)
    {
        JFFERR("reading version failed");
    }
    versionBuf[strlen(_JFFVERSION_)] = '\0';

    int indexType;
    if(fread(&indexType, sizeof(int), 1, jffin) <= 0)
    {
        JFFERR("reading version failed");
    }

    if(indexType != this->getMyIndexInt())
    {
        switch(indexType)
        {
        case 1:
            std::cerr << "Map uses CellVector\n";
            return -1;
            break;
        case 2:
            std::cerr << "Map uses OctTree\n";
            return -2;
            break;
        case 3:
            std::cerr << "Map uses LazyGrid\n";
            return -3;
            break;
        }
    }

    switch(indexType)
    {
    case 1:
    {
        CellVector* cv = dynamic_cast<CellVector * >(index_);
        if(cv->loadFromJFF(jffin) < 0)
        {
            JFFERR("Error loading CellVector");
        }
        break;
    }
#if 0
    case 2:
    {
        OctTree* tr = dynamic_cast<OctTree*>(index_);
        if(tr->loadFromJFF(jffin) < 0)
        {
            JFFERR("Error loading OctTree");
        }
        break;
    }
#endif
    case 3:
    {
	std::cerr << "Map uses LazyGrid\n";
        LazyGrid* gr = dynamic_cast<LazyGrid*>(index_);
        if(gr->loadFromJFF(jffin) < 0)
        {
            JFFERR("Error loading LazyGrid");
        }
        break;
    }
    default:
        JFFERR("error casting index");
    }

    NDTCell *ptCell = new NDTCell();
    index_->setCellType(ptCell);
    delete ptCell;

    fclose(jffin);

   // std::cout << "map loaded successfully " << versionBuf << std::endl;

    isFirstLoad_ = false;

    return 0;

}

/** method to load NDT maps from .jff files
USAGE:	create NDTMap with desired index and PointType (index type is
checked, but Point type is NOT checked) via e.g.

lslgeneric::NDTMap<pcl::PointXYZ> nd1(
new lslgeneric::LazyGrid<pcl::PointXYZ>(0.4)); --> (*)

and then call

nd1.loadFromJFF("map0027.jff");

 *) use this constructor so index is not initialized and attributes
 can be set manually
 */
int NDTMap::loadFromJFF(FILE * jffin)
{


    char versionBuf[16];
    if(fread(&versionBuf, sizeof(char), strlen(_JFFVERSION_), jffin) <= 0)
    {
        JFFERR("reading version failed");
    }
    versionBuf[strlen(_JFFVERSION_)] = '\0';

    int indexType;
    if(fread(&indexType, sizeof(int), 1, jffin) <= 0)
    {
        JFFERR("reading version failed");
    }

    if(indexType != this->getMyIndexInt())
    {
        switch(indexType)
        {
        case 1:
            std::cerr << "Map uses CellVector\n";
            return -1;
            break;
        case 2:
            std::cerr << "Map uses OctTree\n";
            return -2;
            break;
        case 3:
            std::cerr << "Map uses LazyGrid\n";
            return -3;
            break;
        }
    }

    switch(indexType)
    {
    case 1:
    {
        CellVector* cv = dynamic_cast<CellVector * >(index_);
        if(cv->loadFromJFF(jffin) < 0)
        {
            JFFERR("Error loading CellVector");
        }
        break;
    }
#if 0
    case 2:
    {
        OctTree* tr = dynamic_cast<OctTree*>(index_);
        if(tr->loadFromJFF(jffin) < 0)
        {
            JFFERR("Error loading OctTree");
        }
        break;
    }
#endif
    case 3:
    {
	std::cerr << "Map uses LazyGrid\n";
        LazyGrid* gr = dynamic_cast<LazyGrid*>(index_);
        if(gr->loadFromJFF(jffin) < 0)
        {
            JFFERR("Error loading LazyGrid");
        }
        break;
    }
    default:
        JFFERR("error casting index");
    }

    NDTCell *ptCell = new NDTCell();
    index_->setCellType(ptCell);
    delete ptCell;


   // std::cout << "map loaded successfully " << versionBuf << std::endl;

    isFirstLoad_ = false;

    return 0;

}


/// returns the current spatial index as a string (debugging function)
std::string NDTMap::getMyIndexStr() const
{
    CellVector* cl = dynamic_cast<CellVector * >(index_);
    if(cl!=NULL)
    {
        return std::string("CellVector");
    }
#if 0
    OctTree<PointT>* tr = dynamic_cast<OctTree<PointT>*>(index_);
    if(tr!=NULL)
    {
        return std::string("OctTree");
    }
#endif
    LazyGrid *gr = dynamic_cast<LazyGrid*>(index_);
    if(gr!=NULL)
    {
        return std::string("LazyGrid<PointT>");
    }

    return std::string("Unknown index type");
}

/// returns the current spatial index as an integer (debugging function)
int NDTMap::getMyIndexInt() const
{
    CellVector* cl = dynamic_cast<CellVector * >(index_);
    if(cl!=NULL)
    {
        return 1;
    }
#if 0
    OctTree<PointT>* tr = dynamic_cast<OctTree<PointT>*>(index_);
    if(tr!=NULL)
    {
        return 2;
    }
#endif
    LazyGrid *gr = dynamic_cast<LazyGrid*>(index_);
    if(gr!=NULL)
    {
        return 3;
    }

    return -1;
}

//computes the *negative log likelihood* of a single observation
double NDTMap::getLikelihoodForPoint(const double * const pt)
{
    //assert(false);
    double uniform=0.00100;
    NDTCell* ndCell = NULL;
    
#if 0
    OctTree<PointT>* tr = dynamic_cast<OctTree<PointT>*>(index_);

    if(tr==NULL)
    {
#endif
	LazyGrid *gr = dynamic_cast<LazyGrid*>(index_);
	if(gr==NULL)
	{
	    //cout<<"bad index - getLikelihoodForPoint\n";
	    return uniform;
	}
	ndCell = gr->getClosestNDTCell(pt);
#if 0
    }
    else
    {
        ndCell = tr->getClosestNDTCell(pt);
    }
#endif
    if(ndCell == NULL) return uniform;

    double prob = ndCell->getLikelihood(pt);
    prob = (prob<0) ? 0 : prob; //uniform!! TSV
    return prob;
}

/*
//use trilinear interpolation from available immediate neighbors
template<typename PointT>
double NDTMap<PointT>::getLikelihoodForPointWithInterpolation(PointT pt) {

    //ATM only for grid map
    //     tll------tlr
    //     /|       /|
    //    / |      / |
    //  tul------tur |    z
    //   | bll----|-blr   ^  y
    //   | /      | /     | /
    //   |/       |/      |/
    //  bul------bur      ---> x
    double uniform=0;//0.00100;
    Cell* cell = NULL;
    NDTCell<PointT>* ndCell = NULL;
    double cumProb = 0;
    double weight = 0;
    int evals = 1;

    LazyGrid<PointT> *gr = dynamic_cast<LazyGrid<PointT>*>(index_);
    if(gr==NULL) {
	//cout<<"bad index - getLikelihoodForPointWithInterpolation\n";
	return uniform;
    }
    cell = gr->getCellForPoint(pt);
    if(cell == NULL) return uniform;


    //get coordinates of cell
    int indXn, indYn, indZn;
    PointT centerGrid, sizeCell, centerCell;
    int sizeGridX, sizeGridY,sizeGridZ;
    centerCell = cell->getCenter();
    gr->getCenter(centerGrid.x,centerGrid.y,centerGrid.z);
    gr->getGridSize(sizeGridX,sizeGridY,sizeGridZ);
    gr->getCellSize(sizeCell.x,sizeCell.y,sizeCell.z);
    gr->getIndexForPoint(pt,indXn,indYn,indZn);

    double x,y,z;
    x = (pt.x - centerCell.x)/sizeCell.x;
    y = (pt.y - centerCell.y)/sizeCell.y;
    z = (pt.z - centerCell.z)/sizeCell.z;
    if(x <0 ) x = 0;
    if(y <0 ) y = 0;
    if(z <0 ) z = 0;
    if(x >1 ) x = 1;
    if(y >1 ) y = 1;
    if(z >1 ) z = 1;

    //bul
    double prob = 0;
    ndCell = dynamic_cast<NDTCell<PointT>*> (cell);
    if(ndCell != NULL) {
	prob = ndCell->getLikelihood(pt);
	prob = (prob<0) ? uniform : prob;
	weight = (1 - x + 1 - y + 1 - z)/(3.0);
	if(weight < 0) cerr<<weight<<endl;
	cumProb += prob*weight;
	//cout<<"\t"<<weight<<" "<<prob<<" --> "<<cumProb<<endl;
	evals++;
    }

    //tul
    Cell* c = gr->getCellAt(indXn,indYn,indZn+1);
    if(c != NULL) {
	ndCell = dynamic_cast<NDTCell<PointT>*> (c);
	if(ndCell != NULL) {
	    prob = ndCell->getLikelihood(pt);
	    prob = (prob<0) ? uniform : prob;
	    weight = (1 - x + 1 - y + z)/(3.0);
	    if(weight < 0) cerr<<weight<<endl;
	    cumProb += prob*weight;
	    //cout<<"\t"<<weight<<" "<<prob<<" --> "<<cumProb<<endl;
	    evals++;
	}
    }
    //tur
    c = gr->getCellAt(indXn+1,indYn,indZn+1);
    if(c != NULL) {
	ndCell = dynamic_cast<NDTCell<PointT>*> (c);
	if(ndCell != NULL) {
	    prob = ndCell->getLikelihood(pt);
	    prob = (prob<0) ? uniform : prob;
	    weight = (x + 1-y + z )/(3.0);
	    if(weight < 0) cerr<<weight<<endl;
	    cumProb += prob*weight;
	    //cout<<"\t"<<weight<<" "<<prob<<" --> "<<cumProb<<endl;
	    evals++;
	}
    }
    //tll
    c = gr->getCellAt(indXn,indYn+1,indZn+1);
    if(c != NULL) {
	ndCell = dynamic_cast<NDTCell<PointT>*> (c);
	if(ndCell != NULL) {
	    prob = ndCell->getLikelihood(pt);
	    prob = (prob<0) ? uniform : prob;
	    weight = (1-x + y + z )/(3.0);
	    if(weight < 0) cerr<<weight<<endl;
	    cumProb += prob*weight;
	    //cout<<"\t"<<weight<<" "<<prob<<" --> "<<cumProb<<endl;
	    evals++;
	}
    }
    //tlr
    c = gr->getCellAt(indXn+1,indYn+1,indZn+1);
    if(c != NULL) {
	ndCell = dynamic_cast<NDTCell<PointT>*> (c);
	if(ndCell != NULL) {
	    prob = ndCell->getLikelihood(pt);
	    prob = (prob<0) ? uniform : prob;
	    weight = (x + y + z )/(3.0);
	    if(weight < 0) cerr<<weight<<endl;
	    cumProb += prob*weight;
	    //cout<<"\t"<<weight<<" "<<prob<<" --> "<<cumProb<<endl;
	    evals++;
	}
    }
    //bur
    c = gr->getCellAt(indXn+1,indYn,indZn);
    if(c != NULL) {
	ndCell = dynamic_cast<NDTCell<PointT>*> (c);
	if(ndCell != NULL) {
	    prob = ndCell->getLikelihood(pt);
	    prob = (prob<0) ? uniform : prob;
	    weight = (x + 1-y + 1-z )/(3.0);
	    if(weight < 0) cerr<<weight<<endl;
	    cumProb += prob*weight;
	    //cout<<"\t"<<weight<<" "<<prob<<" --> "<<cumProb<<endl;
	    evals++;
	}
    }
    //bll
    c = gr->getCellAt(indXn,indYn+1,indZn);
    if(c != NULL) {
	ndCell = dynamic_cast<NDTCell<PointT>*> (c);
	if(ndCell != NULL) {
	    prob = ndCell->getLikelihood(pt);
	    prob = (prob<0) ? uniform : prob;
	    weight = (1-x + y + 1-z )/(3.0);
	    if(weight < 0) cerr<<weight<<endl;
	    cumProb += prob*weight;
	    //cout<<"\t"<<weight<<" "<<prob<<" --> "<<cumProb<<endl;
	    evals++;
	}
    }
    //blr
    c = gr->getCellAt(indXn+1,indYn+1,indZn);
    if(c != NULL) {
	ndCell = dynamic_cast<NDTCell<PointT>*> (c);
	if(ndCell != NULL) {
	    prob = ndCell->getLikelihood(pt);
	    prob = (prob<0) ? uniform : prob;
	    weight = (x + y + 1-z )/(3.0);
	    if(weight < 0) cerr<<weight<<endl;
	    cumProb += prob*weight;
	    //cout<<"\t"<<weight<<" "<<prob<<" --> "<<cumProb<<endl;
	    evals++;
	}
    }

    //cout<<"== "<<cumProb<<endl;
    return cumProb;
}
*/

std::vector<NDTCell*> NDTMap::getInitializedCellsForPoint(const double* const pt) const
{
    std::vector<NDTCell*> cells;
    LazyGrid *gr = dynamic_cast<LazyGrid*>(index_);
    if(gr==NULL)
    {
        //cout<<"bad index - getCellsForPoint\n";
        return cells;
    }
    cells = gr->getClosestCells(pt);
    return cells;

}

std::vector<NDTCell*> NDTMap::getCellsForPoint(const double* const pt, int n_neigh, bool checkForGaussian) const
{
    //assert(false);
    std::vector<NDTCell*> cells;
    LazyGrid *gr = dynamic_cast<LazyGrid*>(index_);
    if(gr==NULL)
    {
        //cout<<"bad index - getCellsForPoint\n";
        return cells;
    }
    cells = gr->getClosestNDTCells(pt,n_neigh,checkForGaussian);
    return cells;

    //OctTree<PointT>* tr = dynamic_cast<OctTree<PointT>*>(index_);
    //if(tr==NULL) {
    //}
}

bool NDTMap::getCellForPoint(const double* const pt, NDTCell* &out_cell, bool checkForGaussian) const
{

    out_cell = NULL;
    CellVector *cl = dynamic_cast<CellVector*>(index_);
    if(cl!=NULL)
    {
        out_cell = cl->getClosestNDTCell(pt);
        return true;
    }
#if 0
    OctTree<PointT>* tr = dynamic_cast<OctTree<PointT>*>(index_);
    if(tr!=NULL)
    {
        out_cell = tr->getClosestNDTCell(pt);
        return true;
    }
#endif
    LazyGrid *gr = dynamic_cast<LazyGrid*>(index_);
    if(gr!=NULL)
    {
        out_cell = gr->getClosestNDTCell(pt,checkForGaussian);
        return true;
    }
    //cout<<"bad index - getCellForPoint\n";
    return false;
}
void NDTMap::transformNDTMap(Eigen::Transform<double, 3, Eigen::Affine, Eigen::ColMajor> T)
{
    std::vector<NDTCell*> previous_cells = pseudoTransformNDT(T);
    LazyGrid* lz = dynamic_cast<LazyGrid *>(index_);
    Eigen::Vector3d tr = T * lz->translation;
    lz->initialize();
    lz->setSensorPose(tr.data());
    for(std::vector<NDTCell*>::iterator it=previous_cells.begin();it!=previous_cells.end();++it)
    {
        NDTCell* cell = (*it);
        Eigen::Vector3d center=cell->getCenter();
        center=T*center;
        cell->setCenter(center.data());
        lz->addNDTCell(cell);
    }
}
/*
void NDTMap::transformNDTMap(Eigen::Transform<double, 3, Eigen::Affine, Eigen::ColMajor> T)
{
    assert(index_!=NULL);
    std::cerr<<"START_TRANSFORM1"<<std::endl;
    SpatialIndex *si = index_->clone();
    LazyGrid* lzN = dynamic_cast<LazyGrid *>(si);
    LazyGrid* lz = dynamic_cast<LazyGrid *>(index_);
    Eigen::Vector3d tr = T * lz->translation;
	NDTCell *ptCell = new NDTCell();
	si->setCellType(ptCell);
	delete ptCell;
    lzN->setCenter(0,0,0);
    lzN->setSize(map_sizex, map_sizey, map_sizez);
    lzN->setSensorPose(tr.data());
    lzN->cloud_index=lz->cloud_index;
    typename SpatialIndex::CellVectorItr it = index_->begin();
    std::cerr<<"START_TRANSFORM2"<<std::endl;
    while (it != index_->end())
    {
        NDTCell *cell = (*it);
        if(cell!=NULL)
        {
            if(cell->hasGaussian_)
            {
                Eigen::Vector3d mean = cell->getMean();
                Eigen::Matrix3d cov = cell->getCov();
                mean = T*mean;
                ///NOTE: The rotation of the covariance fixed by Jari 6.11.2012
                cov = T.rotation()*cov*T.rotation().transpose();
                NDTCell* nd = (NDTCell*)cell->clone();
                nd->setMean(mean);
                nd->setCov(cov);
                lzN->addNDTCell(cell);
            }
        }
        it++;
    }
    std::cerr<<"END_TRANSFORM2"<<std::endl;
    assert(index_!=NULL);
    delete index_;
    index_ = si;
    std::cerr<<"END_TRANSFORM1"<<std::endl;
}
*/
NDTMap* NDTMap::pseudoTransformNDTMap(Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> T)
{
    NDTMap* map = new NDTMap(new CellVector(), true);
    CellVector* idx = dynamic_cast<CellVector*> (map->getMyIndex());
    typename SpatialIndex::CellVectorItr it = index_->begin();

    while (it != index_->end())
    {
        NDTCell *cell = (*it);
	if(cell->hasGaussian_)
	{

	    Eigen::Vector3d mean = cell->getMean();
	    Eigen::Matrix3d cov = cell->getCov();
	    mean = T*mean;
	    ///NOTE: The rotation of the covariance fixed by Jari 6.11.2012
	    cov = T.rotation()*cov*T.rotation().transpose();
	    NDTCell* nd = (NDTCell*)cell->clone();
	    nd->setMean(mean);
	    nd->setCov(cov);
	    idx->addNDTCell(nd);
	}
        it++;
    }
    return map;
}

std::vector<NDTCell*> NDTMap::pseudoTransformNDT(Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> T) const
{

    std::vector<NDTCell*> ret;
    typename SpatialIndex::CellVectorConstItr it = index_->begin();
    while (it != index_->end())
    {
        NDTCell *cell = (*it);
        if(cell!=NULL)
        {
            if(cell->hasGaussian_)
            {
                Eigen::Vector3d mean = cell->getMean();
                Eigen::Matrix3d cov = cell->getCov();
                mean = T*mean;
                ///NOTE: The rotation of the covariance fixed by Jari 6.11.2012
                cov = T.rotation()*cov*T.rotation().transpose();
                NDTCell* nd = (NDTCell*)cell->clone();
                nd->cloud_index=cell->cloud_index;
                nd->setMean(mean);
                nd->setCov(cov);
                ret.push_back(nd);
            }
        }
        else
        {
            //ERR("problem casting cell to NDT!\n");
        }
        it++;
    }
    return ret;
}

NDTCell* NDTMap::getCellIdx(unsigned int idx) const
{
    CellVector *cl = dynamic_cast<CellVector*>(index_);
    if (cl != NULL)
    {
        return cl->getCellIdx(idx);
    }
    return NULL;
}

std::vector<perception_oru::NDTCell*> NDTMap::getAllCells() const
{

    std::vector<NDTCell*> ret;
    typename SpatialIndex::CellVectorItr it = index_->begin();
    while (it != index_->end())
    {
        NDTCell *cell = (*it);
	if(cell->hasGaussian_)
	{
	    NDTCell* nd = cell->copy();
	    ret.push_back(nd);
	}
        it++;
    }
    return ret;
}

std::vector< boost::shared_ptr<perception_oru::NDTCell> > NDTMap::getAllCellsShared() const
{

    std::vector< boost::shared_ptr< NDTCell > > ret;
    typename SpatialIndex::CellVectorItr it = index_->begin();
    while (it != index_->end())
    {
		
        NDTCell *cell = (*it);
		if(cell->hasGaussian_)
		{
			
			NDTCell* nd = cell->copy();
			boost::shared_ptr< NDTCell > smart_pointer(nd);
// 			NDTCell** ndd = &nd;
			ret.push_back(smart_pointer);
		}
        it++;
    }
    std::cout << "Return " << ret.size() << std::endl;
    return ret;
}

std::vector<perception_oru::NDTCell*> NDTMap::getAllInitializedCells() const
{
    std::vector<NDTCell*> ret;
    typename SpatialIndex::CellVectorItr it = index_->begin();
    while (it != index_->end())
    {
	NDTCell* nd = (*it)->copy();
	ret.push_back(nd);
        it++;
    }
    return ret;
}

std::vector< boost::shared_ptr<perception_oru::NDTCell> > NDTMap::getAllInitializedCellsShared() const
{
    std::vector<boost::shared_ptr<NDTCell> > ret;
    typename SpatialIndex::CellVectorItr it = index_->begin();
    while (it != index_->end())
    {
		NDTCell* nd = (*it)->copy();
		boost::shared_ptr< NDTCell > smart_pointer(nd);
		ret.push_back(smart_pointer);
        it++;
    }
    return ret;
}

int NDTMap::numberOfActiveCells()
{
    int ret = 0;
    if(index_ == NULL) return ret;
    typename SpatialIndex::CellVectorItr it = index_->begin();
    while (it != index_->end())
    {
	if((*it)->hasGaussian_)
	{
	    ret++;
	}
	it++;
    }
    return ret;
}

int NDTMap::numberOfActiveCells() const 
{
    int ret = 0;
    if(index_ == NULL) return ret;
    typename SpatialIndex::CellVectorItr it = index_->begin();
    while (it != index_->end())
    {
	if((*it)->hasGaussian_)
	{
	    ret++;
	}
	it++;
    }
    return ret;
}
NDTCell* NDTMap::getCellAtID(int x,int y,int z) const {
    NDTCell* cell;
    LazyGrid *lz = dynamic_cast<LazyGrid*>(index_);
    lz->getCellAt(x,y,z,cell);
    return cell;
}

bool NDTMap::insertCell(NDTCell cell){
    LazyGrid* gr = dynamic_cast<LazyGrid*>(index_);
    gr->insertCell(*cell.copy());
    return true;
}
std::string NDTMap::ToString(){
   std::stringstream ss;
   if(index_!=NULL)
    ss<<"NDTMap: index= \n"<<index_->ToString()<<std::endl;
return ss.str();


}
}
