#include <Eigen/Eigen>
#include <pcl/point_cloud.h>
#include <pcl/features/feature.h>

#include <string>
#include <climits>

#include <ndt_map/ndt_map_hmt.h>
#include <ndt_map/lazy_grid.h>
#include <ndt_map/cell_vector.h>

#include <cstring>
#include <cstdio>

namespace lslgeneric
{

void NDTMapHMT::initializeGrids() {
   
    if(grids_init) return;

    LazyGrid *proto = dynamic_cast<LazyGrid*>(index_);
    if(proto == NULL) return;
    
    double centerX,centerY,centerZ;
    proto->getCenter(centerX, centerY, centerZ);
    double sizeX,sizeY,sizeZ;
    proto->getGridSizeInMeters(sizeX, sizeY, sizeZ);
    std::cout<<"inti grids: res="<<resolution<<" cen "<<centerX<<" "<<centerY<<" "<<centerZ<<" size "<<sizeX<<" "<<sizeY<<" "<<sizeZ<<std::endl;

    for(int i=-1; i<2; i++) {
	for(int j=-1; j<2; j++) {
	    if(i==0 && j==0) {
		//no need to allocate, this is us
		LazyGrid *lz = dynamic_cast<LazyGrid*>(index_);
		grid_[i+1][j+1] = lz; 
	    } else {
		double cenX,cenY;
		cenX = centerX+(double)i*sizeX;
		cenY = centerY+(double)j*sizeY;
		std::cout<<i<<":"<<j<<" center "<<cenX<<" "<<cenY<<std::endl;
		NDTCell *ptCell = new NDTCell();
		LazyGrid *lz = new LazyGrid(resolution);
		lz->setCellType(ptCell);
		lz->setCenter(cenX,cenY,centerZ);
		lz->setSize(sizeX,sizeY,sizeZ);
		lz->initializeAll();
		grid_[i+1][j+1] = lz; 
		delete ptCell;
	    }
	}
    }
    grids_init = true;

}

bool NDTMapHMT::tryLoadPosition(const Eigen::Vector3d &newPos) {
    //open meta file and read grid centers
    if(my_directory == "" || !grids_init)
    {
	std::cout<<"cannot load from directory!\n";
	return false;
    }
    
    LazyGrid *proto = dynamic_cast<LazyGrid*>(index_);
    if(proto == NULL) return false;
    double sizeX,sizeY,sizeZ;
    double centerX,centerY,centerZ;
    proto->getGridSizeInMeters(sizeX, sizeY, sizeZ);

    std::string meta = my_directory;
    meta += "/metadata.txt";
    //std::cerr<<"metadata file at "<<meta<<std::endl;
    FILE* meta_f = fopen(meta.c_str(),"a+");
    if(meta_f==0) return false;
    char *line = NULL;
    size_t len;
    bool found=false;
    //std::cerr<<"reading metadata file at "<<meta<<std::endl;
    //read in all metadata
    //NOTE: v 2.0 -> Added header:
    //VERSION X.x\n
    //SIZE X.x\n
    //each next line: centerx centery centerz filename\n
    if(getline(&line,&len,meta_f) >0 ) {
	char *tk = strtok(line," ");
	if(tk==NULL) return false;
	if (strncmp(tk,"VERSION",7) == 0) {
	    tk = strtok(NULL," ");
	    if(tk==NULL) return false;
	    if(strncmp(tk,"2.0",3) == 0) {
		if(!getline(&line,&len,meta_f) >0 ) {
		    return false;
		}
		tk = strtok(line," ");
		if(tk==NULL) return false;
		if(strncmp(tk,"SIZE",4) != 0) return false;
		tk = strtok(NULL," ");
		double sizeMeta = atof(tk);
		if(fabsf(sizeMeta - sizeX) > 0.01) {
		    std::cerr<<"cannot load map, different grid size used... reverting to empty map\n";
		    return false;
		}
	    }
	    //add cases ...
	    
	} else {
	    //this is a version 1.0 file, close and re-open and go on
	    std::cerr<<"metafile version 1.0, no protection against different grid size\n";
	    fclose(meta_f);
	    meta_f = fopen(meta.c_str(),"a+");
	}
    }

    while(getline(&line,&len,meta_f) >0 )
    {
	pcl::PointXYZ cen;
	char *token = strtok(line," ");
	if(token==NULL) return -1;
	cen.x = atof(token);
	token = strtok(NULL," ");
	if(token==NULL) return -1;
	cen.y = atof(token);
	token = strtok(NULL," ");
	if(token==NULL) return -1;
	cen.z = atof(token);
	
	token = strtok(NULL," ");
	if(token==NULL) return -1;
	
	if(fabsf(newPos(0)-cen.x) < sizeX/2 && 
		fabsf(newPos(1)-cen.y) < sizeY/2 && fabsf(newPos(2)-cen.z) < sizeZ/2) {
	    found = true;
	    centerX = cen.x; centerY = cen.y; centerZ = cen.z;
	    //std::cerr<<"file already for "<<jffname<<std::endl;
	    break;
	}   
		
    }
    fclose(meta_f);
    if(!found) {
	std::cerr<<"Map file not found!\n";
	return false;
    }
    
    //newPos is inside one of the maps, load that for the center
    LazyGrid *tmp_grid[3][3];
    for(int i=-1; i<2; i++) {
	for(int j=-1; j<2; j++) {
	    double cenX,cenY;
	    cenX = centerX+(double)i*sizeX;
	    cenY = centerY+(double)j*sizeY;
	    std::cout<<i<<":"<<j<<" NEW center "<<cenX<<" "<<cenY<<std::endl;
	    if(this->tryLoad(cenX,cenY,centerZ,tmp_grid[i+1][j+1])) {
		delete grid_[i+1][j+1];
		grid_[i+1][j+1] = tmp_grid[i+1][j+1];
	    } else {	
		grid_[i+1][j+1]->setCenter(cenX,cenY,centerZ);
	    }
	}
    }
    return true;
}

void NDTMapHMT::setInsertPosition(const Eigen::Vector3d &newPos) {

    last_insert = newPos;
    pcl::PointXYZ newPosP;
    newPosP.x = newPos(0);
    newPosP.y = newPos(1);
    newPosP.z = newPos(2);
    //check if newPos is outside current grid
    if(grid_[1][1]->isInside(newPosP)) return;

    std::cout<<"We are outside the central grid, time to switch pointers\n";
    //if yes, save maps
    this->writeTo();

    //find new center grid
    int qi=0,qj=0;
    for(int i=-1; i<2; i++) {
	for(int j=-1; j<2; j++) {
	    if(grid_[i+1][j+1]->isInside(newPosP)) {
		//std::cout<<"switching to grid "<<i+1<<" "<<j+1<<std::endl;
		qi=i; qj=j;
	    }
	}
    }
    LazyGrid *tmp_grid[3][3];
    bool copy[3][3];
    
    double centerX,centerY,centerZ;
    grid_[qi+1][qj+1]->getCenter(centerX, centerY, centerZ);
    double sizeX,sizeY,sizeZ;
    grid_[qi+1][qj+1]->getGridSizeInMeters(sizeX, sizeY, sizeZ);
    for(int i=0; i<3; i++) {
	for(int j=0; j<3; j++) {
	    copy[i][j]=false;
	}
    }
    int oi,oj;
    for(int i=0; i<3; i++) {
	for(int j=0; j<3; j++) {
	    oi = i+qi; 
	    oj = j+qj;
	    if(oi>=0 && oi<3 && oj>=0 &&oj<3) {
		//copy old pointer
		std::cout<<"will copy: "<<oi<<" "<<oj<<" to "<<i<<" "<<j<<std::endl;
		copy[oi][oj] = true;
		tmp_grid[i][j] = grid_[oi][oj];
	    } else {
		//try load of previous map 
		double cenX,cenY;
		cenX = centerX+(double)(i-1)*sizeX;
		cenY = centerY+(double)(j-1)*sizeY;
		if(!this->tryLoad(cenX,cenY,centerZ,tmp_grid[i][j])) {	
		    //de-allocate old pointer, create new one
		    std::cout<<"will allocate new at "<<i<<" "<<j<<std::endl;

		    std::cout<<i<<":"<<j<<" center "<<cenX<<" "<<cenY<<std::endl;
		    NDTCell *ptCell = new NDTCell();
		    LazyGrid *lz = new LazyGrid(resolution);
		    lz->setCellType(ptCell);
		    lz->setCenter(cenX,cenY,centerZ);
		    lz->setSize(sizeX,sizeY,sizeZ);
		    lz->initializeAll();
		    tmp_grid[i][j] = lz; 
		    delete ptCell;
		} else {
		    //std::cout<<"loading previous map\n";
		}
	    }
	}
    }

    //deallocate old border
    for(int i=0; i<3; i++) {
	for(int j=0; j<3; j++) {
	    if(!copy[i][j]) {
		//std::cout<<"DELETE old grid at "<<i<<" "<<j<<std::endl;
		delete grid_[i][j];
	    }
	    grid_[i][j] = tmp_grid[i][j];
	}
    }

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
void NDTMapHMT::loadPointCloud(const pcl::PointCloud<pcl::PointXYZ> &pc, double range_limit)
{

    //std::cout<<"LOAD pc\n";
    typename pcl::PointCloud<pcl::PointXYZ>::const_iterator it = pc.points.begin();
    NDTCell *ptCell;
    
    while(it!=pc.points.end())
    {
        Eigen::Vector3d d;
        if(std::isnan(it->x) ||std::isnan(it->y) ||std::isnan(it->z))
        {
            it++;
            continue;
        }
        /*if(range_limit>0)
        {
            d << it->x, it->y, it->z;
            if(d.norm()>range_limit)
            {
                it++;
                continue;
            }
        }*/
	
	for(int i=0; i<3; ++i) {
	    for(int j=0; j<3; ++j) {
		if(grid_[i][j]->isInside(*it)) {
		    ptCell = (grid_[i][j]->addPoint(*it));
		    if(ptCell!=NULL) {
			update_set.insert(ptCell);
		    }
		    break;
		}
	    }
	}
	it++;
    }
    isFirstLoad_ = false;
}

/**
 * Adds a new cloud: NDT-OM update step
 */
void NDTMapHMT::addPointCloud(const Eigen::Vector3d &origin, const pcl::PointCloud<pcl::PointXYZ> &pc, double classifierTh, double maxz, double sensor_noise, double occupancy_limit)
{
    //std::cerr<<"ADDPointCloud\n";
    if(isFirstLoad_)
    {
	//std::cout<<"First add, loading pcs\n";
	this->loadPointCloud(pc);
	return;
    }
    
    this->setInsertPosition(origin);
    pcl::PointCloud<pcl::PointXYZ>::const_iterator it = pc.points.begin();
    
    LazyGrid *lz = grid_[1][1];
    
    double centerX, centerY, centerZ;
    double sizeXmeters, sizeYmeters, sizeZmeters;
    double eps = 0.001;
    Eigen::Vector3d diff2;
    lz->getGridSizeInMeters(sizeXmeters, sizeYmeters, sizeZmeters);
    lz->getCenter(centerX, centerY, centerZ);
    pcl::PointXYZ po, pt;
    Eigen::Vector3d pf;
    po.x = origin(0); po.y = origin(1); po.z = origin(2);
    NDTCell* ptCell = NULL; 

    std::vector< NDTCell*> cells; 
    bool updatePositive = true;
    //double max_range = 200.;

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

	if(l>max_range_)
	{
	    //fprintf(stderr,"Very long distance (%lf) :( \n",l);
	    it++;
	    continue;
	}

	cells.clear();
	if(!lz->traceLineWithEndpoint(origin,*it,diff,maxz,cells,pf)) {
	    it++;
	    continue;
	}
	//check if pf is within one resolution of diff
	diff2 << it->x-pf(0), it->y-pf(1), it->z-pf(2);
	if((diff2).norm() > eps) {
	    //line continues outside center map, need to trace more!
	    int i=0,j=0;
	    if(it->x > centerX+sizeXmeters/2-resolution/2) {
		i = 1;
	    }
	    if(it->x < centerX-sizeXmeters/2-resolution/2) {
		i = -1;
	    }
	    if(it->y > centerY+sizeYmeters/2-resolution/2) {
		j = 1;
	    }
	    if(it->y < centerY-sizeYmeters/2-resolution/2) {
		j = -1;
	    }
	    std::vector< NDTCell*> cells2; 
	    if(!grid_[i+1][j+1]->traceLineWithEndpoint(pf,*it,diff2,maxz,cells,pf)) {
		std::cerr<<"ray outside central map: at "<<it->x<<" "<<it->y<<" "<<it->z<<" with diff "<<diff.transpose()<<" starts at " <<origin.transpose()
		    <<" stops at "<<pf.transpose()<<" diff2 is "<<diff2.transpose()<<" check grid "<<i<<" "<<j<<std::endl;
		it++;
		continue;
	    }
	    cells.insert(cells.begin(),cells2.begin(),cells2.end());
	    //TODO: we might need to make one more raytrace here
	}

	for(unsigned int i=0; i<cells.size(); i++)
	{
	    ptCell = cells[i];
	    if(ptCell != NULL)
	    {
		double l2target = 0;
		if(ptCell->hasGaussian_)
		{
		    Eigen::Vector3d out, pend,vpt;
		    pend << it->x,it->y,it->z;
		    double lik = ptCell->computeMaximumLikelihoodAlongLine(po, *it, out);
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
		    ptCell->updateOccupancy(logoddlik, occupancy_limit);
		    if(ptCell->getOccupancy()<=0) ptCell->hasGaussian_ = false; 
		}
		else
		{
		    ptCell->updateOccupancy(-0.2, occupancy_limit);
		    if(ptCell->getOccupancy()<=0) ptCell->hasGaussian_ = false; 
		}
	    }
	}
	if(updatePositive) {
	    int i=0,j=0;
	    if(it->x > centerX+sizeXmeters/2-resolution/2) {
		i = 1;
	    }
	    if(it->x < centerX-sizeXmeters/2-resolution/2) {
		i = -1;
	    }
	    if(it->y > centerY+sizeYmeters/2-resolution/2) {
		j = 1;
	    }
	    if(it->y < centerY-sizeYmeters/2-resolution/2) {
		j = -1;
	    }
	    if(grid_[i+1][j+1]->isInside(*it)) {
		ptCell = (grid_[i+1][j+1]->addPoint(*it));
		if(ptCell!=NULL) {
		    update_set.insert(ptCell);
		}
	    }
	}
	it++;
    }
    isFirstLoad_ = false;
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
void NDTMapHMT::addPointCloudMeanUpdate(const Eigen::Vector3d &origin, 
	const pcl::PointCloud<pcl::PointXYZ> &pc, 
	const Eigen::Vector3d &localmapsize,
	unsigned int maxnumpoints, float occupancy_limit,double maxz, double sensor_noise){

    //std::cout<<localmapsize.transpose()<<" "<<maxnumpoints<<" "<<occupancy_limit<<" "<<maxz<<" "<<sensor_noise<<std::endl;
    //fprintf(stderr,"Mean nasty point cloud hmt\n");
    if(isFirstLoad_){
	//fprintf(stderr,"First load... ");
	loadPointCloud(pc);
	return;
    }
    this->setInsertPosition(origin);

    LazyGrid *lz = grid_[1][1];
    ///fprintf(stderr,"UPDATING\n");
    double sizeXmeters, sizeYmeters, sizeZmeters;
    double eps = 0.001;
    lz->getGridSizeInMeters(sizeXmeters, sizeYmeters, sizeZmeters);
    double centerX,centerY,centerZ;
    lz->getCenter(centerX, centerY, centerZ);
    int sizeX,sizeY,sizeZ;
    lz->getGridSize(sizeX, sizeY, sizeZ);

    Eigen::Vector3d old_centroid;
    old_centroid(0) = centerX;
    old_centroid(1) = centerY;
    old_centroid(2) = centerZ;

    ///Lets first create a local ndmap (this really works only if we have lazy grid as well)
    lslgeneric::NDTMap ndlocal(new lslgeneric::LazyGrid(resolution));
    ndlocal.loadPointCloudCentroid(pc,origin, old_centroid, localmapsize, max_range_); 
    ndlocal.computeNDTCells();

    ///TODO: This returns now copies --- should be pointers?
    std::vector<lslgeneric::NDTCell*> ndts;
    ndts = ndlocal.getAllCells();

    NDTCell* ptCell = NULL; 
    std::vector< NDTCell*> cells; 

    pcl::PointXYZ pt;
    pcl::PointXYZ po,mpoint;
    po.x = origin(0); po.y = origin(1); po.z = origin(2);
    //fprintf(stderr,"NUM = %d\n",ndts.size());
    bool updatePositive = true;
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
	Eigen::Vector3d diff,diff2,pf;
	Eigen::Vector3d m = ndts[it]->getMean();
	mpoint.x = m(0);
	mpoint.y = m(1);
	mpoint.z = m(2);
	diff = m-origin;
	double l = diff.norm();

	if(l>max_range_)
	{
	    //fprintf(stderr,"Very long distance (%lf) :( \n",l);
	    it++;
	    continue;
	}
	cells.clear();
	if(!lz->traceLineWithEndpoint(origin,mpoint,diff,maxz,cells,pf)) {
	    it++;
	    continue;
	}
	//check if pf is within one resolution of diff
	diff2 = m-pf;
	if((diff2).norm() > eps) {
	    //line continues outside center map, need to trace more!
	    int i=0,j=0;
	    if(mpoint.x > centerX+sizeXmeters/2-resolution/2) {
		i = 1;
	    }
	    if(mpoint.x < centerX-sizeXmeters/2-resolution/2) {
		i = -1;
	    }
	    if(mpoint.y > centerY+sizeYmeters/2-resolution/2) {
		j = 1;
	    }
	    if(mpoint.y < centerY-sizeYmeters/2-resolution/2) {
		j = -1;
	    }
	    std::vector< NDTCell*> cells2; 
	    if(!grid_[i+1][j+1]->traceLineWithEndpoint(pf,mpoint,diff2,maxz,cells,pf)) {
		std::cerr<<"ray outside central map: at "<<mpoint.x<<" "<<mpoint.y<<" "<<mpoint.z<<" with diff "<<diff.transpose()<<" starts at " <<origin.transpose()
		    <<" stops at "<<pf.transpose()<<" diff2 is "<<diff2.transpose()<<" check grid "<<i<<" "<<j<<std::endl;
		it++;
		continue;
	    }
	    cells.insert(cells.begin(),cells2.begin(),cells2.end());
	    //TODO: we might need to make one more raytrace here
	}

	for(unsigned int i=0; i<cells.size(); i++)
	{
	    ptCell = cells[i];

	    if(ptCell != NULL){
		double l2target = 0;
		if(ptCell->hasGaussian_)
		{
		    Eigen::Vector3d out, pend,vpt;
		    pend = m;
		    double lik = ptCell->computeMaximumLikelihoodAlongLine(po, mpoint, out);
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
		    ptCell->updateOccupancy(numpoints*logoddlik, occupancy_limit);
		    if(ptCell->getOccupancy()<=0) ptCell->hasGaussian_ = false; 
		}
		else
		{
		    ptCell->updateOccupancy(-0.2*numpoints, occupancy_limit);
		    if(ptCell->getOccupancy()<=0) ptCell->hasGaussian_ = false; 
		}
	    }
	    else{
		std::cerr<<"PANIC!!\n";
		index_->addPoint(pt); ///Add fake point to initialize!
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

/**
 * Add a distribution to the map
 */
void NDTMapHMT::addDistributionToCell(const Eigen::Matrix3d &ucov, const Eigen::Vector3d &umean, unsigned int numpointsindistribution, 
	float r, float g,float b,  unsigned int maxnumpoints, float max_occupancy)
{
    pcl::PointXYZ pt;
    pt.x = umean[0];
    pt.y = umean[1];
    pt.z = umean[2];
    
    double centerX,centerY,centerZ;
    grid_[1][1]->getCenter(centerX, centerY, centerZ);
    double sizeXmeters, sizeYmeters, sizeZmeters;
    grid_[1][1]->getGridSizeInMeters(sizeXmeters, sizeYmeters, sizeZmeters);
    
    int i=0, j=0;
    if(pt.x > centerX+sizeXmeters/2-resolution/2) {
	i = 1;
    }
    if(pt.x < centerX-sizeXmeters/2-resolution/2) {
	i = -1;
    }
    if(pt.y > centerY+sizeYmeters/2-resolution/2) {
	j = 1;
    }
    if(pt.y < centerY-sizeYmeters/2-resolution/2) {
	j = -1;
    }

    if(grid_[i+1][j+1]->isInside(pt)) {
	NDTCell *ptCell = NULL; 
	grid_[i+1][j+1]->getNDTCellAt(pt,ptCell);
	if(ptCell != NULL)
	{
	    ptCell->updateSampleVariance(ucov, umean, numpointsindistribution, true, max_occupancy,maxnumpoints);
	    ptCell->setRGB(r,g,b);
	} 
    }
    
}

/** Helper function, computes the  NDTCells
 */
void NDTMapHMT::computeNDTCells(int cellupdatemode, unsigned int maxnumpoints, float occupancy_limit, Eigen::Vector3d origin, double sensor_noise)
{
    conflictPoints.clear();
    typename std::set<NDTCell*>::iterator it = update_set.begin();
    while (it != update_set.end())
    {
        NDTCell *cell = *it;
        if(cell!=NULL)
        {
            cell->computeGaussian(cellupdatemode,maxnumpoints, occupancy_limit, origin,sensor_noise);
            ///Process the conflict points
            if(cell->points_.size()>0)
            {
                for(unsigned int i=0; i<cell->points_.size(); i++) conflictPoints.push_back(cell->points_[i]);
                cell->points_.clear();
            }
        }
        it++;
    }
    update_set.clear();
}



/** output methods for saving the map in the jff format
 */
int NDTMapHMT::writeTo()
{

    if(my_directory == "" || !grids_init)
    {
	std::cout<<"provide directory name\n";
	return -1;
    }
    
    std::cerr<<"writing to "<<my_directory<<std::endl;
    char fname[500];
    std::string jffnames[3][3];
    bool already[3][3];
    double centx, centy, centz, sizeX, sizeY, sizeZ;

    for(int i=0; i<3; i++) {
	for(int j=0; j<3; j++) {
	    if(grid_[i][j]==NULL) return -1;
	    grid_[i][j]->getCenter(centx,centy,centz);
	    snprintf(fname,499,"lz_%05lf_%05lf_%05lf.jff",centx,centy,centz);
	    jffnames[i][j] = fname;
	    already[i][j] = false;
	    //std::cerr<<"fname "<<jffnames[i][j]<<std::endl;
	}
    }
    grid_[1][1]->getGridSizeInMeters(sizeX, sizeY, sizeZ);
    //open metadata file
    std::string meta = my_directory;
    meta += "/metadata.txt";

    //std::cerr<<"metadata file at "<<meta<<std::endl;
    FILE* meta_f = fopen(meta.c_str(),"a+");
    if(meta_f==0) return -1;
    //read in all metadata
    //each line: centerx centery centerz filename\n
    char *line = NULL;
    size_t len;
    bool first=true;
    size_t length = 0;
    double epsilon = 1e-5;
    //std::cerr<<"reading metadata file at "<<meta<<std::endl;
    //------------------------------------------------------------------------//
    //read in all metadata
    //NOTE: v 2.0 -> Added header:
    //VERSION X.x\n
    //SIZE X.x\n
    //each next line: centerx centery centerz filename\n
    if(getline(&line,&len,meta_f) >0 ) {
	char *tk = strtok(line," ");
	if(tk==NULL) return false;
	if (strncmp(tk,"VERSION",7) == 0) {
	    tk = strtok(NULL," ");
	    if(tk==NULL) return false;
	    if(strncmp(tk,"2.0",3) == 0) {
		if(!getline(&line,&len,meta_f) >0 ) {
		    return false;
		}
		tk = strtok(line," ");
		if(tk==NULL) return false;
		if(strncmp(tk,"SIZE",4) != 0) return false;
		tk = strtok(NULL," ");
		double sizeMeta = atof(tk);
		if(fabsf(sizeMeta - sizeX) > 0.01) {
		    std::cerr<<"cannot write map, different grid size used...\n";
		    char ndir[500];
		    snprintf(ndir,499,"%s_%5d",my_directory.c_str(),rand());
		    my_directory = ndir;
		    std::cerr<<"SWITCHING DIRECTORY! "<<my_directory<<std::endl;
		    int res = mkdir(my_directory.c_str(),S_IRWXU);
		    if(res <0 ) return false;
		    fclose(meta_f);
		    meta   = my_directory+"/metadata.txt";
		    meta_f = fopen(meta.c_str(),"a+");
		    fprintf(meta_f,"VERSION 2.0\nSIZE %lf\n",sizeX);
		}
	    }
	    //add cases ...
	    
	} else {
	    //this is a version 1.0 file, close and re-open and go on
	    std::cerr<<"metafile version 1.0, no protection against different grid size\n";
	    fclose(meta_f);
	    meta_f = fopen(meta.c_str(),"a+");
	}
    } else {
	//empty metadata file -> write in version and grid sizes!
	fprintf(meta_f,"VERSION 2.0\nSIZE %lf\n",sizeX);
    }
    //------------------------------------------------------------------------//
    while(getline(&line,&len,meta_f) >0 )
    {
	pcl::PointXYZ cen;
	char *token = strtok(line," ");
	if(token==NULL) return -1;
	cen.x = atof(token);
	token = strtok(NULL," ");
	if(token==NULL) return -1;
	cen.y = atof(token);
	token = strtok(NULL," ");
	if(token==NULL) return -1;
	cen.z = atof(token);
	
	token = strtok(NULL," ");
	if(token==NULL) return -1;
	
	for(int i=0; i<3; i++) {
	    for(int j=0; j<3; j++) {
		if(grid_[i][j]==NULL) return -1;
		grid_[i][j]->getCenter(centx,centy,centz);
		if(fabsf(cen.x-centx) < epsilon && 
		   fabsf(cen.y-centy) < epsilon && fabsf(cen.z-centz) < epsilon) {
		    already[i][j] = true;
		    token[strlen(token)-1]='\0';
		    jffnames[i][j] = token;
		    //std::cerr<<"file already for "<<jffnames[i][j]<<std::endl;
		}   
	    }
	}
		
    }
    //fclose(meta_f);

    //now open for appending, in case we need to add new metadata
    //meta_f = fopen(meta.c_str(),"a");
    //std::cerr<<"appending metadata file at "<<meta<<std::endl;
    for(int i=0; i<3; i++) {
	for(int j=0; j<3; j++) {
	    if(grid_[i][j]==NULL) return -1;
	    grid_[i][j]->getCenter(centx,centy,centz);
	    if(!already[i][j]) {
		fprintf(meta_f,"%05lf %05lf %05lf %s\n",centx,centy,centz,jffnames[i][j].c_str());
		//std::cerr<<"appending metadata for "<<jffnames[i][j]<<std::endl;
	    }
	}
    }
    fclose(meta_f);

    //now pass through all grids and write to disk
    std::string sep = "/";
    for(int i=0; i<3; i++) {
	for(int j=0; j<3; j++) {
	    std::string path = my_directory+sep+jffnames[i][j];
	    //std::cout<<"saving "<<path<<std::endl;
	    FILE * jffout = fopen(path.c_str(), "w+b");
	    fwrite(_JFFVERSION_, sizeof(char), strlen(_JFFVERSION_), jffout);
	    int indexType[1] = {3};
	    fwrite(indexType, sizeof(int), 1, jffout);
	    double sizeXmeters, sizeYmeters, sizeZmeters;
	    double cellSizeX, cellSizeY, cellSizeZ;
	    double centerX, centerY, centerZ;
	    LazyGrid *ind = grid_[i][j];
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
	    typename SpatialIndex::CellVectorItr it = ind->begin();
	    while (it != ind->end())
	    {
		NDTCell *cell = (*it);
		if(cell!=NULL)
		{
		    if(cell->writeToJFF(jffout) < 0)	return -1;
		}
		else
		{
		    // do nothing
		}
		it++;
	    }
	    fclose(jffout);
	}
    }

    return 0;
}

bool NDTMapHMT::tryLoad(const double &cx, const double &cy, const double &cz, LazyGrid *&grid) {
std::cout<<"trying to load at "<<cx<<" "<<cy<<" "<<cz<<std::endl;
    if(my_directory == "" || !grids_init)
    {
	std::cout<<"provide directory name\n";
	return false;
    }

    std::string jffname;
    std::string meta = my_directory;
    meta += "/metadata.txt";
    //std::cerr<<"metadata file at "<<meta<<std::endl;
    FILE* meta_f = fopen(meta.c_str(),"a+");
    if(meta_f==0) return -1;
    //read in all metadata
    //each line: centerx centery centerz filename\n
    char *line = NULL;
    size_t len;
    bool found=false;
    size_t length = 0;
    double epsilon = 1e-2;
    //std::cerr<<"reading metadata file at "<<meta<<std::endl;
    // ---------------------------------------------------------------------------- //
    double sizeX,sizeY,sizeZ;
    grid_[1][1]->getGridSizeInMeters(sizeX, sizeY, sizeZ);
    //read in all metadata
    //NOTE: v 2.0 -> Added header:
    //VERSION X.x\n
    //SIZE X.x\n
    //each next line: centerx centery centerz filename\n
    if(getline(&line,&len,meta_f) >0 ) {
	char *tk = strtok(line," ");
	if(tk==NULL) return false;
	if (strncmp(tk,"VERSION",7) == 0) {
	    tk = strtok(NULL," ");
	    if(tk==NULL) return false;
	    if(strncmp(tk,"2.0",3) == 0) {
		if(!getline(&line,&len,meta_f) >0 ) {
		    return false;
		}
		tk = strtok(line," ");
		if(tk==NULL) return false;
		if(strncmp(tk,"SIZE",4) != 0) return false;
		tk = strtok(NULL," ");
		double sizeMeta = atof(tk);
		if(fabsf(sizeMeta - sizeX) > 0.01) {
		    std::cerr<<"cannot load map, different grid size used... reverting to empty map\n";
		    return false;
		}
	    }
	    //add cases ...
	    
	} else {
	    //this is a version 1.0 file, close and re-open and go on
	    std::cerr<<"metafile version 1.0, no protection against different grid size\n";
	    fclose(meta_f);
	    meta_f = fopen(meta.c_str(),"a+");
	}
    }
    // ---------------------------------------------------------------------------- //

    while(getline(&line,&len,meta_f) >0 )
    {
	pcl::PointXYZ cen;
	char *token = strtok(line," ");
	if(token==NULL) return -1;
	cen.x = atof(token);
	token = strtok(NULL," ");
	if(token==NULL) return -1;
	cen.y = atof(token);
	token = strtok(NULL," ");
	if(token==NULL) return -1;
	cen.z = atof(token);
	
	token = strtok(NULL," ");
	if(token==NULL) return -1;
	
	if(fabsf(cen.x-cx) < epsilon && 
		fabsf(cen.y-cy) < epsilon && fabsf(cen.z-cz) < epsilon) {
	    token[strlen(token)-1]='\0';
	    jffname = token;
	    found = true;
	    //std::cerr<<"file already for "<<jffname<<std::endl;
	    break;
	}   
		
    }
    fclose(meta_f);
    if(!found) return false;
    
    FILE * jffin;
    std::string sep = "/";
    std::string path = my_directory+sep+jffname;
    std::cout<<"reading file "<<path<<std::endl;
    jffin = fopen(path.c_str(),"r+b");
    if(jffin==NULL) return false;

    char versionBuf[16];
    if(fread(&versionBuf, sizeof(char), strlen(_JFFVERSION_), jffin) <= 0)
    {
	std::cerr<<"reading version failed";
	return false;
    }
    versionBuf[strlen(_JFFVERSION_)] = '\0';
    int indexType;
    if(fread(&indexType, sizeof(int), 1, jffin) <= 0)
    {
	std::cerr<<"reading index type failed";
	return false;
    }

    NDTCell *ptCell = new NDTCell();
    grid = new LazyGrid(resolution);
    grid->setCellType(ptCell);
    delete ptCell;
    if(grid->loadFromJFF(jffin) < 0)
    {
	std::cerr<<"loading grid failed";
	return false;
    }
    fclose(jffin);

    return true;
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
int NDTMapHMT::loadFrom()
{

    if(my_directory == "" || !grids_init)
    {
	std::cout<<"provide directory name\n";
	return -1;
    }
   
    //read in metadata
    //find the grid centered at center and the grids around
    //read them in
   /* 

    if(filename == NULL)
    {
        JFFERR("problem outputing to jff");
    }

    FILE * jffin;
    jffin = fopen(filename,"r+b");
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
    if(gr->loadFromJFF(jffin) < 0)
    {
	JFFERR("Error loading LazyGrid");
    }
    NDTCell *ptCell = new NDTCell();
    index_->setCellType(ptCell);
    delete ptCell;
    fclose(jffin);

   // std::cout << "map loaded successfully " << versionBuf << std::endl;

    isFirstLoad_ = false;

    return 0;
    */
}

///Get the cell for which the point fall into (not the closest cell)
bool NDTMapHMT::getCellAtPoint(const pcl::PointXYZ &refPoint, NDTCell *&cell)
{
    if(grid_[1][1]->isInside(refPoint)) {
	grid_[1][1]->getNDTCellAt(refPoint,cell);
    } else {
	for(int i=0; i<3; ++i) {
	    for(int j=0; j<3; ++j) {
		if(grid_[i][j]->isInside(refPoint)) {
		    grid_[i][j]->getNDTCellAt(refPoint,cell);
		    break;
		}
	    }
    	}
    }
    return (cell != NULL);
}

//computes the *negative log likelihood* of a single observation
double NDTMapHMT::getLikelihoodForPoint(pcl::PointXYZ pt)
{
    double uniform=0.00100;
    NDTCell* ndCell = NULL;
    this->getCellAtPoint(pt,ndCell); 
    if(ndCell == NULL) return uniform;
    double prob = ndCell->getLikelihood(pt);
    prob = (prob<0) ? 0 : prob; //uniform!! TSV
    return prob;
}

std::vector<NDTCell*> NDTMapHMT::getInitializedCellsForPoint(const pcl::PointXYZ pt) const
{
    std::vector<NDTCell*> cells, tmpcells;
    for(int i=0; i<3; ++i) {
	for(int j=0; j<3; ++j) {
	    if(grid_[i][j]->isInside(pt)) {
		tmpcells = grid_[i][j]->getClosestCells(pt);
		cells.insert(cells.begin(),tmpcells.begin(),tmpcells.end());
	    }
	}
    }
    return cells;
}

std::vector<NDTCell*> NDTMapHMT::getCellsForPoint(const pcl::PointXYZ pt, int n_neigh, bool checkForGaussian) const
{
    std::vector<NDTCell*> cells, tmpcells;
    for(int i=0; i<3; ++i) {
	for(int j=0; j<3; ++j) {
	    if(grid_[i][j]->isInside(pt)) {
		tmpcells = grid_[i][j]->getClosestNDTCells(pt,n_neigh,checkForGaussian);
		cells.insert(cells.begin(),tmpcells.begin(),tmpcells.end());
	    }
	}
    }
    return cells;
}

bool NDTMapHMT::getCellForPoint(const pcl::PointXYZ &pt, NDTCell* &out_cell, bool checkForGaussian) const
{
    out_cell = NULL;
    if(grid_[1][1]->isInside(pt)) {
	out_cell = grid_[1][1]->getClosestNDTCell(pt,checkForGaussian);
	return true;
    } else {
	for(int i=0; i<3; ++i) {
	    for(int j=0; j<3; ++j) {
		if(grid_[i][j]->isInside(pt)) {
		    out_cell = grid_[i][j]->getClosestNDTCell(pt,checkForGaussian);
		    return true;
		}
	    }
    	}
    }
    return false;
}

std::vector<NDTCell*> NDTMapHMT::pseudoTransformNDT(Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> T)
{

    std::vector<NDTCell*> ret;
    for(int i=0; i<3; ++i) {
	for(int j=0; j<3; ++j) {
	    typename SpatialIndex::CellVectorItr it = grid_[i][j]->begin();
	    while (it != grid_[i][j]->end())
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
			ret.push_back(nd);
		    }
		}
		else
		{
		    //ERR("problem casting cell to NDT!\n");
		}
		it++;
	    }
	}
    }
    return ret;
}

std::vector<lslgeneric::NDTCell*> NDTMapHMT::getAllCells() const
{

    std::vector<NDTCell*> ret;
    for(int i=0; i<3; ++i) {
	for(int j=0; j<3; ++j) {
	    //if(i==1 && j==1) continue;
	    typename SpatialIndex::CellVectorItr it = grid_[i][j]->begin();
	    while (it != grid_[i][j]->end())
	    {
		NDTCell *cell = (*it);
		if(cell!=NULL)
		{
		    if(cell->hasGaussian_)
		    {
			NDTCell* nd = (NDTCell*)cell->copy();
			ret.push_back(nd);
		    }
		}
		else
		{
		}
		it++;
	    }
	}
    }
    return ret;
}

std::vector<lslgeneric::NDTCell*> NDTMapHMT::getAllInitializedCells()
{

    std::vector<NDTCell*> ret;
    for(int i=0; i<3; ++i) {
	for(int j=0; j<3; ++j) {
	    typename SpatialIndex::CellVectorItr it = grid_[i][j]->begin();
	    while (it != grid_[i][j]->end())
	    {
		NDTCell *cell = (*it);
		if(cell!=NULL)
		{
		    NDTCell* nd = (NDTCell*)cell->copy();
		    ret.push_back(nd);
		}
		else
		{
		}
		it++;
	    }
	}
    }
    return ret;
}

int NDTMapHMT::numberOfActiveCells()
{
    int ret = 0;
    for(int i=0; i<3; ++i) {
	for(int j=0; j<3; ++j) {
	    typename SpatialIndex::CellVectorItr it = grid_[i][j]->begin();
	    while (it != grid_[i][j]->end())
	    {
		NDTCell *cell = (*it);
		if(cell!=NULL)
		{
		    if(cell->hasGaussian_)
		    {
			ret++;
		    }
		}
		it++;
	    }
	}
    }
    return ret;
}



}
