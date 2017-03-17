/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2010, AASS Research Center, Orebro University.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of AASS Research Center nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 */
#ifndef NDT_MAP_HMT_HH
#define NDT_MAP_HMT_HH

#include <ndt_map/ndt_map.h>
#include <ndt_map/spatial_index.h>
#include <ndt_map/ndt_cell.h>
//#include <ndt_map/depth_camera.h>
#include <ndt_map/lazy_grid.h>

#include <cstdlib>
#include <cv.h>
#include <sys/stat.h>
#include <sys/dir.h>

#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

namespace lslgeneric
{

/**
*  \brief Implements an NDT based spatial index on top of an hybrid metric topological grid
*  \author Jari Saarinen (jari.saarinen@aalto.fi) and Todor Stoyanov (todor.stoyanov@oru.se)
*  \version 0.1
*  \details This class derives from NDTMap and extends with functionalities for tracking several
* overlapping grids. On top of NDTMap, this class assumes updates are coming from a vehicle that moves
* through an environment and scans with a limited range sensor. Border areas are tracked as separate 
* LazyGrid maps and stored and loaded from disk when the robot switches the map region
*/
class NDTMapHMT : public NDTMap
{
public:

    /**
     * Construct with given centroid and sizes
     * @param cenx, ceny, cenz; (x,y,z) of the center of the map
     * @param sizex, sizey, sizez: The size of the map in each respective direction
     * NOTE: Implementation only for the lazy gridSpatialIndex<PointT> *idx
     **/
    NDTMapHMT(double resolution_, float cenx, float ceny, float cenz, float sizex, float sizey, float sizez, double max_range, std::string directory, bool _saveOnDelete=false)
    {
	max_range_ = max_range;
	grids_init=false;
	my_directory = directory;
	resolution = resolution_;

	//check if we can stat my_directory and try to create it if we can't
	DIR *mdir = opendir(my_directory.c_str());
	if(mdir == NULL) {
	    std::cerr<<"Error accessing map directory "<<my_directory<<" will attempt to recover:\n";
	    int res = mkdir(my_directory.c_str(),S_IRWXU);
	    if (res<0) {
		std::cerr<<"CRITICAL! FAILED to create directory\n";
		return;
	    }
	} else {
	    closedir(mdir);
	}

	LazyGrid *lz = new LazyGrid(resolution);
	this->index_ = lz;

	//this is used to prevent memory de-allocation of the *si
	//si was allocated outside the NDT class and should be deallocated outside
	this->isFirstLoad_=true;

	NDTCell *ptCell = new NDTCell();
	this->index_->setCellType(ptCell);
	delete ptCell;
	this->index_->setCenter(cenx,ceny,cenz);
	this->index_->setSize(sizex,sizey,sizez);
	this->map_sizex = sizex;
	this->map_sizey = sizey;
	this->map_sizez = sizez;
	this->is3D=true;
	lz->initializeAll();
	this->guess_size_ = false;
	this->saveOnDelete = _saveOnDelete;

	//create the grid of LazyGrids and set their centers
	initializeGrids();
    }
    NDTMapHMT(const NDTMapHMT& other)
    {
        if(other.index_ != NULL)
        {
            this->index_ = other.index_->copy();
            isFirstLoad_ = false;
	    this->max_range_ = other.max_range_;
	    this->grids_init = false;
	    this->saveOnDelete = other.saveOnDelete;
	    //TODO copy all LazyGrids
        }
    }
    /**
    * destructor
    	*/
    virtual ~NDTMapHMT()
    {
	if(saveOnDelete) {
	    this->writeTo();
	}
	for(int i=0; i<3; ++i) {
	    for(int j=0; j<3; ++j) {
		if(grid_[i][j]!=NULL) {
		    delete grid_[i][j];
		}
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
    virtual void loadPointCloud(const pcl::PointCloud<pcl::PointXYZ> &pc, double range_limit = -1);
    /**
    	* Add new pointcloud to map - This is the main interface for NDT-OM!
    	* Performs raytracing, updates conflicts and adds points to cells
    	* computeNDTCells must be called after calling this
    	*
    	* @param &origin is the position of the sensor, from where the scan has been taken from.
    	* @param &pc is the pointcloud to be added
    	* @param classifierTh A treshold to judge if the ray passes through a gaussian (obsolete)
    	* @param maxz threshold for the maximum z-coordinate value for the measurement point_cloud
    	* @param sensor_noise The expected standard deviation of the sensor noise
    	*/
    virtual void addPointCloud(const Eigen::Vector3d &origin, const pcl::PointCloud<pcl::PointXYZ> &pc, double classifierTh=0.06, 
				double maxz = 100.0, double sensor_noise = 0.25, double occupancy_limit = 255);

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
    virtual void addPointCloudMeanUpdate(const Eigen::Vector3d &origin, 
	    const pcl::PointCloud<pcl::PointXYZ> &pc, 
	    const Eigen::Vector3d &localmapsize,
	    unsigned int maxnumpoints = 1e9, float occupancy_limit=255 ,double maxz = 100.0, double sensor_noise = 0.25);


    /**
    * Adds a sample mean and covariance to the map
    * @param &ucov The covariance matrix to be added
    * @param &umean The mean of the normal distribution
    * @param numpointsindistribution The number of points used in computation of the sample mean and covariance
    * @param r,g,b -- optional color parameters
    * @param maxnumpoints -- optional adaptation of the gaussians
    */
    virtual void addDistributionToCell(const Eigen::Matrix3d &ucov,const Eigen::Vector3d &umean, unsigned int numpointsindistribution, 
															float r=0, float g=0,float b=0, unsigned int maxnumpoints=1e9, float max_occupancy=1024);

    /**
    * Computes the NDT-cells after a measurement has been added
    * @param cellupdatemode Defines the update mode (default CELL_UPDATE_MODE_SAMPLE_VARIANCE)
    * @param maxnumpoints Defines the forgetting factor (default 100000) the smaller the value the faster the adaptation
    */
    virtual void computeNDTCells(int cellupdatemode = CELL_UPDATE_MODE_SAMPLE_VARIANCE, unsigned int maxnumpoints = 1e9, float occupancy_limit=255, Eigen::Vector3d origin = Eigen::Vector3d(0,0,0), double sensor_noise=0.1);

    /**
    * Stuff for saving things to the directory of the map
    */
    int writeTo();
    int loadFrom();
    //cehck if there is a grid at this center and try to load it
    bool tryLoad(const double &cx, const double &cy, const double &cz, LazyGrid *&grid);

///------------------- non-essential stuff --------------///
    //computes the likelihood of a single observation
    virtual double getLikelihoodForPoint(pcl::PointXYZ pt);

    ///Get the cell for which the point fall into (not the closest cell)
    virtual bool getCellAtPoint(const pcl::PointXYZ &refPoint, NDTCell *&cell);

    virtual bool getCentroid(double &cx, double &cy, double &cz)
    {
        LazyGrid *lz = grid_[1][1];
        if(lz == NULL) return false;
        lz->getCenter(cx, cy, cz);
        return true;
    }
    /**
     * returns the closest cell to refPoint
     * Does not work with NDT-OM
     */
    virtual bool getCellForPoint(const pcl::PointXYZ &refPoint, NDTCell *&cell, bool checkForGaussian=true) const;
    /**
     * Returns all the cells within radius
     * Does not work with NDT-OM
     */
    virtual std::vector<NDTCell*> getCellsForPoint(const pcl::PointXYZ pt, int n_neighbours, bool checkForGaussian=true) const;
    /**
     * Returns all the cells within radius
     */
    virtual std::vector<NDTCell*> getInitializedCellsForPoint(const pcl::PointXYZ pt) const;

    /**
    * Returns a transformed NDT as a vector of NDT cells
    */
    virtual std::vector<NDTCell*> pseudoTransformNDT(Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> T);

    /**
     * Returns a transformed NDT as an NDT map with a CellVector data structure
    NDTMap<PointT>* pseudoTransformNDTMap(Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> T);
     */
    /**
     * Returns all computed cells from the map
     * This method gives all the vectors that contain a gaussian within a cell (hasGaussian is true).
     */
    virtual std::vector<lslgeneric::NDTCell*> getAllCells() const;
    /**
     * Returns all cells that have been initialized (including ones that do not contain gaussian at the moment).
     * This is useful if you want to use the empty cells or dynamic cells
     */
    virtual std::vector<lslgeneric::NDTCell*> getAllInitializedCells();

    int numberOfActiveCells();
    void setInsertPosition(const Eigen::Vector3d &newPos);
    bool tryLoadPosition(const Eigen::Vector3d &newPos);
protected:
    bool is3D;
    bool saveOnDelete;
    SpatialIndex *index_;
    bool isFirstLoad_;
    float map_sizex;
    float map_sizey;
    float map_sizez;
    float centerx,centery,centerz;
    bool grids_init;
    Eigen::Vector3d last_insert;
    std::string my_directory;
    double resolution;
    std::set<NDTCell*> update_set;

    double max_range_;
    lslgeneric::LazyGrid* grid_[3][3];
    //helper functions
    void initializeGrids();
    

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    pcl::PointCloud<pcl::PointXYZ> conflictPoints; ///< points that were conflicting during update
		

};

} // end namespace


#endif
