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

/**
* This is an extension of NDTMap, which allows the update of
* occupancy on a cell level as well.
*
*/
#ifndef NDT_OCCUPANCY_MAP_HH
#define NDT_OCCUPANCY_MAP_HH

#include <ndt_map/spatial_index.h>
#include <ndt_map/ndt_cell.h>
#include <ndt_map/depth_camera.h>

#include <cstdlib>

#include <cv.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>

namespace lslgeneric
{

/** \brief Implements an NDT based spatial index
    \details This is an interface to a SpatialIndex (custom defined)
	that contains NDT cells. Provides methods to create from a PointCloud
*/
template <typename PointT>
class NDTOccupancyMap
{
public:

    NDTOccupancyMap()
    {
        index_ = NULL;
    }
    /** default constructor. The SpatialIndex sent as a paramter
    is used as a factory every time that loadPointCloud is called.
    it can/should be deallocated outside the class after the destruction of the NDTMap
    */
    NDTOccupancyMap(SpatialIndex<PointT> *idx, float _resolution)
    {
        //std::cout<<"STORE INDEX PROTOTYPE\n";
        index_ = idx;
        resolution = _resolution;
        //this is used to prevent memory de-allocation of the *si
        //si was allocated outside the NDT class and should be deallocated outside
        isFirstLoad_=true;
    }

    NDTOccupancyMap(const NDTOccupancyMap& other)
    {
        //std::cout<<"COPY MAP\n";
        if(other.index_ != NULL)
        {
            //std::cout<<"COPY INDEX\n";
            this->index_ = index_->copy();
            isFirstLoad_ = false;
        }
        this->resolution = other.resolution;
    }

    virtual ~NDTOccupancyMap()
    {
        //std::cout<<"DELETE MAP\n";
        if(index_ !=NULL && !isFirstLoad_)
        {
            //std::cout<<"DELETE INDEX\n";
            delete index_;
        }
    }

    /**
    	* Add new pointcloud to map
    	* @param &origin is the position of the sensor, from where the scan has been taken from.
    	* @param &pc is the pointcloud to be added
    	*/
    void addPointCloud(const Eigen::Vector3d &origin, const pcl::PointCloud<PointT> &pc);



    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /**
    * These load methods are not currently supported by the NDTOccupancyMap!
    *
    */
    /// each entry in the indices vector contains a set of indices to a NDC cell.
    void loadPointCloud(const pcl::PointCloud<PointT> &pc, const std::vector<std::vector<size_t> > &indices);
    void loadDepthImage(const cv::Mat& depthImage, DepthCamera<PointT> &cameraParams);
    pcl::PointCloud<PointT> loadDepthImageFeatures(const cv::Mat& depthImage, std::vector<cv::KeyPoint> &keypoints,
            size_t &supportSize, double maxVar, DepthCamera<PointT> &cameraParams, bool estimateParamsDI=false, bool nonMean = false);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void computeNDTCells(int cellupdatemode = CELL_UPDATE_MODE_SAMPLE_VARIANCE_WITH_RESET);

    virtual void writeToVRML(FILE* fout);
    virtual void writeToVRML(FILE* fout, Eigen::Vector3d col);

    inline SpatialIndex<PointT>* getMyIndex() const
    {
        return index_;
    }
    /// return the spatial index used as a string
    std::string getMyIndexStr() const;

    //computes the likelihood of a single observation
    double getLikelihoodForPoint(PointT pt);

    //returns the covariance matrix of the closest cell to refPoint
    bool getCellForPoint(const PointT &refPoint, NDTCell<PointT> *&cell);
    std::vector<NDTCell<PointT>*> getCellsForPoint(const PointT pt, double radius);

    ///return the cell using a specific index (not available for all spatialindexes), will return NULL if the idx is not valid.
    NDTCell<PointT>* getCellIdx(unsigned int idx);

    std::vector<NDTCell<PointT>*> pseudoTransformNDT(Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> T);

    /**
    * Returns all computed cells from the map
    */
    std::vector<lslgeneric::NDTCell<PointT>*> getAllCells();

    /**
    * Returns all computed cells from the map
    */
    std::vector<lslgeneric::NDTCell<PointT>*> getDynamicCells(unsigned int Timescale, float threshold);


    int numberOfActiveCells();

    //tsv: temporary debug function
    void debugToVRML(const char* fname, pcl::PointCloud<PointT> &pc);
protected:
    SpatialIndex<PointT> *index_;
    bool isFirstLoad_;
    void loadPointCloud(const Eigen::Vector3d &origin, const pcl::PointCloud<PointT> &pc);
    float resolution;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW


};

} // end namespace

#include <ndt_map/impl/ndt_occupancy_map.hpp>

#endif
