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
#ifndef DEPTH_CAMERA
#define DEPTH_CAMERA

//#include <cv.h>
#include <opencv/cv.hpp>
#include <pcl/point_cloud.h>
#include <iostream>

namespace lslgeneric
{

template <typename PointT>
class DepthCamera
{
private:
    cv::Mat _camMat, _dist;
    cv::Mat _lookupTable;

    inline cv::Mat getCameraMatrix(double fx, double fy, double cx, double cy)
    {
        cv::Mat ret = cv::Mat::zeros(3,3, CV_64F);
        ret.at<double>(0,0) = fx;
        ret.at<double>(1,1) = fy;
        ret.at<double>(0,2) = cx;
        ret.at<double>(1,2) = cy;
        ret.at<double>(2,2) = 1.;
        return ret;
    }

    inline cv::Mat getDistVector(double d0, double d1, double d2, double d3, double d4)
    {
        cv::Mat ret = cv::Mat(5,1, CV_64F);
        ret.at<double>(0) = d0;
        ret.at<double>(1) = d1;
        ret.at<double>(2) = d2;
        ret.at<double>(3) = d3;
        ret.at<double>(4) = d4;
        return ret;
    }


public:

    double fx, fy, cx, cy, ds, scale_;
    std::vector<double> dist;
    bool isFloatImg;

    DepthCamera() { }
    DepthCamera (double &_fx, double &_fy, double &_cx, double &_cy, std::vector<double> &_distances, double _ds, bool _isFloatImg = false):
        fx(_fx),fy(_fy),cx(_cx),cy(_cy),ds(_ds),isFloatImg(_isFloatImg)
    {
        dist = _distances;
        _camMat = getCameraMatrix(fx, fy, cx, cy);
        _dist = getDistVector(dist[0], dist[1], dist[2], dist[3], dist[4]);
    }
    DepthCamera (const DepthCamera &other)
    {
        fx = other.fx;
        fy = other.fy;
        cx = other.cx;
        cy = other.cy;
        ds = other.ds;
        dist = other.dist;
        isFloatImg = other.isFloatImg;
        _camMat = getCameraMatrix(fx, fy, cx, cy);
        _dist = getDistVector(dist[0], dist[1], dist[2], dist[3], dist[4]);
        _lookupTable = other._lookupTable;
    }

    inline void setLookupTable(cv::Mat lookup)
    {
        _lookupTable = lookup;
    }

    inline void setupDepthPointCloudLookUpTable(const cv::Size &size) //, const cv::Mat &camMat, const cv::Mat &distVec, const double &dsFactor)
    {

        cv::Mat pixels = cv::Mat(size.height * size.width,1, CV_64FC2);
        // Fill the tmp values simply with the image coordinates
        {
            cv::Mat_<cv::Vec2d> _I = pixels;
            size_t iter = 0;
            for (int y = 0; y < size.height; y++)
            {
                for (int x = 0; x < size.width; x++)
                {
                    _I(iter)[0] = x;
                    _I(iter)[1] = y;
                    iter++;
                }
            }
        }
        cv::Mat normpixels = cv::Mat(pixels.size(), CV_64FC2); // normalized undistorted pixels
        cv::undistortPoints(pixels, normpixels, _camMat, _dist);

        _lookupTable = cv::Mat(normpixels.size(), CV_64FC3); // "normpixelsxyz"
        {
            cv::Mat_<cv::Vec2d> _J = normpixels;
            cv::Mat_<cv::Vec3d> _I = _lookupTable;
            size_t iter = 0;
            for (int y = 0; y < size.height; y++)
            {
                for (int x = 0; x < size.width; x++)
                {
                    _I(iter)[0] = _J(iter)[0]*ds;
                    _I(iter)[1] = _J(iter)[1]*ds;
                    _I(iter)[2] = ds;
                    iter++;
                }
            }
        }
    }

    inline size_t convertDepthImageToPointCloud(const cv::Mat &depthImg, pcl::PointCloud<PointT> &pc)
    {
        if (depthImg.depth() != CV_16U && !(isFloatImg) )
        {
            std::cerr<<"wrong depth image format - expected raw 16bit data\n";
            return 0;
        }
        if (depthImg.depth() != CV_32F && isFloatImg )
        {
            std::cerr<<"wrong depth image format - expected 32 bit floats\n";
            return 0;
        }
        float nan = std::numeric_limits<float>::quiet_NaN();

        size_t width = depthImg.size().width;
        size_t height = depthImg.size().height;
        size_t size = width*height;
        //std::cout<<"depth image "<<width<<"x"<<height<<" = "<<size<<std::endl;
        if (pc.size() != size || pc.width != width || pc.height != height || pc.is_dense != true)
        {
            pc.resize(size);
            pc.is_dense = true;
            pc.width = width;
            pc.height = height;
        }
        if(_lookupTable.empty())
        {
            //std::cout<<"setting up lookup table\n";
            this->setupDepthPointCloudLookUpTable(depthImg.size());
        }
        cv::Mat_<cv::Vec3d> _I = _lookupTable;
        double depth;
        if(!isFloatImg)
        {
            const unsigned short* pd = depthImg.ptr<unsigned short>(0);
            for (size_t i = 0; i < size; i++)
            {
                if (*pd == 0)
                {
                    //pc[i] = PointT(nan,nan,nan);
                    pc[i].x = nan;
                    pc[i].y = nan;
                    pc[i].z = nan;
                }
                else
                {
                    depth = *pd;
                    pc[i].x = depth * _I(i)[0];
                    pc[i].y = depth * _I(i)[1];
                    pc[i].z = depth * _I(i)[2];
//			      depth * _I(i)[1],
//			      depth * _I(i)[2]);
                    //	  std::cout<<"depth "<<depth<<" -> "<<pc[i].x<<" "<<pc[i].y<<" "<<pc[i].z<<std::endl;
                }
                pd++;
            }
        }
        else
        {
            const float* pd = depthImg.ptr<float>(0);
            for (size_t i = 0; i < size; i++)
            {
                if (*pd == 0)
                {
                    pc[i] = PointT(nan,nan,nan);
                }
                else
                {
                    depth = *pd;
                    pc[i] = PointT(depth * _I(i)[0],
                                   depth * _I(i)[1],
                                   depth * _I(i)[2]);
                    //	  std::cout<<"depth "<<depth<<" -> "<<pc[i].x<<" "<<pc[i].y<<" "<<pc[i].z<<std::endl;
                }
                pd++;
            }
        }
        return size;
    }

    inline size_t computePointsAtIndex(const cv::Mat &depthImg, cv::KeyPoint &keyPointCenter, size_t &support_size, pcl::PointCloud<PointT> &pc, PointT &center)
    {
        if (depthImg.depth() != CV_16U && !(isFloatImg) )
        {
            std::cerr<<"wrong depth image format - expected raw 16bit data\n";
            return 0;
        }
        if (depthImg.depth() != CV_32F && isFloatImg )
        {
            std::cerr<<"wrong depth image format - expected 32 bit floats\n";
            return 0;
        }
        float nan = std::numeric_limits<float>::quiet_NaN();

        size_t halfSize = support_size/2;
        int halfSizeI = (int)halfSize;

        size_t width = 2*halfSize+1;
        size_t height = 2*halfSize+1;
        size_t size = width*height;
        //std::cout<<"keypoint region "<<width<<"x"<<height<<" = "<<size<<std::endl;
        if (pc.size() != size || pc.width != width || pc.height != height || pc.is_dense != true)
        {
            pc.resize(size);
            pc.is_dense = true;
            pc.width = width;
            pc.height = height;
        }

        if(_lookupTable.empty())
        {
            //std::cout<<"setting up lookup table\n";
            this->setupDepthPointCloudLookUpTable(depthImg.size());
        }

        cv::Mat_<cv::Vec3d> _I = _lookupTable;

        int uKey = static_cast<int>(keyPointCenter.pt.x+0.5);
        int vKey = static_cast<int>(keyPointCenter.pt.y+0.5);
        //int indexKey = vKey * depthImg.size().width + uKey;

        int index, index_cloud, u, v;
        double depth;
        center = PointT(nan,nan,nan);


        if(!isFloatImg)
        {
            const unsigned short* pd = depthImg.ptr<unsigned short>(0);
            const unsigned short* pd_here;
            for (int i = -halfSizeI; i < halfSizeI+1; i++)
            {
                for (int j = -halfSizeI; j < halfSizeI+1; j++)
                {
                    v = vKey+j;
                    u = uKey+i;
                    if( u < 0 || v < 0 || u >= depthImg.size().width || v >= depthImg.size().height)
                    {
                        continue;
                    }
                    index = v * depthImg.size().width + u;
                    index_cloud = (j+halfSizeI) * width + (i+halfSizeI);
                    pd_here = pd+index;
                    if (*pd_here == 0)
                    {
                        pc[index_cloud] = PointT(nan,nan,nan);
                    }
                    else
                    {
                        depth = *pd_here;
                        pc[index_cloud] = PointT(depth * _I(index)[0],
                                                 depth * _I(index)[1],
                                                 depth * _I(index)[2]);
                        //	      std::cout<<"depth "<<depth<<" -> "<<pc[index_cloud].x<<" "<<pc[index_cloud].y<<" "<<pc[index_cloud].z<<std::endl;
                    }
                }
            }
        }
        else
        {
            const float* pd = depthImg.ptr<float>(0);
            const float* pd_here;
            for (int i = -halfSizeI; i < halfSizeI+1; i++)
            {
                for (int j = -halfSizeI; j < halfSizeI+1; j++)
                {
                    v = vKey+j;
                    u = uKey+i;
                    if( u < 0 || v < 0 || u >= depthImg.size().width || v >= depthImg.size().height)
                    {
                        continue;
                    }
                    index = v * depthImg.size().width + u;
                    index_cloud = (j+halfSizeI) * width + (i+halfSizeI);
                    pd_here = pd+index;

                    if (*pd_here == 0)
                    {
                        pc[index_cloud] = PointT(nan,nan,nan);
                    }
                    else
                    {

                        depth = *pd_here;
                        pc[index_cloud] = PointT(depth * _I(index)[0],
                                                 depth * _I(index)[1],
                                                 depth * _I(index)[2]);
                        //	      std::cout<<"depth "<<depth<<" -> "<<pc[index_cloud].x<<" "<<pc[index_cloud].y<<" "<<pc[index_cloud].z<<std::endl;
                    }

                }
            }
        }
        // Get the center (i == 0 && j = 0)
        center = pc[(halfSizeI) * width + (halfSizeI)];
        return size;
    }

    inline size_t computeParamsAtIndex(const cv::Mat &depthImg, cv::KeyPoint &keyPointCenter, size_t &support_size, Eigen::Vector3d &mean, Eigen::Matrix3d &cov)
    {
        if (depthImg.depth() != CV_16U && !(isFloatImg) )
        {
            std::cerr<<"wrong depth image format - expected raw 16bit data\n";
            return 0;
        }
        if (depthImg.depth() != CV_32F && isFloatImg )
        {
            std::cerr<<"wrong depth image format - expected 32 bit floats\n";
            return 0;
        }
        //float nan = std::numeric_limits<float>::quiet_NaN(); //?

        size_t halfSize = support_size/2;
        int halfSizeI = (int)halfSize;

        size_t width = 2*halfSize+1;
        size_t height = 2*halfSize+1;
        size_t size = width*height;
        //std::cout<<"keypoint region "<<width<<"x"<<height<<" = "<<size<<std::endl;

        if(_lookupTable.empty())
        {
            //std::cout<<"setting up lookup table\n";
            this->setupDepthPointCloudLookUpTable(depthImg.size());
        }
        cv::Mat_<cv::Vec3d> _I = _lookupTable;

        int uKey = static_cast<int>(keyPointCenter.pt.x+0.5);
        int vKey = static_cast<int>(keyPointCenter.pt.y+0.5);
        //int indexKey = vKey * depthImg.size().width + uKey;

        int index, u, v;
        double depth = 0;

        double depth_mean = 0, depth_var = 0;

        const unsigned short* pd = depthImg.ptr<unsigned short>(0);
        const unsigned short* pd_here;
        const float* pd_float = depthImg.ptr<float>(0);
        const float* pd_here_float;

        double depthsize = support_size*support_size;

        for (int q = 0; q<2; q++)
        {
            for (int i = -halfSizeI; i < halfSizeI+1; i++)
            {
                for (int j = -halfSizeI; j < halfSizeI+1; j++)
                {
                    v = vKey+j;
                    u = uKey+i;
                    if( u < 0 || v < 0 || u >= depthImg.size().width || v >= depthImg.size().height)
                    {
                        continue;
                    }
                    index = v * depthImg.size().width + u;
                    pd_here = pd+index;
                    pd_here_float = pd_float+index;

                    if(!isFloatImg)
                    {
                        if (*pd_here != 0)
                        {
                            depth = *pd_here;
                        }
                    }
                    else
                    {
                        if (*pd_here_float != 0)
                        {
                            depth = *pd_here_float;
                        }
                    }
                    if(q==0)
                    {
                        depth_mean += depth/depthsize;
                    }
                    else
                    {
                        depth_var += pow((depth-depth_mean),2)/depthsize;
                    }
                }
            }
        }
        double var = 0.001;

        index = vKey * depthImg.size().width + uKey;
        mean<<depth_mean*_I(index)[0],depth_mean*_I(index)[1],depth_mean*_I(index)[2];

        double x1 = depth_mean*_I(index)[0];
        double y1 = depth_mean*_I(index)[1];

        v = vKey+halfSize;
        u = uKey+halfSize;
        if( !(u < 0 || v < 0 || u >= depthImg.size().width || v >= depthImg.size().height) )
        {
            index = v * depthImg.size().width + u;
            pd_here = pd+index;
            pd_here_float = pd_float+index;
            if(!isFloatImg)
            {
                if (*pd_here != 0)
                {
                    depth = *pd_here;
                }
            }
            else
            {
                if (*pd_here_float != 0)
                {
                    depth = *pd_here_float;
                }
            }
            double x2 = depth*_I(index)[0];
            double y2 = depth*_I(index)[1];
            var = sqrt(pow(x1-x2,2)+pow(y1-y2,2))/4;
        }

        cov<<var,0,0,0,var,0,0,0,depth_var*ds*ds;

        Eigen::Vector3d eZ (0,0,1);
        double alpha = acos(eZ.dot(mean)/mean.norm());
        Eigen::AngleAxis<double> tr (-alpha, mean);

        //orient along beam direction
        cov = tr.inverse()*cov*tr;

        return size;
    }
};


};

#endif

