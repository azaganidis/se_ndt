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
#ifndef NDT_HISTOGRAM_HH
#define NDT_HISTOGRAM_HH

#include <ndt_map/ndt_map.h>
#include <vector>

namespace lslgeneric{
  class NDTHistogram{
  private:
    std::vector<int> histogramBinsFlat;
    std::vector<int> histogramBinsLine;
    std::vector<int> histogramBinsSphere;
    
    int N_LINE_BINS;
    int N_FLAT_BINS;
    int N_SPHERE_BINS;
    double D1, D2;
    bool inited;

    std::vector< Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>,Eigen::aligned_allocator<Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> > > topThree;
    double topThreeS[3];

    std::vector<int> dist_histogramBinsFlat[3];
    std::vector<int> dist_histogramBinsLine[3];
    std::vector<int> dist_histogramBinsSphere[3];

    std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > averageDirections;
    void constructHistogram(NDTMap &map);
    void incrementLineBin(double d);
    void incrementFlatBin(Eigen::Vector3d &normal, double d);
    void incrementSphereBin(double d);

    void computeDirections();
    void closedFormSolution(pcl::PointCloud<pcl::PointXYZ> &src, pcl::PointCloud<pcl::PointXYZ> &trgt,Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> &T);
  public:
    NDTHistogram();
    NDTHistogram (NDTMap &map);
    NDTHistogram (const NDTHistogram& other);
    
    //get the transform that brings me close to target
    void bestFitToHistogram(NDTHistogram &target, Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> &T, bool bound_transform = true);
    void printHistogram(bool bMatlab=false);

    //call this to get the 1/2/3 best option, AFTER a call to bestFitToHistogram
    double getTransform(size_t FIT_NUMBER, Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> &T){
      double ret = -1;
      T.setIdentity();
      if(FIT_NUMBER >=0 && FIT_NUMBER<3){
        T = topThree[FIT_NUMBER];
        ret = topThreeS[FIT_NUMBER];
      }
      return ret;
    }

    pcl::PointCloud<pcl::PointXYZ> getDominantDirections(int nDirections);
    double getSimilarity(NDTHistogram &other);
    double getSimilarity(NDTHistogram &other, Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> &T);

    std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > directions;
  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      };
}
#endif
