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
#pragma once
#ifndef NDT_HISTOGRAM_HH
#define NDT_HISTOGRAM_HH

namespace perception_oru{

  /**
   * \short
   * Histogram of NDT cells, to be used as a global appearance
   * descriptor of a scene.
   *
   * This is a reimplementation of the code used in the following
   * paper.
   *
   * - Magnusson, M. , Andreasson, H. , Nüchter, A. & Lilienthal,
   *   A. J. (2009). Automatic appearance-based loop detection from
   *   three-dimensional laser data using the normal distributions
   *   transform. Journal of field robotics, 26 (11-12), 892-914.
   */
  class NDTHistogram{
  private:
      mutable std::mutex cerr_mutex;
      Eigen::MatrixXi histogramBinsFlat; ///< The flat (planar, according to the planarity threshold) histogram bins.
      Eigen::MatrixXi histogramBinsLine; ///< The linear (according to the linearity threshold) histogram bins.
      Eigen::MatrixXi histogramBinsSphere; ///< The histogram bins that are not "linear" or "flat".
    
    int N_LINE_BINS; ///< Number of linear bins to use.
    int N_FLAT_BINS; ///< Number of flat bins to use.
    int N_SPHERE_BINS; ///< Number of spherical bins to use.
    int N_CLASSES;
    double D1, D2; 
    bool inited; ///< Whether the histogram is initialised or not.

    std::vector< Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>,Eigen::aligned_allocator<Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> > > topThree; 
    double topThreeS[3]; 

    Eigen::MatrixXi dist_histogramBinsFlat[3];
    Eigen::MatrixXi dist_histogramBinsLine[3];
    Eigen::MatrixXi dist_histogramBinsSphere[3];

    std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > averageDirections;

    /** 
     * Build a histogram descriptor given an NDT map.
     * 
     * The cells in the NDT map are sorted into histogram bins
     * depending on the eigenvalues of the distribution in the cell
     * and the thresholds linear_factor and flat_factor. Each NDT cell
     * has 3 eigenvalues (in the 3D case),
     * sorted as follows: \f$e_1 \leq e_2 \leq e_3\f$.
     *
     * In Magnusson et al. (2009), \f$ l = f = 10 = 1/t_e\f$.
     * 
     * @param map The input data.
     * @param l Threshold for what is a "linear" distribution. If \f$e_3 > l \cdot e_2\f$ then it is classified as linear.
     * @param f Threshold for what is a "flat" distribution. If \f$e_2 > f \cdot e_1\f$, and it is not linear, then it is classified as flat.
     */
    void constructHistogram(NDTMap **map,
                            double l = 50,
                            double f = 50
                            );

    void incrementLineBin(double d, int c);
    void incrementFlatBin(Eigen::Vector3d &normal, double d, int c);
    void incrementSphereBin(double d, int c);

    void computeDirections();
    void closedFormSolution(pcl::PointCloud<pcl::PointXYZ> &src, pcl::PointCloud<pcl::PointXYZ> &trgt,Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> &T);

  public:
    /** 
     * Default constructor.
     * 
     * @param linear_classes 
     * @param flat_classes 
     * @param spherical_classes 
     */
    NDTHistogram (int linear_classes = 1,
                  int flat_classes = 40,
                  int spherical_classes = 10,
                  int n_classes=0,
		  double _D1 = 5.0,
		  double _D2 = 10.0
		 );
    void init (int linear_classes = 1,
                  int flat_classes = 40,
                  int spherical_classes = 10,
                  int n_classes=0,
		  double _D1 = 5.0,
		  double _D2 = 10.0
		 );

    /** 
     * Construct a histogram from an NDTMap. 
     * 
     * @param map The 3D data to generate a histogram from. This typically comes from a single 3D scan.
     * @param linear_classes 
     * @param flat_classes 
     * @param spherical_classes 
     */
    NDTHistogram (NDTMap **map,
                  int linear_classes = 1,
                  int flat_classes = 40,
                  int spherical_classes = 10,
                  int n_classes=0,
		  double _D1 = 5.0,
		  double _D2 = 10.0
		 );

    /** 
     * Copy constructor.
     */
    NDTHistogram (const NDTHistogram& other);
    double calculateEntropy();
    double ENTROPY;
    omp_lock_t writelock;
    
    /**
     * Get the transform that brings me close to target.
     */
    void bestFitToHistogram(NDTHistogram &target, Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> &T, bool bound_transform = true);

    /** 
     * Print out histogram data on stdout. 
     * 
     * @param bMatlab Whether to print in a format suitable for matlab
     * plotting.
     */
    void printHistogram(bool bMatlab=false);

    
    /**
     * Call this to get the 1/2/3 best option, _after_ a call to
     * bestFitToHistogram.
     */
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

    /**
     * Compute similarity measure w.r.t. another histogram.
     */
    double getSimilarity(NDTHistogram &other);

    double getSimilarity(NDTHistogram &other, Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> &T);

    std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> > directions;

  public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  };
}
#endif
