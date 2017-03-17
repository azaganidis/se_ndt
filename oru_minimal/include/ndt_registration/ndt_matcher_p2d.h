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
#ifndef NDT_MATCHER_HH
#define NDT_MATCHER_HH

#include "ndt_map/ndt_map.h"
#include "ndt_map/pointcloud_utils.h"
#include "pcl/point_cloud.h"
#include "Eigen/Core"

namespace lslgeneric
{

/**
 * This class implements NDT registration for 3D point cloud scans.
 */
class NDTMatcherP2D
{
public:
    NDTMatcherP2D(std::vector<double> _resolutions)
    {
        this->init(false,_resolutions);
    }
    NDTMatcherP2D()
    {
        this->init(true,std::vector<double>());
    }
    NDTMatcherP2D(const NDTMatcherP2D& other)
    {
        this->init(false,other.resolutions);
    }

    /**
     * Register two point clouds. This method builds an NDT
     * representation of the "fixed" point cloud and uses that for
     * registering the "moving" point cloud.
     * \param  fixed
     *   Reference data. NDT structure is built for this point cloud.
     * \param  moving
     *   The output transformation registers this point cloud to \c fixed.
     * \param  T
     *   This is an input/output parameter. The initial value of \c T
     *   gives the initial pose estimate of \c moving. When the
     *   algorithm terminates, \c T holds the registration result.
     * \return
     *   False if registration terminated after the max number of
     *   iterations (before convergence), true if registration
     *   terminated because the convergence criteria were fulfilled.
     */
    bool match( pcl::PointCloud<pcl::PointXYZ>& target,
                pcl::PointCloud<pcl::PointXYZ>& source,
                Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T );


  /** 
   * Check whether two point clouds are aligned, using the NDT score and Hessian.
   * 
   * \param  fixed
   *   Reference data. NDT structure is built for this point cloud (if needed).
   * \param  moving
   *   \c T will be applied to this point cloud.
   * \param T 
   *   Optional transformation that can be applied to \c moving. Defaults to zero.
   * \return 
   */
  void check( pcl::PointCloud<pcl::PointXYZ>& fixed,
              pcl::PointCloud<pcl::PointXYZ>& moving,
              Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T
              //= Eigen::Translation3f( 0.0, 0.0, 0.0 ) * Eigen::AngleAxisf( 0.0, Eigen::Vector3d( 1.0, 0.0, 0.0 ) )
              );
  
  
    /**
     * Registers a point cloud to an NDT structure.
     * \param  fixed
     *   Reference data.
     * \param  moving
     *   The output transformation registers this point cloud to \c fixed.
     * \param  T
     *   This is an input/output parameter. The initial value of \c T
     *   gives the initial pose estimate of \c moving. When the
     *   algorithm terminates, \c T holds the registration result.
     * \return
     *   False if registration terminated after the max number of
     *   iterations (before convergence), true if registration
     *   terminated because the convergence criteria were fulfilled.
     */
    bool match( NDTMap& target,
                pcl::PointCloud<pcl::PointXYZ>& source,
                Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T );



  
    /**
      * computes the covariance of the match between moving and fixed, at T.
      * note --- computes NDT distributions based on the resolution in res
      * result is returned in cov
      */
    bool covariance( pcl::PointCloud<pcl::PointXYZ>& target,
                     pcl::PointCloud<pcl::PointXYZ>& source,
                     Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T,
                     Eigen::Matrix<double,6,6> &cov
                   );

    //compute the score of a point cloud to an NDT
    double scorePointCloud(pcl::PointCloud<pcl::PointXYZ> &source,
                           NDTMap &target);

    //compute the score gradient & hessian of a point cloud + transformation to an NDT
    // input: moving, fixed, tr, computeHessian
    //output: score_gradient, Hessian
    void derivativesPointCloud(pcl::PointCloud<pcl::PointXYZ> &source,
                               NDTMap &target,
                               Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> &transform,
                               Eigen::Matrix<double,6,1> &score_gradient,
                               Eigen::Matrix<double,6,6> &Hessian,
                               bool computeHessian);

    void generateScoreDebug(const char* out, pcl::PointCloud<pcl::PointXYZ>& target,
                            pcl::PointCloud<pcl::PointXYZ>& source);

    double finalscore;
private:

    Eigen::Matrix<double,3,6> Jest;
    Eigen::Matrix<double,18,6> Hest;
    //lf = likelihood function d1 and d2 from the paper
    double lfd1,lfd2;
    bool useSimpleDerivatives;
    double current_resolution;
    bool isIrregularGrid;

    //pre-computes the multipliers of the derivatives for all points
    void precomputeAngleDerivatives(Eigen::Vector3d &eulerAngles);

    //iteratively update the score gradient
    bool update_score_gradient(Eigen::Matrix<double,6,1> &score_gradient,
                               Eigen::Vector3d &transformed,
                               Eigen::Matrix3d &Cinv);
    //iteratively update the hessian matrix
    void update_hessian(Eigen::Matrix<double,6,6> &Hessian,
                        Eigen::Vector3d &transformed,
                        Eigen::Matrix3d &Cinv);

    //pre-computes the derivative matrices Jest and Hest
    void computeDerivatives(pcl::PointXYZ &pt);

    //perform line search to find the best descent rate (naive case)
    double lineSearch(double score_here,
                      Eigen::Matrix<double,6,1> &score_gradient,
                      Eigen::Matrix<double,6,1> &increment,
                      pcl::PointCloud<pcl::PointXYZ> &source,
                      NDTMap &target) ;

    //perform line search to find the best descent rate (Mohre&Thuente)
    //adapted from NOX
    double lineSearchMT( Eigen::Matrix<double,6,1> &score_gradient_init,
                         Eigen::Matrix<double,6,1> &increment,
                         pcl::PointCloud<pcl::PointXYZ> &source,
                         Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> &globalT,
                         NDTMap &target) ;

    //compute finited difference derivative, for debug
    /*
        void fdd(pcl::PointCloud<PointSource> &source,
    	     NDTMap<PointTarget> &target,
    	     Eigen::Matrix<double,6,1> &score_gradient);
    */

    //auxiliary functions for MoreThuente line search
    struct MoreThuente
    {
        static double min(double a, double b);
        static double max(double a, double b);
        static double absmax(double a, double b, double c);
        static int cstep(double& stx, double& fx, double& dx,
                         double& sty, double& fy, double& dy,
                         double& stp, double& fp, double& dp,
                         bool& brackt, double stmin, double stmax);
    }; //end MoreThuente

    //perform a subsampling depending on user choice
    pcl::PointCloud<pcl::PointXYZ> subsample(pcl::PointCloud<pcl::PointXYZ>& original);
    int NUMBER_OF_POINTS;
    int NUMBER_OF_ACTIVE_CELLS;

private:
    //storage for pre-computed angular derivatives
    Eigen::Vector3d jest13, jest23, jest04, jest14, jest24, jest05, jest15, jest25;
    Eigen::Vector3d a2,a3, b2,b3, c2,c3, d1,d2,d3, e1,e2,e3, f1,f2,f3;

    std::vector<double> resolutions;
    //initializes stuff;
    void init(bool useDefaultGridResolutions, std::vector<double> _resolutions);
    double normalizeAngle(double a);
public:
    int ITR_MAX;
    double subsample_size;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} // end namespace

#endif
