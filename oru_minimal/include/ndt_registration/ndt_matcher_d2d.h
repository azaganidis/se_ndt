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
#ifndef NDT_MATCHER_D2D_HH
#define NDT_MATCHER_D2D_HH

#include "ndt_map/ndt_map.h"
#include "pcl/point_cloud.h"
#include "Eigen/Core"

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

namespace perception_oru
{

/**
 * This class implements NDT registration for 3D point cloud scans.
 */
class NDTMatcherD2D
{
public:
    /**
     parametrized constructor. A default set is (false,false,true,empty_vector). parameters are:
    \param _isIrregularGrid --- experimental single pass through an irregular grid. also unstable
    \param useDefaultGridResolutions --- if set, the following parameter is set to a preset value
    \param _resolutions --- if previous bool is not set, these are the resolutions (in reverse order) that we will go through
    */
    NDTMatcherD2D(bool _isIrregularGrid,
                  bool useDefaultGridResolutions, std::vector<double> _resolutions)
    {
        this->init(_isIrregularGrid,useDefaultGridResolutions,_resolutions);
    }
    NDTMatcherD2D()
    {
        this->init(false,true,std::vector<double>());
    }
    NDTMatcherD2D(const NDTMatcherD2D& other)
    {
        this->init(other.isIrregularGrid,false,other.resolutions);
    }

    /**
     * Register two point clouds. This method builds an NDT
     * representation of the "target" point cloud and uses that for
     * registering the "source" point cloud.
     * \param  target
     *   Reference data. NDT structure is built for this point cloud.
     * \param  source
     *   The output transformation registers this point cloud to \c target.
     * \param  T
     *   This is an input/output parameter. The initial value of \c T
     *   gives the initial pose estimate of \c moving. When the
     *   algorithm terminates, \c T holds the registration result.
     */
    bool match( pcl::PointCloud<pcl::PointXYZ>& target,
                pcl::PointCloud<pcl::PointXYZ>& source,
                Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T,
                bool useInitialGuess = false);

    /**
     * Registers a point cloud to an NDT structure.
     * \param  target
     *   Reference data.
     * \param  source
     *   The output transformation registers this point cloud to \c source.
     * \param  T
     *   This is an input/output parameter. The initial value of \c T
     *   gives the initial pose estimate of \c moving. When the
     *   algorithm terminates, \c T holds the registration result.
     */
    bool match( NDTMap& target,
                NDTMap& source,
                Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T,
                bool useInitialGuess = false);

    /**
      * @brief computes the covariance of the match between target and source, at T.
      * note --- computes NDT distributions based on the resolution in res
      * result is returned in cov
      */
    bool covariance( pcl::PointCloud<pcl::PointXYZ>& target,
                     pcl::PointCloud<pcl::PointXYZ>& source,
                     Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T,
                     Eigen::MatrixXd &cov
                   );

    /**
      * computes the covariance of the match between moving and fixed, at T.
      * result is returned in cov
      */
    bool covariance( NDTMap& target,
                     NDTMap& source,
                     Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T,
                     Eigen::MatrixXd &cov
                   );

    //compute the score of a point cloud to an NDT
    virtual double scoreNDT(std::vector<NDTCell*> &source, NDTMap &target);

    virtual double scoreNDT_OM(NDTMap &source, NDTMap &target);


    virtual double scoreNDTPositive(std::vector<NDTCell*> &sourceNDT, NDTMap &targetNDT,
                                    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T);

    //compute the score gradient & hessian of a point cloud + transformation to an NDT
    // input: moving, fixed, tr, computeHessian
    //output: score_gradient, Hessian, returns: score!
    virtual double derivativesNDT(
        const std::vector<NDTCell*> &source,
        const NDTMap &target,
        Eigen::MatrixXd &score_gradient,
        Eigen::MatrixXd &Hessian,
        bool computeHessian
    );
    /*
        void generateScoreDebug(const char* out, pcl::PointCloud<pcl::PointXYZ>& fixed,
    			    pcl::PointCloud<pcl::PointXYZ>& moving);
    */
    double finalscore;
    double current_resolution;
    ///max iterations, set in constructor 
    int ITR_MAX;
    ///sets step control on/off. set to true in constructor
    bool step_control;
    ///the change in score after which we converge. Set to 1e-3 in constructor
    double DELTA_SCORE;

    //should we try to regularize the hessian or just give up?
    bool regularize;
    //how many neighbours to use in the objective
    int n_neighbours;

    //perform line search to find the best descent rate (Mohre&Thuente)
    //adapted from NOX
    double lineSearchMT(
        Eigen::Matrix<double,6,1> &increment,
        std::vector<NDTCell*> &source,
        NDTMap &target) ;


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


protected:

    Eigen::Matrix<double,3,6> Jest;
    Eigen::Matrix<double,18,6> Hest;
    Eigen::Matrix<double,3,18> Zest;
    Eigen::Matrix<double,18,18> ZHest;

    Eigen::Matrix<double,3,3> dRdx, dRdy, dRdz;
    Eigen::Matrix<double,3,3> dRdxdx, dRdxdy, dRdxdz, dRdydy, dRdydz, dRdzdz;
    //lf = likelihood function d1 and d2 from the paper
    double lfd1,lfd2;
    //bool useSimpleDerivatives;
    int iteration_counter_internal;

    bool isIrregularGrid;
    std::vector<double> resolutions;

    //initializes stuff;
    void init(bool _isIrregularGrid,
              bool useDefaultGridResolutions, std::vector<double> _resolutions);

    //iteratively update the score gradient and hessian
    virtual bool update_gradient_hessian(
        Eigen::MatrixXd &score_gradient,
        Eigen::MatrixXd &Hessian,
        const Eigen::Vector3d &m1,
        const Eigen::Matrix3d &C1,
        const double &likelihood,
        bool computeHessian);

    //pre-computes the derivative matrices Jest, Hest, Zest, ZHest
    void computeDerivatives(Eigen::Vector3d &m1, Eigen::Matrix3d C1, bool computeHessian=true);

    //iteratively update the score gradient and hessian
    virtual bool update_gradient_hessian_local(
        Eigen::MatrixXd &score_gradient,
        Eigen::MatrixXd &Hessian,
        const Eigen::Vector3d &m1,
        const Eigen::Matrix3d &C1,
        const double &likelihood,
        const Eigen::Matrix<double,3,6> &_Jest,
        const Eigen::Matrix<double,18,6> &_Hest,
        const Eigen::Matrix<double,3,18> &_Zest,
        const Eigen::Matrix<double,18,18> &_ZHest,
        bool computeHessian);

    //pre-computes the derivative matrices Jest, Hest, Zest, ZHest
    void computeDerivativesLocal(Eigen::Vector3d &m1, Eigen::Matrix3d C1,
                                 Eigen::Matrix<double,3,6> &_Jest,
                                 Eigen::Matrix<double,18,6> &_Hest,
                                 Eigen::Matrix<double,3,18> &_Zest,
                                 Eigen::Matrix<double,18,18> &_ZHest,
                                 bool computeHessian);

    //perform a subsampling depending on user choice
    int NUMBER_OF_POINTS;
    /*
    void gradient_numeric(
        std::vector<NDTCell*> &moving,
        NDTMap &fixed,
        Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> &transform,
        //Eigen::Vector3d &eulerAngles,
        Eigen::Matrix<double,6,1> &score_gradient
        );
    */

protected:
    //storage for pre-computed angular derivatives
    //Eigen::Vector3d jest13, jest23, jest04, jest14, jest24, jest05, jest15, jest25;
    //Eigen::Vector3d a2,a3, b2,b3, c2,c3, d1,d2,d3, e1,e2,e3, f1,f2,f3;
    int NUMBER_OF_ACTIVE_CELLS;
    double normalizeAngle(double a);

    //vars for gradient
    Eigen::Matrix<double,6,1> xtBJ, xtBZBx, Q;
    //vars for hessian
    Eigen::Matrix<double,6,6> JtBJ, xtBZBJ, xtBH, xtBZBZBx, xtBZhBx;
    Eigen::Matrix<double,1,3> TMP1, xtB;

public:
    int nb_match_calls;
    int nb_success_reg;


public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} // end namespace

//#include <ndt_registration/impl/ndt_matcher_d2d.hpp>
#endif
