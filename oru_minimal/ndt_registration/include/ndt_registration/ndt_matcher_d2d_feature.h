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
#ifndef NDTMATCHERFEATUREF2F_HH
#define NDTMATCHERFEATUREF2F_HH

#include <ndt_registration/ndt_matcher_d2d.h>
namespace lslgeneric
{
/**
 * This class implements NDT / NDT registration with a priory known correspondances.
 */
class NDTMatcherFeatureD2D : public lslgeneric::NDTMatcherD2D
{
public:
    NDTMatcherFeatureD2D(const std::vector<std::pair<int, int> > &corr, double trimFactor = 1.) : _corr(corr), _trimFactor(trimFactor)
    {
        _goodCorr.resize(corr.size());
        std::fill(_goodCorr.begin(), _goodCorr.end(), true);
    }

    /**
     * Registers a point cloud to an NDT structure.
     * \param  target
     *   Reference data.
     * \param  source
     *   The output transformation registers this point cloud to \c target.
     * \param  T
     *   This is an input/output parameter. The initial value of \c T
     *   gives the initial pose estimate of \c source. When the
     *   algorithm terminates, \c T holds the registration result.
    bool match( lslgeneric::NDTMap<PointTarget>& target,
          lslgeneric::NDTMap<PointSource>& source,
          Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T);
     */
    /**
     * computes the covariance of the match between moving and fixed, at T.
     * result is returned in cov
     */
    bool covariance( lslgeneric::NDTMap& target,
                     lslgeneric::NDTMap& source,
                     Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T,
                     Eigen::Matrix<double,6,6> &cov
                   );

    //compute the score of a point cloud to an NDT
    virtual double scoreNDT(std::vector<lslgeneric::NDTCell*> &source,
                            lslgeneric::NDTMap &target,
                            Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T);

    //compute the score gradient & hessian of a point cloud + transformation to an NDT
    // input: moving, fixed, tr, computeHessian
    //output: score_gradient, Hessian
    virtual double derivativesNDT(
	    const std::vector<NDTCell*> &sourceNDT,
	    const NDTMap &targetNDT,
	    Eigen::MatrixXd &score_gradient,
	    Eigen::MatrixXd &Hessian,
	    bool computeHessian
    );
#if 0
    virtual bool update_gradient_hessian(
        Eigen::Matrix<double,6,1> &score_gradient,
        Eigen::Matrix<double,6,6> &Hessian,

        Eigen::Vector3d &m1,
        Eigen::Matrix3d &C1);
#endif

    using NDTMatcherD2D::Jest;
    using NDTMatcherD2D::Hest;
    using NDTMatcherD2D::Zest;
    using NDTMatcherD2D::ZHest;
    using NDTMatcherD2D::lfd1;
    using NDTMatcherD2D::lfd2;
    using NDTMatcherD2D::normalizeAngle;
    using NDTMatcherD2D::NUMBER_OF_ACTIVE_CELLS;
protected:
    const std::vector<std::pair<int, int> > & _corr;
    double _trimFactor;
    std::vector<bool> _goodCorr;
};
} // namespace


#endif
