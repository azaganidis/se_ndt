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
#ifndef NDT_CELL_HH
#define NDT_CELL_HH
#define BIN_ARCHIVE
#include <ndt_map/impl/EventCounterData.hpp>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <vector>
#include <cstdio>
#include <fstream>

#include <Eigen/Eigen>
#include "boost/serialization/serialization.hpp"
#include <boost/serialization/array.hpp>
#include <boost/serialization/wrapper.hpp>
#ifdef BIN_ARCHIVE
#include <boost/archive/binary_oarchive.hpp>
#else
#include <boost/archive/text_oarchive.hpp>
#endif
/// A rather unsophisticated way of determining the
/// update method for a cell
/// Covariance intersection based on estimation
#define CELL_UPDATE_MODE_COVARIANCE_INTERSECTION 		0
/// Recursive Sample variance method [Chan, Gene, Randall, Updating Formulae and pairwise algorithm for computing sample variances, tech report Standford, 1979]
#define CELL_UPDATE_MODE_SAMPLE_VARIANCE 				1
/// Yguel, Vasquez, Aycard, Siegward, Laugier, Error-Driven Refinement of Multi-scale gaussian maps
#define CELL_UPDATE_MODE_ERROR_REFINEMENT				2
/// Estimate the surface (reduce the sensor noise)
#define CELL_UPDATE_MODE_SAMPLE_VARIANCE_SURFACE_ESTIMATION 3
///Student-T 
#define CELL_UPDATE_MODE_STUDENT_T 4

#define REFACTORED

#define JFFERR(x) std::cerr << x << std::endl; return -1;
#define _JFFVERSION_ "#JFF V0.50"

#define KEEP_POINTS
//#define ICRA_2013_NDT_OM_SIMPLE_MODE 
namespace perception_oru
{

/** \brief implements a normal distibution cell
    \details The base class for all NDT indeces, contains
mean, covariance matrix as well as eigen decomposition of covariance
*/

//template<typename PointT>
class NDTCell //: public Cell<PointT>
{
public:
  bool hasGaussian_;	///< indicates if the cell has a gaussian in it
  double cost; 	/// ndt_wavefront planner
  char isEmpty;	///<based on the most recent observation, is the cell seen empty (1), occupied (-1) or not at all (0)
  double consistency_score;
  unsigned int cloud_index;
  enum CellClass {HORIZONTAL=0, VERTICAL, INCLINED, ROUGH, UNKNOWN};
  std::vector<pcl::PointXYZ,Eigen::aligned_allocator<pcl::PointXYZ> > points_; ///The points falling into the cell - cleared after update
  void InitializeVariables(){
    cl_=UNKNOWN;
    hasGaussian_=false;
    R=G=B     =0;
    N         =0;
    occ       =0;   			///Occupancy value stored as "Log odds" (if you wish)
    emptyval  =0;			///The number of times a cell was observed empty (using ray casting)
    cost=INT_MAX;
    isEmpty=0;
    emptylik  =0;
    emptydist =0;
    max_occu_ =1;
    consistency_score=0;
    center_.x = center_.y = center_.z =0;
    xsize_    = ysize_    = zsize_    =0;
    cov_   =  Eigen::MatrixXd::Identity(3, 3);
    icov_  =  Eigen::MatrixXd::Identity(3, 3);
    evecs_ =  Eigen::MatrixXd::Identity(3, 3);
    mean_  =  Eigen::Vector3d(0,0,0);
    evals_ =  Eigen::Vector3d(0,0,0);
    d1_= d2_  =0;
    if(!parametersSet_)
      setParameters();

    ///Number of points used for Normal distribution estimation so far
    ///RGB values [0..1] - Special implementations for PointXYZRGB & PointXYZ

  }
  NDTCell()
  {
   InitializeVariables();
  }

  virtual ~NDTCell()
  {
    points_.clear();
  }

  NDTCell(pcl::PointXYZ &center, double &xsize, double &ysize, double &zsize)
  {
    InitializeVariables();
    center_ = center;
    xsize_ = xsize;
    ysize_ = ysize;
    zsize_ = zsize;

  }

  NDTCell(const NDTCell& other)
  {
    InitializeVariables();
    this->center_ = other.center_;
    this->xsize_ = other.xsize_;
    this->ysize_ = other.ysize_;
    this->zsize_ = other.zsize_;
    this->hasGaussian_ = other.hasGaussian_;
    this->R = other.R;
    this->G = other.G;
    this->B = other.B;
    this->N = other.N;
    this->occ = other.occ;
    this->max_occu_= other.getMaxOccupancy();
    this->emptyval = other.getEmptyval();
    this->edata = other.edata;
    this->consistency_score=other.consistency_score;
    this->isEmpty = other.isEmpty;
    this->emptydist=other.getEmptyDist();
    this->emptylik=other.getEmptyLik();
    this->setMean(other.getMean());
    this->setCov(other.getCov());
    this->icov_=other.getInverseCov();
    this->evals_=other.getEvals();
    this->evecs_=other.getEvecs();
    this->cost=other.cost;
  }

  virtual NDTCell* clone() const;
  virtual NDTCell* copy() const;

  inline void setCenter(const pcl::PointXYZ &cn)
  {
    center_ = cn;
  }
  inline void setDimensions(const double &xs, const double &ys, const double &zs)
  {
    xsize_ = xs;
    ysize_ = ys;
    zsize_ = zs;
  }

  inline pcl::PointXYZ getCenter() const
  {
    return center_;
  }

  inline void getCenter(double &cx,double &cy,double &cz) const
  {
    cx=center_.x;
    cy=center_.y;
    cz=center_.z;

  }




  inline void getDimensions(double &xs, double &ys, double &zs) const
  {
    xs = xsize_;
    ys = ysize_;
    zs = zsize_;
  }
  inline bool isInside(const pcl::PointXYZ pt) const
  {
    if(pt.x < center_.x-xsize_/2 || pt.x > center_.x+xsize_/2)
    {
      return false;
    }
    if(pt.y < center_.y-ysize_/2 || pt.y > center_.y+ysize_/2)
    {
      return false;
    }
    if(pt.z < center_.z-zsize_/2 || pt.z > center_.z+zsize_/2)
    {
      return false;
    }
    return true;
  }
  virtual double getDiagonal() const
  {
    return std::sqrt(xsize_*xsize_+ysize_*ysize_+zsize_*zsize_);
  }
  /**
    * Updates the current Sample mean and covariance based on
    * give new sample mean @m2 and covariance @cov2,
    * which have been computed from @numpointsindistribution number of points
    */

  void updateSampleVariance(const Eigen::Matrix3d &cov2,const Eigen::Vector3d &m2, unsigned int numpointsindistribution,
                            bool updateOccupancyFlag=true,  float max_occu=1024, unsigned int maxnumpoints=1e9);

  /**
    * Fits and updates the sample mean and covariance for the cell after the scan has been added.
    * This function updates the occupancy, it uses the given method for the update.
    * The behavior of the update can be altered with axnumpoints and occupancy_limit
    *
    * @param mode Determines the mode of the cell update - the mode defines are in the beginning of the header
    * @param maxnumpoints This adapts the cell content in the presence of non-stationary distribution. The lower the value the faster the adaptation
    * @param occupancy_limit This sets the limit for the confidence (in log-odds) that the cell can obtain.
    * @param origin The CELL_UPDATE_MODE_SAMPLE_VARIANCE_SURFACE_ESTIMATION requires to know the sensor origin, so in case You use that, you should provide the
    * 							position of the sensor, from where the measurement was taken for the method
    * @param sensor_noise A standard deviation of the sensor noise, used only by method CELL_UPDATE_MODE_SAMPLE_VARIANCE_SURFACE_ESTIMATION
    * if (maxnumpoints<=0) then the cell adaptation strategy is not used
    */
  void computeGaussian(int mode=CELL_UPDATE_MODE_SAMPLE_VARIANCE, unsigned int maxnumpoints=1e9, float occupancy_limit=255, Eigen::Vector3d origin = Eigen::Vector3d(0,0,0), double sensor_noise=0.1);

  /**
     * just updates the parameters based on points and leaves the points to cell
     */
  void computeGaussianSimple();
  /**
    * Calculates the average color for cell if the point type is pcl::PointXYZI or pcl::PointXYZRGB
    */
  void updateColorInformation();


  void rescaleCovariance();
  /**
    * Rescales the covariance to protect against near sigularities
    * and computes the inverse - This does not change class member values
    * @return true if success, false if eigen values were negative
    */
  bool rescaleCovariance(Eigen::Matrix3d &cov, Eigen::Matrix3d &invCov);

  //void updateObservation();
  void classify();

  int writeToJFF(FILE * jffout);
  int loadFromJFF(FILE * jffin);

  inline CellClass getClass() const
  {
    return cl_;
  }
  inline Eigen::Matrix3d getCov() const
  {
    return cov_;
  }
  inline Eigen::Matrix3d getInverseCov() const
  {
    return icov_;
  }
  inline Eigen::Vector3d getMean() const
  {
    return mean_;
  }
  inline Eigen::Matrix3d getEvecs() const
  {
    return evecs_;
  }
  inline Eigen::Vector3d getEvals() const
  {
    return evals_;
  }


  void setCov(const Eigen::Matrix3d &cov);

  inline void setMean(const Eigen::Vector3d &mean)
  {
    mean_ = mean;
  }
  inline void setEvals(const Eigen::Vector3d &ev)
  {
    evals_ = ev;
  }

  ///FOR TESTING ONLY
  inline void setEvecs(const Eigen::Matrix3d& ev){
    evecs_ = ev;
  }


  ///use this to set the parameters for the NDTCell. \note be careful, remember that the parameters are static, thus global
  static void setParameters(double _EVAL_ROUGH_THR   =0.1,
                            double _EVEC_INCLINED_THR=8*M_PI/18,
                            double _EVAL_FACTOR      =1000,
                            int _MIN_NB_POINTS_FOR_GAUSSIAN = 3,
                            bool _CLEAR_MIN_NB_POINTS = true,
                            bool _MIN_NB_POINTS_SET_UNIFORM = false);

    
    /**
    * Get likelihood for a given point
    */
  double getLikelihood(const pcl::PointXYZ &pt) const;

  /**
    * Adds a new point to distribution (does not update the distribution)
    * Call computeGaussian() to update the content
    */
  virtual void addPoint(const pcl::PointXYZ &pt)
  {
    points_.push_back(pt);
  }

  virtual void addPoints(pcl::PointCloud<pcl::PointXYZ> &pt)
  {
    points_.insert(points_.begin(),pt.points.begin(),pt.points.end());
  }

  /// Set the RGB value for this cell
  void setRGB(float r, float g, float b)
  {
    R=r;
    G = g;
    B = b;
  }
  void getRGB(float &r, float &g, float &b)
  {
    r = R;
    g = G;
    b = B;
  }

  /**
    * Updates the occupancy value of the cell by summing @occ_val to
    * class variable
    */
  void updateOccupancy(float occ_val, float max_occu=255.0)
  {
    occ+=occ_val;
    if(occ>max_occu) occ = max_occu;
    if(occ<-max_occu) occ = -max_occu;
    max_occu_= max_occu;
  }
  const double getMaxOccupancy()const {
    return max_occu_;
  }

  /**
    * Returns the current accumulated occupancy value
    */
  float getOccupancy()
  {
    return occ;
  }
  /**
    * Returns the current accumulated occupancy value (rescaled)
    */
  float getOccupancyRescaled()
  {
    float occupancy = 1 - 1/(1+exp(occ));
    return occupancy > 1 ? 1 : (occupancy < 0 ? 0 : occupancy);
  }

  void updateEmpty(double elik, double dist)
  {
    emptyval++;
    emptylik+=elik;
    emptydist+=dist;
  }
   const double getEmptyDist()const{return emptydist;}
   const double getEmptyLik()const{return emptylik;}
  float getDynamicLikelihood(unsigned int N)
  {
    return edata.computeSemiStaticLikelihood(N);
  }
  void setOccupancy(float occ_)
  {
    occ = occ_;
  }
  void setEmptyval(int emptyval_)
  {
    emptyval=emptyval_;
  }
  void setEventData(TEventData _ed)
  {
    edata = _ed;
  }
  const int getEmptyval() const
  {
    return emptyval;
  }
  TEventData getEventData()
  {
    return edata;
  }
  void setN(int N_)
  {
    N = N_;
  }
  int getN()
  {
    return N;
  }
  /**
    * Computes the maximum likelihood that a point moving along a line
    * defined by two points p1 and p2, gets measured agains the normaldistribution that
    * is within this cell.
    * This is used in raytracing to check if a measurement ray passes through a previously
    * mapped object (thus provides evidence of inconsistency)
    *
    * @param p1 One point along the ray
    * @param p2 second point along the ray (it must hold that p1 != p2);
    * @param &out Gives out the exact maximum likelihood point
    */
  double computeMaximumLikelihoodAlongLine(const pcl::PointXYZ &p1, const pcl::PointXYZ &p2, Eigen::Vector3d &out);

  double computeMaximumLikelihoodAlongLine(const Eigen::Vector3d &ep1, const Eigen::Vector3d &ep2, Eigen::Vector3d &out) {
    pcl::PointXYZ p1(ep1(0), ep1(1), ep1(2));
    pcl::PointXYZ p2(ep2(0), ep2(1), ep2(2));
    return computeMaximumLikelihoodAlongLine(p1, p2, out);
  }


  std::string ToString();
private:
  pcl::PointXYZ center_;
  double xsize_, ysize_, zsize_;
  Eigen::Matrix3d cov_;		/// Contains the covariance of the normal distribution
  Eigen::Matrix3d icov_;  /// Precomputed inverse covariance (updated every time the cell is updated)
  Eigen::Matrix3d evecs_; /// Eigen vectors
  Eigen::Vector3d mean_;  /// Mean of the normal distribution
  Eigen::Vector3d evals_; /// Eigen values
  CellClass cl_;
  static bool parametersSet_;													// ???
  static double EVAL_ROUGH_THR;		// = 0.1;								// ???
  static double EVEC_INCLINED_THR; 	// = cos(8*M_PI/18);//10 degree slope;	// ???
  static double EVAL_FACTOR;													// ???
  static int    MIN_NB_POINTS_FOR_GAUSSIAN;
  static bool   CLEAR_MIN_NB_POINTS;
  static bool   MIN_NB_POINTS_SET_UNIFORM;
  double d1_,d2_;
  unsigned int N; 	///Number of points used for Normal distribution estimation so far
  int emptyval;			///The number of times a cell was observed empty (using ray casting)
  double emptylik;
  double emptydist;
  float R,G,B; 			///RGB values [0..1] - Special implementations for PointXYZRGB & PointXYZI
  float occ;   			///Occupancy value stored as "Log odds" (if you wish)
  float max_occu_;
  TEventData edata;

  /**
    * Cell estimation using student-t
    */
  void studentT();

  double squareSum(const Eigen::Matrix3d &C,const Eigen::Vector3d &x){
    double sum;
    sum = C(0,0)*x(0)*x(0) + C(1,1)*x(1)*x(1) + C(2,2)*x(2)*x(2);
    sum += 2.0*C(0,1)*x(0)*x(1)   + 2.0*C(0,2)*x(0)*x(2) + 2.0*C(1,2)*x(1)*x(2);
    return sum;
  }


  void writeJFFMatrix(FILE * jffout, Eigen::Matrix3d &mat);
  void writeJFFVector(FILE * jffout, Eigen::Vector3d &vec);
  void writeJFFEventData(FILE * jffout, TEventData &evdata);
  int loadJFFMatrix(FILE * jffin, Eigen::Matrix3d &mat);
  int loadJFFVector(FILE * jffin, Eigen::Vector3d &vec);
  int loadJFFEventData(FILE * jffin, TEventData &evdata);

  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)//In order to clal this you need to register it to boost using "ar.template register_type<LazyGrid>();"
  {
      ar & cloud_index;
    ar & hasGaussian_;
    ar & cost; 	/// ndt_wavefront planner
    ar & isEmpty;	///<based on the most recent observation, is the cell seen empty (1), occupied (-1) or not at all (0)
    ar & consistency_score;
    ar & center_.x & center_.y & center_.z;
    ar & xsize_ & ysize_ & zsize_;
    ar & boost::serialization::make_array(cov_.data(), 9);
    ar & boost::serialization::make_array(icov_.data(), 9);
    ar & boost::serialization::make_array(evecs_.data(), 9);
    ar & boost::serialization::make_array(mean_.data(), 3 );
    ar & boost::serialization::make_array(evals_.data(), 3 );
    ar & parametersSet_;
    ar & EVAL_ROUGH_THR;
    ar & EVEC_INCLINED_THR;
    ar & EVAL_FACTOR;
    ar & d1_ & d2_;
    ar & N;
    ar & emptyval;
    ar & emptylik;
    ar & emptydist;
    ar & R & G & B;
    ar & occ & max_occu_;
    //ar & edata;
    ar & cl_;
  }
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW

};
};


#endif
