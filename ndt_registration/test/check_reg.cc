/**
 * \file Test program for checking registration quality. Reads two
 * point clouds and a transformation, and computes a measure of
 * whether they are aligned or not.
 * 
 * \author Martin Magnusson
 * \date   Created 2014-10-14
 */

#include "pcl/io/pcd_io.h"
#include "pcl/point_cloud.h"
#include <Eigen/Geometry> // for eulerAngles
#include <iostream>
#include <ndt_registration/ndt_matcher_d2d.h>
#include <ndt_registration/ndt_matcher_p2d.h>
#include <string>
//#include <boost/program_options.hpp>
//#include <cstdlib>


using namespace std;
using namespace lslgeneric;

int main( int argc, char** argv )
{
  if (argc < 8+1)
  {
    std::cerr << "Usage: " << argv[0] << " <pc1> <pc2> <3xtranslation> <3xrotation>\n";
    return -1;
  }
  int i = 1;
  std::string file_fixed  (argv[i++]);
  std::string file_offset (argv[i++]);

  double tx (atof(argv[i++]));
  double ty (atof(argv[i++]));
  double tz (atof(argv[i++]));
  double rx (atof(argv[i++]));
  double ry (atof(argv[i++]));
  double rz (atof(argv[i++]));

  Eigen::Transform<double,3,Eigen::Affine> T;
  T = Eigen::Translation3d( tx, ty, tz );
   //Eigen::AngleAxisd( 0.0, Eigen::Vector3d( 1.0, 0.0, 0.0 ) );
  //  Eigen::eulerAngles( rx, ry, rz );
    //  T.setIdentity();

  pcl::PointCloud<pcl::PointXYZ> cloud_fixed;
  if (pcl::io::loadPCDFile (file_fixed, cloud_fixed) == -1)
  {
    cerr << "Was not able to open file \"" << file_fixed << "\".\n";
    return 1;
  }
  pcl::PointCloud<pcl::PointXYZ> cloud_offset;
  if (pcl::io::loadPCDFile (file_offset, cloud_offset) == -1)
  {
    cerr << "Was not able to open file \"" << file_offset << "\".\n";
    return 1;
  }

  lslgeneric::transformPointCloudInPlace( T, cloud_offset );
  
  std::vector<double> resolutions;
  resolutions.push_back( 1 );
  NDTMatcherP2D matcherP2D(resolutions);
  NDTMatcherD2D matcherD2D(false, false, resolutions);

  Eigen::Matrix<double,6,6> covP2D;
  Eigen::MatrixXd covD2D(6,6);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,6,6> > Sol ;
  //Eigen::SelfAdjointEigenSolver<Eigen::Matrix6d > Sol ();
  Eigen::Matrix<double,6,1> evals;
  
  // Check the alignment before registration
  //matcherP2D.check( cloud_fixed, cloud_offset, T );
  matcherP2D.covariance( cloud_fixed, cloud_offset, T, covP2D );
  cout << "Cov(P2D)\n" << covP2D << "\n";
  Sol.compute( covP2D );
  evals = Sol.eigenvalues().real();
  cout << "Cov(P2D) eigenvalues "
       << evals[0] << " "
       << evals[1] << " "
       << evals[2] << " "
       << evals[3] << " "
       << evals[4] << " "
       << evals[5] << " "
       << "\n";
  matcherD2D.covariance( cloud_fixed, cloud_offset, T, covD2D );
  Sol.compute( covD2D );
  evals = Sol.eigenvalues().real();
  cout << "Cov(D2D) eigenvalues "
       << evals[0] << " "
       << evals[1] << " "
       << evals[2] << " "
       << evals[3] << " "
       << evals[4] << " "
       << evals[5] << " "
       << "\n";

  // (optionally) perform registration
  matcherP2D.ITR_MAX = 0;
  matcherP2D.subsample_size = 0.4;
  bool converged = matcherP2D.match( cloud_fixed, cloud_offset, T );
  if (not converged)
    cerr << "warning: matcher terminated before convergence\n";

  cout << "Result T: " << T.matrix() << "\n";
  
  // Now check the alignment
  //matcherP2D.check( cloud_fixed, cloud_offset, T );
  matcherP2D.covariance( cloud_fixed, cloud_offset, T, covP2D );
  cout << "Cov(P2D)\n" << covP2D << "\n";
  Sol.compute(covP2D);
  evals = Sol.eigenvalues().real();
  cout << "Cov(P2D) eigenvalues "
       << evals[0] << " "
       << evals[1] << " "
       << evals[2] << " "
       << evals[3] << " "
       << evals[4] << " "
       << evals[5] << " "
       << "\n";
  matcherD2D.covariance( cloud_fixed, cloud_offset, T, covD2D );
  Sol.compute( covD2D );
  evals = Sol.eigenvalues().real();
  cout << "Cov(D2D) eigenvalues "
       << evals[0] << " "
       << evals[1] << " "
       << evals[2] << " "
       << evals[3] << " "
       << evals[4] << " "
       << evals[5] << " "
       << "\n";

  // pcl::io::savePCDFile ("test_pcd_fixed.pcd", cloud_fixed);
  // pcl::io::savePCDFile ("test_pcd_trans.pcd", cloud_trans);
  
  return 0;
}


