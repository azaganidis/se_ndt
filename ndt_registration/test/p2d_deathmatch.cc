#include <ndt_registration/ndt_matcher_p2d.h>
//#include <ndt_registration/ndt_matcher_d2d_2d.h>
//#include <ndt_registration/ndt_matcher_d2d.h>
#include <ndt_map/ndt_map.h>
#include <ndt_map/pointcloud_utils.h>

#include "pcl/point_cloud.h"
#include "pcl/io/pcd_io.h"
#include <cstdio>
#include <Eigen/Eigen>
#include <Eigen/Geometry>

#include <iostream>
#include <sstream>

using namespace std;

#define DODBG

double getDoubleTime()
{
    struct timeval time;
    gettimeofday(&time,NULL);
    return time.tv_sec + time.tv_usec * 1e-6;
}


int
main (int argc, char** argv)
{
    if(argc<22) {
	std::cout<<"usage "<<argv[0]<<" cloud_fixed.pcd cloud_offset.pcd T00 T01 T02 T03 T10 T11 T12 T13 T20 T21 T22 T23 T30 T31 T32 T33 max_resolution n_iterations subsample_size\n";
	return -1;
    }

    pcl::PointCloud<pcl::PointXYZ> cloud_fixed, cloud_offset, cloud_trans;
    if (pcl::io::loadPCDFile (argv[1], cloud_fixed) == -1)
    {
	cerr << "Was not able to open file \""<<argv[1]<<"\".\n";
	return 1;
    }

    if (pcl::io::loadPCDFile (argv[2], cloud_offset) == -1)
    {
	cerr << "Was not able to open file \""<<argv[2]<<"\".\n";
	return 1;
    }
    cloud_trans = cloud_offset;

    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> Tinit;
    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> Tout;
    Tout.setIdentity();

    Tinit.matrix()(0,0) = atof(argv[3]);
    Tinit.matrix()(0,1) = atof(argv[4]);
    Tinit.matrix()(0,2) = atof(argv[5]);
    Tinit.matrix()(0,3) = atof(argv[6]);
   
    Tinit.matrix()(1,0) = atof(argv[7]);
    Tinit.matrix()(1,1) = atof(argv[8]);
    Tinit.matrix()(1,2) = atof(argv[9]);
    Tinit.matrix()(1,3) = atof(argv[10]);
   
    Tinit.matrix()(2,0) = atof(argv[11]);
    Tinit.matrix()(2,1) = atof(argv[12]);
    Tinit.matrix()(2,2) = atof(argv[13]);
    Tinit.matrix()(2,3) = atof(argv[14]);

    Tinit.matrix()(3,0) = 0;
    Tinit.matrix()(3,1) = 0;
    Tinit.matrix()(3,2) = 0;
    Tinit.matrix()(3,3) = 1;
#ifdef DODBG
    std::cerr<<"Tinit:\n"<<Tinit.matrix()<<std::endl;
#endif

    double tnow, tend;
    lslgeneric::transformPointCloudInPlace(Tinit,cloud_offset);

    tnow = getDoubleTime();
    double res_max = atof(argv[19]);
    int itr_max = atoi(argv[20]);
    double subsample = atof(argv[21]);

    double __res[] = {0.5, 1, 2};
    std::vector<double> resolutions (__res, __res+sizeof(__res)/sizeof(double));
    resolutions.push_back(res_max);

    lslgeneric::NDTMatcherP2D matcherP2D(resolutions);
    matcherP2D.ITR_MAX = itr_max;
    matcherP2D.subsample_size = subsample;
    bool ret = matcherP2D.match(cloud_fixed,cloud_offset,Tout);
    
    tend = getDoubleTime();
    Tout = Tout*Tinit;
    Eigen::Matrix<double, 4, 4> m = Tout.matrix();
    std::cout<<tend-tnow<<", "<<m(0,0)<<", "<<m(0,1)<<", "<<m(0,2)<<", "<<m(0,3)<<", "
	<<m(1,0)<<", "<<m(1,1)<<", "<<m(1,2)<<", "<<m(1,3)<<", "
	<<m(2,0)<<", "<<m(2,1)<<", "<<m(2,2)<<", "<<m(2,3)<<", "
	<<m(3,0)<<", "<<m(3,1)<<", "<<m(3,2)<<", "<<m(3,3)<<std::endl;

#ifdef DODBG
    std::cerr<<"Transform: \n"<<Tout.matrix()<<std::endl;
    std::cerr<<"Time: "<<tend-tnow<<std::endl;
    //rgb point clouds
    //pcl::PointCloud<pcl::PointXYZ> cloud_trans = cloud_offset;
    lslgeneric::transformPointCloudInPlace(Tout,cloud_trans);
    pcl::PointCloud<pcl::PointXYZRGB> cloud_comb;
    pcl::PointXYZRGB red(255,0,0);
    for(int i=0; i<cloud_fixed.points.size(); ++i ) {
	red.x = cloud_fixed.points[i].x;
	red.y = cloud_fixed.points[i].y;
	red.z = cloud_fixed.points[i].z;
	cloud_comb.points.push_back(red);
    }
    pcl::PointXYZRGB green(0,200,0);
    for(int i=0; i<cloud_offset.points.size(); ++i ) {
	green.x = cloud_offset.points[i].x;
	green.y = cloud_offset.points[i].y;
	green.z = cloud_offset.points[i].z;
	cloud_comb.points.push_back(green);
    }
    pcl::PointXYZRGB blue(10,20,200);
    for(int i=0; i<cloud_trans.points.size(); ++i ) {
	blue.x = cloud_trans.points[i].x;
	blue.y = cloud_trans.points[i].y;
	blue.z = cloud_trans.points[i].z;
	cloud_comb.points.push_back(blue);
    }
    cloud_comb.width=1;
    cloud_comb.height=cloud_comb.points.size();
    cloud_comb.is_dense = false;
    pcl::io::savePCDFileBinary ("test_pcd.pcd", cloud_comb);
    pcl::io::savePCDFile ("test_pcd_fixed.pcd", cloud_fixed);
    pcl::io::savePCDFile ("test_pcd_trans.pcd", cloud_trans);
#endif

    return 0;
}
