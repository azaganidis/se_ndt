#include <ndt_registration/ndt_matcher_p2d.h>
#include <ndt_registration/ndt_matcher_d2d_2d.h>
#include <ndt_registration/ndt_matcher_d2d.h>
#include <ndt_map/ndt_map.h>
#include <pointcloud_vrml/pointcloud_utils.h>

#include "pcl/point_cloud.h"
#include <cstdio>
#include <Eigen/Eigen>
#include <Eigen/Geometry>

#include <iostream>
#include <sstream>

using namespace std;

int
main (int argc, char** argv)
{
    if(argc!=9) {
	std::cout<<"usage "<<argv[0]<<" x y z r p t cloud1.wrl cloud2.wrl\n";
	return -1;
    }

    istringstream roll_c(argv[4]),pitch_c(argv[5]),yaw_c(argv[6]),xoffset_c(argv[1]),yoffset_c(argv[2]),zoffset_c(argv[3]);
    double roll,pitch,yaw,xoffset,yoffset,zoffset;
    roll_c >> roll;
    pitch_c >> pitch;
    yaw_c >> yaw;
    xoffset_c >> xoffset;
    yoffset_c >> yoffset;
    zoffset_c >> zoffset;

    printf("X %f Y %f Z %f Roll %f Pitch %f Yaw %f \n",xoffset,yoffset,zoffset,roll,pitch,yaw);	

    char fname[50];
    FILE *fout;
    double __res[] = {0.5, 1, 2, 4};
    std::vector<double> resolutions (__res, __res+sizeof(__res)/sizeof(double));

    struct timeval tv_start,tv_end,tv_reg_start,tv_reg_end;

    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> Tout;
    //Tout.setIdentity();
    Tout =  Eigen::Translation<double,3>(xoffset,yoffset,zoffset)*
	Eigen::AngleAxis<double>(roll,Eigen::Vector3d::UnitX()) *
	Eigen::AngleAxis<double>(pitch,Eigen::Vector3d::UnitY()) *
	Eigen::AngleAxis<double>(yaw,Eigen::Vector3d::UnitZ()) ;


    gettimeofday(&tv_start,NULL);
    //we do a single scan to scan registration
    lslgeneric::NDTMap<pcl::PointXYZI> ndmap(new lslgeneric::LazyGrid<pcl::PointXYZI>(0.4)), local_map(new lslgeneric::LazyGrid<pcl::PointXYZI>(0.4));
    ndmap.loadFromJFF(argv[7]);
    local_map.loadFromJFF(argv[8]);

    //lslgeneric::NDTMatcherD2D_2D<pcl::PointXYZI,pcl::PointXYZI> matcherD2D(false, false, resolutions);
    lslgeneric::NDTMatcherD2D<pcl::PointXYZI,pcl::PointXYZI> matcherD2D(false, false, resolutions);
    bool ret = matcherD2D.match(ndmap,local_map,Tout,true);

    std::cout<<"Transform: \n"<<Tout.matrix()<<std::endl;

}
