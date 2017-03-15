//#include <ndt_registration/ndt_matcher_p2d.h>
//#include <ndt_registration/ndt_matcher_d2d_2d.h>
#include <ndt_registration/ndt_matcher_d2d.h>
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
    pcl::PointCloud<pcl::PointXYZ> cloud, cloud_offset, cloud_trans;
    char fname[50];
    FILE *fout;
    double __res[] = {0.5, 1, 2, 4};
    std::vector<double> resolutions (__res, __res+sizeof(__res)/sizeof(double));

    struct timeval tv_start,tv_end,tv_reg_start,tv_reg_end;

    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> Tout;
    Tout.setIdentity();
    //printf("Working");	
    if(argc == 9)
    {

        gettimeofday(&tv_start,NULL);
        //we do a single scan to scan registration
        //TODO fix these to load pcd files
	if (pcl::io::loadPCDFile<pcl::PointXYZ> (argv[7], cloud) == -1) //* load the file
	{
	    std::cerr<<"Couldn't read file\n";
	    return (-1);
	}
	if (pcl::io::loadPCDFile<pcl::PointXYZ> (argv[8], cloud_offset) == -1) //* load the file
	{
	    std::cerr<<"Couldn't read file\n";
	    return (-1);
	}
        
        Tout =  Eigen::Translation<double,3>(xoffset,yoffset,zoffset)*
            Eigen::AngleAxis<double>(roll,Eigen::Vector3d::UnitX()) *
            Eigen::AngleAxis<double>(pitch,Eigen::Vector3d::UnitY()) *
            Eigen::AngleAxis<double>(yaw,Eigen::Vector3d::UnitZ()) ;
        
	//lslgeneric::NDTMatcherD2D_2D<pcl::PointXYZ,pcl::PointXYZ> matcherD2D(false, false, resolutions);
	lslgeneric::NDTMatcherD2D matcherD2D(false, false, resolutions);
	cloud_trans = cloud_offset;
        bool ret = matcherD2D.match(cloud,cloud_offset,Tout,true);

	std::cout<<"Transform: \n"<<Tout.matrix()<<std::endl;

	lslgeneric::transformPointCloudInPlace(Tout,cloud_trans);
	pcl::PointCloud<pcl::PointXYZRGB> cloud_comb;
	pcl::PointXYZRGB red(255,0,0);
	for(int i=0; i<cloud.points.size(); ++i ) {
	    red.x = cloud.points[i].x;
	    red.y = cloud.points[i].y;
	    red.z = cloud.points[i].z;
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

    }
}
