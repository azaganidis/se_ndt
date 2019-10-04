#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <pcl/filters/crop_box.h>
#include <se_ndt/se_ndt.hpp>
#include <pcl/common/io.h>
#include <pcl/common/transforms.h>
#include <pcl/filters/voxel_grid.h>
#include <ndt_registration/ndt_matcher_d2d.h>
#include <ctime>

#include <ros/ros.h>
#include <pcl_ros/point_cloud.h>
#include <pcl/point_types.h>
#include <tf2_ros/transform_broadcaster.h>
#include <geometry_msgs/TransformStamped.h>
#include <tf2_eigen/tf2_eigen.h>
#include <pcl_conversions/pcl_conversions.h>

#include <pcl/registration/icp.h>
#include "rviz_ndt.h"
#include "se_ndt/ndt_histogram.h"


using namespace std;
namespace po = boost::program_options;
void send_transform(Eigen::Affine3d &T)
{
    static tf2_ros::TransformBroadcaster br;
    geometry_msgs::TransformStamped tS=tf2::eigenToTransform(T);
    tS.header.stamp=ros::Time::now();
    tS.header.frame_id="world";
    tS.child_frame_id="velodyne";
    br.sendTransform(tS);
}
void print_transform(Eigen::Affine3d &T)
{
        fprintf(stdout, "%e %e %e %e %e %e %e %e %e %e %e %e\n", 
                T(0,0),T(0,1),T(0,2),T(0,3),
                T(1,0),T(1,1),T(1,2),T(1,3),
                T(2,0),T(2,1),T(2,2),T(2,3));
        fflush(stdout);
}
class MapLoader{
    public:
        NDTMatch_SE matcher;
        ros::Time r_time;
        ndt_rviz rvNDT;
        MapLoader(ros::NodeHandle& nh):
            matcher({4,0.8},{0,1},{50,50,10},{1,1,1,1,1,1,1,1},5),
            rvNDT(nh, 2)
        {
            matcher.setNeighbours((int )2);
        }
        void load_pose(int start_index, int stop_index, double* poseD)
        {
            clock_t begin_time = clock();
            r_time = ros::Time::now();
            pcl::PointXYZL pose;
            pose.x=poseD[0];
            pose.y=poseD[1];
            pose.z=poseD[2];
            matcher.loadMap(start_index,stop_index, pose);
            //T=matcher.matchFaster_OM(T,cloud_mv);
            rvNDT.plotNDTs(matcher.toRVIZ);
            //cerr<<float( clock() -begin_time ) / CLOCKS_PER_SEC<<endl;begin_time=clock();
        }
};
int main(int argc, char** argv)
{
    if(argc!=6)
    {
        std::cerr<<"syntax: lt start_index stop_index posex posey posez"<<std::endl;
        return 0;
    }
	ros::init (argc,argv,"pub_sendt_r");
	ros::NodeHandle nh;
    MapLoader loader(nh);
    double pose[3];
    int start_index=atoi(argv[1]);
    int stop_index=atoi(argv[2]);
    pose[0]=atof(argv[3]);
    pose[1]=atof(argv[4]);
    pose[2]=atof(argv[5]);
    std::cerr<<"GOOD"<<std::endl;
    loader.load_pose(start_index, stop_index, pose);
    return 0;
}
    
