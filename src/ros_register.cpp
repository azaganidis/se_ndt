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
#include <std_msgs/Float32MultiArray.h>
#include <tf/transform_listener.h>
#include <tf_conversions/tf_eigen.h>
#include <message_filters/subscriber.h>
#include <message_filters/time_sequencer.h>

#include <pcl/registration/icp.h>
#include "rviz_ndt.h"


using namespace std;
namespace po = boost::program_options;
void send_transform(Eigen::Affine3d T)
{
    static tf2_ros::TransformBroadcaster br;
    geometry_msgs::TransformStamped tS=tf2::eigenToTransform(T);
    tS.header.stamp=ros::Time::now();
    tS.header.frame_id="world";
    tS.child_frame_id="velodyne";
    br.sendTransform(tS);
}
void print_transform_KITTI(Eigen::Affine3d T)
{
        fprintf(stdout, "%e %e %e %e %e %e %e %e %e %e %e %e\n", 
                T(0,0),T(0,2),T(0,1),-T(1,3),
                T(2,0),T(2,2),T(2,1),T(2,3),
                T(1,0),T(1,2),T(1,1),T(0,3));
        fflush(stdout);
}
void print_transform(Eigen::Affine3d T)
{
        fprintf(stdout, "%e %e %e %e %e %e %e %e %e %e %e %e\n", 
                T(0,0),T(0,1),T(0,2),T(0,3),
                T(1,0),T(1,1),T(1,2),T(1,3),
                T(2,0),T(2,1),T(2,2),T(2,3));
        fflush(stdout);
}
class Registration{
    public:
        static NDTMatch_SE matcher;
        Eigen::Affine3d Td,T;
        ros::Publisher pub;
        Eigen::Affine3d calib;
        pcl::IterativeClosestPoint<pcl::PointXYZI,pcl::PointXYZI> icp;
        pcl::PointCloud<pcl::PointXYZI>::Ptr prev_cloud;
//        int numreg=0;
        ros::Time r_time;
        ndt_rviz rvNDT;
        bool use_gfeat=false;
//#define STOP_AFTERN
#ifdef STOP_AFTERN
        int stop_down_timer=2;
#endif
#ifdef GROUND_TRUTH_POSE 
        tf::TransformListener tf_listener;
#endif
        //Registration(ros::Publisher& pub_):matcher({4,1,0.5},{0,1,2},{50,50,10},{'=','=','=','=','=','=','=','=','=','=','=','=','=','=','='},{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14},0.01,5)
//        Registration(ros::NodeHandle& nh, std::string topic_out, bool use_gfeat=false):
        Registration(ros::NodeHandle& nh, bool use_gfeat=false):
            //matcher({4,0.8},{0,1},{200,200,50},{'=','=','=','=','=','=','=','=','=','='},{0,1,2,3,4,5,6,7,8,9},0.01,5),
            rvNDT(nh, 2),
            use_gfeat(use_gfeat)
        {
//            pub = nh.advertise<pcl::PointCloud<pcl::PointXYZI> >(topic_out, 1000);
            matcher.setNeighbours((int )2);
            Td.setIdentity();
            T.setIdentity();
 //           calib.matrix() << 4.276802385584e-04,-9.999672484946e-01,-8.084491683471e-03,-1.198459927713e-02,-7.210626507497e-03,8.081198471645e-03,-9.999413164504e-01,-5.403984729748e-02,9.999738645903e-01,4.859485810390e-04,-7.206933692422e-03,-2.921968648686e-01,0,0,0,1;
#ifdef GL_VISUALIZE
            matcher.visualize();
#endif
            //matcher.initV(&nh);
        }
        static void pri(int i)
        {
            std::cerr<<"Writing descriptor stats."<<std::endl;
            matcher.print_vals();
            std::cerr<<"Done."<<std::endl;
        }
        void callback(const pcl::PointCloud<pcl::PointXYZI>::Ptr& cloud_msg)
        {
            pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_mv(new pcl::PointCloud<pcl::PointXYZI>);
            pcl::copyPointCloud(*cloud_msg, *cloud_mv);
//            perception_oru::transformPointCloudInPlace(calib, *cloud_mv);
            /*
            if(use_gfeat)
            {
                Eigen::Matrix<double, 1024, 1> gfeat_t;
                for(int i=0;i<256;i++)
                {
                    pcl::PointXYZI p=cloud_msg->points.at(i);
                    gfeat_t(i*4)=p.x;
                    gfeat_t(i*4+1)=p.y;
                    gfeat_t(i*4+2)=p.z;
                    gfeat_t(i*4+3)=p.intensity;
                }
                matcher.gfeat.push_back(gfeat_t);
                cloud_mv->erase(cloud_mv->begin(),cloud_mv->begin()+256);   
            }
            */
            r_time = ros::Time::now();
            //std::thread loop_check(&Registration::loop_close_check, this);
            //Td=matcher.matchFaster(Td,cloud_mv);
            //T=T*Td;
#ifdef GROUND_TRUTH_POSE 
            ros::Time t;
            tf::StampedTransform stf;
            pcl_conversions::fromPCL(cloud_msg->header.stamp, t);
            tf_listener.lookupTransform("world","velodyne", t,stf);
            tf::transformTFToEigen(stf, matcher.T);
#endif
            T=matcher.mapUpdate(cloud_mv,false);
            rvNDT.plotNDTs(matcher.toRVIZ);
//            numreg++;
            //rvNDT.plotNDTs(matcher.map, matcher.resolutions.size(), matcher.NumInputs, r_time); 
            //rvNDT.plotNDTs(matcher.map, 2, matcher.NumInputs, r_time); 
#ifdef STOP_AFTERN
            if(stop_down_timer--==0)
                ros::shutdown();
#endif
//            if(numreg%4==0){
//                pcl_conversions::toPCL(r_time,cloud_mv->header.stamp);
//                cloud_mv->header.frame_id="velodyne";
//                pub.publish(cloud_mv);
//            }
#ifndef GROUND_TRUTH_POSE 
            send_transform(T);
#endif
            ///print_transform(ts);
            print_transform_KITTI(T);
        }
};
NDTMatch_SE Registration::matcher = NDTMatch_SE(8);
int main(int argc, char** argv)
{
	string transforms;
	string p_dir;
    string topic_in="semantic";
//    string topic_out="registered";

	po::options_description desc("Allowed options");
	desc.add_options()
	("help,h", "produce help message")
	("gfeat,g", "Use gfeat")
 //   ("topic_out,t", po::value<string>(&topic_out),"Topic to publish")
	("topic_in,p", po::value<string>(&topic_in), "Point cloud files");

	po::parsed_options parsed_options = po::command_line_parser(argc, argv)
      .options(desc).style(po::command_line_style::unix_style ).run();
    po::variables_map vm;
	po::store(parsed_options, vm);
    po::notify(vm);
	if(vm.count("help"))
        { cout<<desc; return 0; }

	ros::init (argc,argv,"pub_sendt");
	ros::NodeHandle nh;
//    Registration registration(nh, topic_out,vm.count("gfeat"));
    Registration registration(nh, vm.count("gfeat"));
    signal(SIGTRAP, Registration::pri);
    message_filters::Subscriber<pcl::PointCloud<pcl::PointXYZI> > sub(nh,topic_in,200);
    message_filters::TimeSequencer<pcl::PointCloud<pcl::PointXYZI> > seq(sub,ros::Duration(5),ros::Duration(0.05),1500);
    seq.registerCallback(&Registration::callback, &registration);
//    ros::Subscriber sub = nh.subscribe(topic_in, 500, &Registration::callback, &registration);
    //ros::spin();
    ros::MultiThreadedSpinner spinner(1);
    spinner.spin();
    return 0;
}
    
