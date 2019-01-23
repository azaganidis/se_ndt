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
class Registration{
    public:
        NDTMatch_SE matcher;
        Eigen::Affine3d Td,T;
        ros::Publisher pub;
        Eigen::Affine3d calib;
        pcl::IterativeClosestPoint<pcl::PointXYZI,pcl::PointXYZI> icp;
        pcl::PointCloud<pcl::PointXYZI>::Ptr prev_cloud;
        int numreg=0;
        ros::Time r_time;
        ndt_rviz rvNDT;
//#define STOP_AFTERN
#ifdef STOP_AFTERN
        int stop_down_timer=5;
#endif
        //Registration(ros::Publisher& pub_):matcher({4,1,0.5},{0,1,2},{50,50,10},{'=','=','=','=','=','=','=','=','=','=','=','=','=','=','='},{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14},0.01,5)
        Registration(ros::NodeHandle& nh, std::string topic_out):
            //matcher({4,0.8},{0,1},{200,200,50},{'=','=','=','=','=','=','=','=','=','='},{0,1,2,3,4,5,6,7,8,9},0.01,5),
            matcher({4, 0.8},{0,1},{50,50,50},{1,1,1,1,1,1,1,1},5),
            rvNDT(nh, 2)
        {
            pub = nh.advertise<pcl::PointCloud<pcl::PointXYZI> >(topic_out, 1000);
            matcher.setNeighbours((int )2);
            Td.setIdentity();
            T.setIdentity();
            calib.matrix() << 4.276802385584e-04,-9.999672484946e-01,-8.084491683471e-03,-1.198459927713e-02,-7.210626507497e-03,8.081198471645e-03,-9.999413164504e-01,-5.403984729748e-02,9.999738645903e-01,4.859485810390e-04,-7.206933692422e-03,-2.921968648686e-01,0,0,0,1;
#ifdef GL_VISUALIZE
            matcher.visualize();
#endif
            //matcher.initV(&nh);
        }
        /*
        void loop_close_check(){
            if(poses.size()<1)
                return;
            Eigen::Vector3d pose_last=poses.back();
            perception_oru::NDTHistogram hist_last=hists.back();
            int num_poses=poses.size();
            for(int i=0;i<num_poses-50;i++)
            {
                auto pD=pose_last-poses.at(i);
                double dist=pD.dot(pD);
                dist=sqrt(dist);
                if(dist<55)
                {
                    double similarity=hist_last.getSimilarity(hists.at(i));
                    if(similarity>1)
                        continue;
                    cerr<<"LC\t"<<similarity<<"\t"<<dist<<endl;
                }
            }

        }
        */
        void callback(const pcl::PointCloud<pcl::PointXYZI>::Ptr& cloud_msg)
        {
            pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_mv(new pcl::PointCloud<pcl::PointXYZI>);
            pcl::copyPointCloud(*cloud_msg, *cloud_mv);
            clock_t begin_time = clock();
            r_time = ros::Time::now();
            //std::thread loop_check(&Registration::loop_close_check, this);
            //Td=matcher.matchFaster(Td,cloud_mv);
            //T=T*Td;
            T=matcher.mapUpdate(cloud_mv, true);
            rvNDT.plotNDTs(matcher.map, matcher.resolutions.size(), matcher.NumInputs, r_time); 
            //rvNDT.plotNDTs(matcher.map, 2, matcher.NumInputs, r_time); 
            //hists.push_back(perception_oru::NDTHistogram (matcher.map[1], 1, 40, 10, 8,2, 5));
            //poses.push_back((T*Td).translation());
            //loop_check.join();
            //T=matcher.matchFaster_OM(T,cloud_mv);
            //cerr<<"1 "<<float( clock() -begin_time ) / CLOCKS_PER_SEC<<endl;begin_time=clock();
#ifdef STOP_AFTERN
            if(stop_down_timer--==0)
                ros::shutdown();
#endif
            ///ICP
            /*
            if(numreg>0){
                pcl::PointCloud<pcl::PointXYZI> cloud_registered;
                icp.setInputSource(cloud_mv);
                icp.setInputTarget(prev_cloud);
                icp.setMaxCorrespondenceDistance(0.1);
                icp.setMaximumIterations(5);
                //Td.setIdentity();
                Eigen::Matrix4f tm_ = Td.matrix().cast<float>();
                icp.align(cloud_registered, tm_);
                Td.matrix()=tm_.cast<double>();
                Td.matrix()=icp.getFinalTransformation().cast<double>();
            }
            cout<<"2 "<<float( clock() -begin_time ) / CLOCKS_PER_SEC<<endl;begin_time=clock();
            prev_cloud=cloud_mv;
            numreg++;
            */




            pcl_conversions::toPCL(r_time,cloud_mv->header.stamp);
            cloud_mv->header.frame_id="velodyne";
            Eigen::Affine3d ts=T;
            pub.publish(cloud_mv);
            send_transform(ts);
            ts.matrix()=calib.matrix()*ts.matrix();
            //print_transform(ts);
        }
};
int main(int argc, char** argv)
{
	string transforms;
	string p_dir;
    string topic_in="semantic";
    string topic_out="registered";

	po::options_description desc("Allowed options");
	desc.add_options()
	("help,h", "produce help message")
    ("topic_out,t", po::value<string>(&topic_out),"Topic to publish")
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
    Registration registration(nh, topic_out);
    ros::Subscriber sub = nh.subscribe(topic_in, 1000, &Registration::callback, &registration);
    ros::spin();
    return 0;
}
    
