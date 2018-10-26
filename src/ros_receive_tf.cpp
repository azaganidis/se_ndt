#include <iostream>
#include <string>
#include <ros/ros.h>
#include <pcl/registration/icp.h>
#include <std_msgs/Float32MultiArray.h>
#include <std_msgs/MultiArrayDimension.h>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "rviz_ndt.h"
#include <vector>


using namespace std;
class RvizShow{
    public:
        ros::Publisher pub;
        Eigen::Affine3d calib;
        pcl::IterativeClosestPoint<pcl::PointXYZI,pcl::PointXYZI> icp;
        pcl::PointCloud<pcl::PointXYZI>::Ptr prev_cloud;
        int numreg=0;
        ros::Time r_time;
        ndt_rviz rvNDT;
        RvizShow(ros::NodeHandle& nh): rvNDT(nh, 1)
        {
        }
        void callback(const std_msgs::Float32MultiArray::ConstPtr& msg)
        {
            rvNDT.dur=ros::Duration(5);
            auto r_time = ros::Time::now();
            vector<float> data=msg->data;
            int NC = data.size()/12;
            for(int i=0;i<NC;i++)
            {
                float *v = &data[i*12];
                Eigen::Map<Eigen::Matrix3f> mat(v, 3, 3);
                Eigen::Matrix3d C = mat.cast<double> ();
                float *e = &data[i*12+9];
                Eigen::Map<Eigen::Vector3f> mat2(e, 3);
                Eigen::Vector3d M = mat2.cast<double> ();
                //cout<<M.transpose()<<endl;
                rvNDT.show_cell(C,M,125,r_time,0,0,i);
            }

            //rvNDT.show_cell(matcher.map, 2, matcher.NumInputs, r_time); 
        }
};
int main(int argc, char** argv)
{
	string transforms;
	string p_dir;
    string topic_in="NDTs";

	ros::init (argc,argv,"show_rviz");
	ros::NodeHandle nh;
    RvizShow rs(nh);
    ros::Subscriber sub = nh.subscribe(topic_in, 1000, &RvizShow::callback, &rs);
    ros::spin();
    return 0;
}
    
