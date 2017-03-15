#include <NDTMatcherF2F.hh>
#include <NDTMatcher.hh>
#include <NDTMap.hh>
#include <OctTree.hh>
#include <AdaptiveOctTree.hh>
#include <PointCloudUtils.hh>

#include "ros/ros.h"
#include "pcl/point_cloud.h"
#include "sensor_msgs/PointCloud2.h"
#include "pcl/io/pcd_io.h"
#include "pcl/features/feature.h"
#include <cstdio>
#include <Eigen/Eigen>
#include <fstream>
#include "message_filters/subscriber.h"
#include "tf/message_filter.h"
#include <tf/transform_broadcaster.h>
#include <tf_conversions/tf_eigen.h>
#include <boost/circular_buffer.hpp>

class NDTMatcherNode
{
protected:
    // Our NodeHandle
    ros::NodeHandle nh_;

    // Components for tf::MessageFilter
    ros::Subscriber points2_sub_;


    // Components for publishing
    tf::TransformBroadcaster tf_;
    ros::Publisher output_pub_;

    // Use the vector as a cyclic buffer (increment with std::rotate).
    std::vector<pcl::PointCloud<pcl::PointXYZ>,Eigen::aligned_allocator<pcl::PointCloud<pcl::PointXYZ> > > pcl_buffer_;
    unsigned int nb_added_clouds_;
    boost::mutex m;
    lslgeneric::NDTMatcherF2F *matcher;

    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> incremental_pose_;
    void TransformEigenToTF(const Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> &k, tf::Transform &t)
    {
        t.setOrigin(tf::Vector3(k.matrix()(0,3), k.matrix()(1,3), k.matrix()(2,3)));
        t.setBasis(btMatrix3x3(k.matrix()(0,0), k.matrix()(0,1),k.matrix()(0,2),k.matrix()(1,0), k.matrix()(1,1),k.matrix()(1,2),k.matrix()(2,0), k.matrix()(2,1),k.matrix()(2,2)));
    };


public:
    // Constructor
    NDTMatcherNode() : nb_added_clouds_(0)

    {
        double __res[] = {0.2, 0.4, 1, 2};
        std::vector<double> resolutions (__res, __res+sizeof(__res)/sizeof(double));
        matcher = new lslgeneric::NDTMatcherF2F(false, false, false, resolutions);
        points2_sub_ = nh_.subscribe("kinect_head/camera/rgb/points", 10, &NDTMatcherNode::points2Callback, this);
        pcl_buffer_.resize(2);
        incremental_pose_.setIdentity();
    }

    ~NDTMatcherNode()
    {
        delete matcher;
    }


    // Callback
    void points2Callback(const sensor_msgs::PointCloud2::ConstPtr& msg_in)
    {
        // Add to a queue

        ROS_INFO("Got points");
        m.lock ();
        pcl::fromROSMsg (*msg_in, pcl_buffer_[0]);
        m.unlock ();

        if (nb_added_clouds_ < pcl_buffer_.size())
        {
            nb_added_clouds_++;
            std::rotate(pcl_buffer_.begin(), pcl_buffer_.begin()+1, pcl_buffer_.end());
            return;
        }

        // The most recent cloud is in [0], the secound in [1], etc.
        Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> T;
        bool ret = matcher->match(pcl_buffer_[1], pcl_buffer_[0],T);
        if (!ret)
            ROS_INFO("Registration failed!");

        incremental_pose_ = T * incremental_pose_; // * T;

        // Publish the transformation...
        tf::Transform transform;
        TransformEigenToTF(T, transform);
        //tf::TransformEigenToTF(incremental_pose_, transform);
        tf_.sendTransform(tf::StampedTransform(transform, ros::Time::now(), "world", "test"));

//      std::rotate(pcl_buffer_.begin(), pcl_buffer_.begin()+1, pcl_buffer_.end());
    }
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

int main(int argc, char **argv)
{
    ros::init(argc, argv, "ndt_matcher_node");

    NDTMatcherNode t;
    ros::spin();

    return 0;
}
