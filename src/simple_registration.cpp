#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <se_ndt/se_ndt.hpp>
#include <pcl/common/io.h>
#include <pcl/common/transforms.h>

using namespace std;
namespace po = boost::program_options;
int main(int argc,char** argv)
{
    NDTMatch_SE matcher({100,20,4,2,1},{0,1,0,2,4,3,4},{100,100,100},8,50);
    matcher.setNeighbours((int )2);//A value of 2 sets neighboors to 8. From oru
    pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_st(new pcl::PointCloud<pcl::PointXYZI>);
    pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_mv(new pcl::PointCloud<pcl::PointXYZI>);
    //Class as integer intensity value, 0--7
    Eigen::Affine3d T=matcher.simple_match(cloud_st,cloud_mv);
}
    
