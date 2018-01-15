#include <pcl/kdtree/kdtree_flann.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <chrono>
#include <ctime>

#include <opencv/cv.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/common/transforms.h>
#include <pcl/filters/crop_box.h>
#include <pcl/io/pcd_io.h>
#include <pcl/console/parse.h>
#include <se_ndt/se_ndt.hpp>
#include <se_ndt/ndt_matcher_d2d_se.h>
#include <omp.h>

#define n_threads 12
using namespace std;

int occluded(pcl::PointXYZI a, pcl::PointXYZI b, float d)
{
	float d1 = sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
	float d2 = sqrt(b.x*b.x+b.y*b.y+b.z*b.z);
	if(d1>d2)
	{
		float dx=b.x-a.x*d2/d1;
		float dy=b.y-a.y*d2/d1;
		float dz=b.z-a.z*d2/d1;
		if(sqrt(dx*dx+dy*dy+dz*dz)/d2<d)
			return 1;
	}
	return 0;
}
pcl::PointCloud<pcl::PointXYZI>::Ptr cropIt(pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloud)
{
	pcl::PointCloud<pcl::PointXYZI>::Ptr outC(new pcl::PointCloud<pcl::PointXYZI>);
	pcl::CropBox<pcl::PointXYZI> cropBoxFilter(false);
	cropBoxFilter.setInputCloud(laserCloud);
	cropBoxFilter.setMin(Eigen::Vector4f (-0.35, -0.35, 0, 1));
	cropBoxFilter.setMax(Eigen::Vector4f ( 0.2,  0.3, 0.65, 1));
	cropBoxFilter.setNegative(true);
	cropBoxFilter.filter(*outC);
	return outC;
}

std::vector<Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> > getPose(string filename)
{
	ifstream infile(filename); // for example
	string line = "";
	std::vector< Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> >poses;
	getline(infile, line);
	while (getline(infile, line)){
	//cin>>line;
	//while (cin>>line){
		stringstream strstr(line);
		string word = "";
		Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> poseI;
		getline(strstr,word, ',');
		getline(strstr,word, ',');
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
			{
				getline(strstr,word, ',');
				poseI(i,j)=stof(word);
			}
		poses.push_back(poseI);
	}
	return poses;
}
typedef Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> ET;



