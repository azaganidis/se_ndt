#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <ndt_registration/ndt_matcher_d2d.h>

typedef Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> ET;
using namespace std;
namespace po = boost::program_options;
pcl::PointCloud<pcl::PointXYZ>::Ptr  getCloudXYZ(string filename)
{
	ifstream infile(filename); // for example
	string line = "";
	pcl::PointCloud<pcl::PointXYZ>::Ptr laserCloud(new pcl::PointCloud<pcl::PointXYZ>);
	while (getline(infile, line)){
		stringstream strstr(line);
		string word = "";
		pcl::PointXYZ point;
		getline(strstr,word, ',');
		point.x=stof(word);
		getline(strstr,word, ',');
		point.y=stof(word);
		getline(strstr,word, ',');
		point.z=stof(word);
		getline(strstr,word);
		(*laserCloud).points.push_back(point);
	}
    return laserCloud;
}
int main(int argc, char** argv)
{
	string point_cloud_dir=string(argv[1]);
	string fname1=string(argv[2]);
	string fname2=string(argv[3]);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1=getCloudXYZ(point_cloud_dir+'/'+fname1);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2=getCloudXYZ(point_cloud_dir+'/'+fname2);
	size_t num_res=2;
	float size_x=50;
	float size_y=50;
	float size_z=10;
		lslgeneric::NDTMatcherD2D matcher(false,false,{0.5,0.25});
		matcher.ITR_MAX =5;
		matcher.step_control=true;
		ET T;
		T.setIdentity();
		matcher.match(*cloud1,*cloud2,T,false);
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				cout<<T(i,j)<<", ";
		cout<<endl;

	return 0;
}

