#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <se_ndt/features_CI.hpp>

using namespace std;
namespace po = boost::program_options;
pcl::PointCloud<pcl::PointXYZI>::Ptr  getCloud2(string filename)
{
	ifstream infile(filename); // for example
	string line = "";
	pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloud(new pcl::PointCloud<pcl::PointXYZI>);
	getline(infile, line);
	while (getline(infile, line)){
		stringstream strstr(line);
		string word = "";
		pcl::PointXYZI point;
		getline(strstr,word, ',');
		point.x=stof(word);
		getline(strstr,word, ',');
		point.y=stof(word);
		getline(strstr,word, ',');
		point.z=stof(word);
		getline(strstr,word);
		point.intensity=stof(word);
		(*laserCloud).points.push_back(point);
	}
    return laserCloud;
}
//call with: smoothness name_of_pointcloud_file.csv 
int main(int argc, char** argv)
{
	pcl::PointCloud<pcl::PointXYZI>::Ptr cloud=getCloud2(string(argv[1]));
	std::vector<double> smoothness;
	if(string(argv[2])=="R") smoothness=getRSD(cloud,std::stof(std::string(argv[3])),std::stof(std::string(argv[4])),std::stof(std::string(argv[5])),std::stof(std::string(argv[6])));
	else if(string(argv[2])=="2") smoothness=getCornerness2(cloud,10);
	else smoothness=getCornerness(cloud,10);
	for(auto a:smoothness)
		std::cout<<a<<std::endl;

	return 0;
}

