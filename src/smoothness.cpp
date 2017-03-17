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

int main(int argc, char** argv)
{
	pcl::PointCloud<pcl::PointXYZI>::Ptr cloud=getCloud2(string(argv[1]));
	std::vector<double> smoothness;
	if(argc==3&&string(argv[2])=="2") smoothness=getCornerness2(cloud,10);
	else smoothness=getCornerness(cloud,10);
	for(auto a:smoothness)
		std::cout<<a<<std::endl;

	return 0;
}

