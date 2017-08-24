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
	string in_p;
	float norm_radius,rsd_radius,plane_radius,r_neighbours;
	int descr,n_neighbours;
	po::options_description desc("Allowed options");
	desc.add_options()
	("help", "Produce help message.")
	 ("rsd", "Calculate RSD.")
	 ("normal-neighbors",po::value<float>(&norm_radius)->default_value(6), "RSD: Normal estimation number of neighbors.")
	 ("normal-radius",po::value<float>(&norm_radius), "RSD: Normal estimation radius-not implemented.")
	 ("rsd-radius",po::value<float>(&rsd_radius)->default_value(0.1), "RSD: Search radius for neighbors.")
	 ("plane-radius",po::value<float>(&plane_radius)->default_value(4), "RSD: Higher radius than that is considered a plane.")
	 ("descr",po::value<int>(&descr)->default_value(0), "RSD: Output descriptor: \n\t0: Minimum radius.\n1: Maximum radius.\n2: Max/Min.\n3: Min/Max.")
	 ("smoothness", "Calculate Smoothness.")
	 ("n-neighbors",po::value<int>(&n_neighbours)->default_value(10), "Smoothness: Number of neighbors to consider.")
	 ("r-neighbors",po::value<float>(&r_neighbours)->default_value(0.2), "Smoothness: Radius where neighbors are considered.")
	 ("pointcloud", po::value<std::string>(&in_p), "Point cloud file.");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
    po::notify(vm);
	if(vm.count("help"))
	{
		cout<<desc;
		return 0;
	}

	pcl::PointCloud<pcl::PointXYZI>::Ptr cloud=getCloud2(in_p);

	std::vector<double> measure_out;
	if(vm.count("rsd")) 
		measure_out=getRSD(cloud,norm_radius,rsd_radius,plane_radius,descr);
	else if(vm.count("smoothness"))
	{
		if(vm.count("r-neighbors"))
			measure_out=getCornerness(cloud,r_neighbours);
		else 
			measure_out=getCornerness2(cloud,n_neighbours);
	}
	for(auto a:measure_out)
		std::cout<<a<<std::endl;

	return 0;
}

