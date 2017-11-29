#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <se_ndt/features_CI.hpp>

using namespace std;
namespace po = boost::program_options;
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
	 ("rsd-radius",po::value<float>(&rsd_radius)->default_value(0.1), "RSD: Search radius for neighbors.")
	 ("plane-radius",po::value<float>(&plane_radius)->default_value(4), "RSD: Higher radius than that is considered a plane.")
	 ("smoothness", "Calculate Smoothness.")
	 ("n-neighbors",po::value<int>(&n_neighbours)->default_value(10), "Smoothness: Number of neighbors to consider.")
	 ("r-neighbors",po::value<float>(&r_neighbours)->default_value(0.2), "Smoothness: Radius where neighbors are considered.")
	 ("pointcloud,p", po::value<std::string>(&in_p), "Point cloud file.");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
    po::notify(vm);
	if(vm.count("help"))
	{
		cout<<desc;
		return 0;
	}

	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud=getCloud<pcl::PointXYZ>(in_p,' ',false);
	getRSD(cloud,norm_radius,rsd_radius,plane_radius);
/*	std::vector<double> measure_out;
	if(vm.count("rsd")) 
	else if(vm.count("smoothness"))
	{
		if(vm.count("r-neighbors"))
			measure_out=getCornerness(cloud,r_neighbours);
		else 
			measure_out=getCornerness2(cloud,n_neighbours);
	}
	for(auto a:measure_out)
		std::cout<<a<<std::endl;
*/
	return 0;
}

