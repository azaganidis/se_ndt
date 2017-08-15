#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <se_ndt/ndt_fuser_hmt_se.h>
#include <pcl/common/io.h>

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
vector<double>  getMeasure(string filename)
{
	ifstream infile(filename); // for example
	string line = "";
	vector<double> measure;
	while (getline(infile, line)){
		measure.push_back(stod(line));
	}
    return measure;
}


int main(int argc, char** argv)
{
	po::options_description desc("Allowed options");
    desc.add_options()
	("help", "produce help message")
	 ("pointclouds", po::value<std::vector<string> >()->multitoken(), "Point cloud files")
	 ("rsd", po::value<std::vector<string> >()->multitoken(), "RSD files");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
	vector<string> rsd_files;
	vector<string> pointcloud_files;
	if (!vm["pointclouds"].empty() && (pointcloud_files= vm["pointclouds"].as<vector<string> >()).size() >= 2) {
		///cout<<"success pointcloud read";
	}else {cout<<"pointclouds read failure";};
	if (!vm["rsd"].empty() && (rsd_files= vm["rsd"].as<vector<string> >()).size() >= 2) {
		///cout<<"success poles read";
	}else {cout<<"rsd read failure";};
	int num_files=min(rsd_files.size(),pointcloud_files.size());


	Eigen::Affine3d T;
	Eigen::Affine3d Tt;
	T.setIdentity();
	Tt.setIdentity();
	NDTMatch_SE matcher ({0.5,0.1,0.05},{0,1,0,1,2},{25,25,10},{3},{-1},0.60,25);
	//lslgeneric::NDTFuserHMT_SE matcher (the_initial_pose,{the_resolutions},{the_order_with which_the_resolutions_are_used},{the_size_of_the_map},{the_tail_segments},{ignore_values},reject_percentage,number_of_iterations);
	for(int i=0;i<num_files;i++)
	{
		pcl::PointCloud<pcl::PointXYZI>::Ptr cloud3=getCloud2(pointcloud_files[i]);
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1(new pcl::PointCloud<pcl::PointXYZ>);
		pcl::copyPointCloud(*cloud3,*cloud1);
		std::vector<double> rsd=getMeasure(rsd_files[i]);
		Tt=matcher.match(cloud1,{rsd});
		T=T*Tt;

		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				cout<<T(i,j)<<", ";
		cout<<endl;
	}

	return 0;
}

