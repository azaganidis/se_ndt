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
	getline(infile, line);
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
	po::options_description desc("Allowed options");
    desc.add_options()
	("help", "produce help message")
	 ("pointclouds", po::value<std::vector<string> >()->multitoken(), "Point cloud files");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
	vector<string> pointcloud_files;
	if (!vm["pointclouds"].empty() && (pointcloud_files= vm["pointclouds"].as<vector<string> >()).size() >= 2) {
		///cout<<"success pointcloud read";
	}else {cout<<"pointclouds read failure";};
	int num_files=pointcloud_files.size();
	Eigen::Affine3d T;
	T.setIdentity();

		lslgeneric::NDTMatcherD2D matcher(false,false,{1,2,1,0.5});
		matcher.ITR_MAX =5;
		matcher.step_control=true;
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloudp=getCloudXYZ(pointcloud_files[0]);
	for(int i=1;i<num_files;i++)
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloudc=getCloudXYZ(pointcloud_files[i]);
		Eigen::Affine3d T_pred;
		T_pred.setIdentity();
		matcher.match(*cloudp,*cloudc,T_pred,false);
		T=T_pred*T;
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				cout<<T_pred(i,j)<<", ";
		cout<<endl;
		cloudp=cloudc;
	}

	return 0;
}

