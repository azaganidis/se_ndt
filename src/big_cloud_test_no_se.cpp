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
	float p;
    desc.add_options()
	("help", "produce help message")
	 ("p", po::value<float>(&p), "P")
	 ("pointclouds", po::value<std::vector<string> >()->multitoken(), "Point cloud files");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
	vector<string> pointcloud_files;
	if (!vm["pointclouds"].empty() && (pointcloud_files= vm["pointclouds"].as<vector<string> >()).size() >= 2) {
		///cout<<"success pointcloud read";
	}else {cout<<"pointclouds read failure";};
	int num_files=pointcloud_files.size();

		lslgeneric::NDTMatcherD2D matcher(false,false,{5,4,20,10,180});
		//lslgeneric::NDTMatcherD2D matcher(false,false,{1,5,4,20,10,90});
	//NDTMatch_SE matcher ({200,60,200,70,50,5},{0,1,2,3,4,5},{200,200,200},{'u'},{0},0.01,50);// :-D
	//NDTMatch_SE matcher ({100,20,100,4,1,2},{0,1,2,3,4,5},{100,100,100},{'u'},{0},0.01,50);// :-D
		matcher.ITR_MAX =50;
		matcher.step_control=true;
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloudp=getCloudXYZ(pointcloud_files[0]);
	for(int i=0;i<num_files;i++)
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloudc=getCloudXYZ(pointcloud_files[i]);
		Eigen::Affine3d T_pred;
		T_pred.setIdentity();
		matcher.match(*cloudc,*cloudp,T_pred,false);
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				cout<<T_pred(i,j)<<", ";
		cout<<endl;
	}

	return 0;
}

