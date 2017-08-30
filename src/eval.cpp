#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <Eigen/Geometry>
#include <pcl/point_cloud.h>
#include <pcl/common/transforms.h>

using namespace std;
namespace po = boost::program_options;
using namespace Eigen;
pcl::PointCloud<pcl::PointXYZI>::Ptr  getCloud2(string filename,bool skip=false)
{
	ifstream infile(filename); // for example
	string line = "";
	pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloud(new pcl::PointCloud<pcl::PointXYZI>);
	if(skip)getline(infile, line);
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
std::vector<Eigen::Affine3d>  readTransform(istream &infile)
{
	std::vector<Eigen::Affine3d> transforms;
	string line = "";
	while(getline(infile,line))
	{
		stringstream strstr(line);
		Eigen::Affine3d T;
		T.setIdentity();
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
			{
				string word="";
				getline(strstr,word,',');
				T(i,j)=stof(word);
			}
		transforms.push_back(T);
	}
    return transforms;
}
int main(int argc, char** argv)
{
	int num;
	float dev;
	float deva;
	std::vector<std::string> pointcloud_files;
	std::vector<std::string> semantics_files;
	string out_name, transforms_name;
	po::options_description desc("Allowed options");
    desc.add_options()
	("help", "produce help message")
	 ("pointclouds", po::value<std::vector<string> >()->multitoken(), "Point cloud files")
	 ("transforms", po::value<std::string>(&transforms_name), "Transforms file")
	 ("out", po::value<std::string>(&out_name), "Output file prefix")
	("skip", "skip point cloud first line");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
	if(vm.count("help"))
	{
		cout<<desc;
		return 0;
	}

	bool skip=false;
	if(vm.count("skip"))
		skip=true;

	if (!vm["pointclouds"].empty() && (pointcloud_files= vm["pointclouds"].as<vector<string> >()).size() > 0) {
	}else {cout<<"pointclouds read failure";};

	ifstream in_trans;
	if(vm.count("transforms")) { in_trans.open(transforms_name, ifstream::in); }else{cout<<"no transform file"<<endl;return -1;};


	Eigen::Matrix4d I;
	Eigen::Matrix3f I_f;
	I.setIdentity();
	I_f.setIdentity();

	std::vector<Eigen::Affine3d> transforms = readTransform(in_trans);
	ofstream transf_out(out_name+"tr.csv");
	for(int i=0;i<transforms.size();i++)
	{
		Eigen::Matrix4d TM= transforms[i].matrix();
		Eigen::Matrix4d TI= I.colPivHouseholderQr().solve(TM);
		for(int k=0;k<4;k++)
			for(int j=0;j<4;j++)
				transf_out<<TI(k,j)<<", ";
		transf_out<<endl;
		pcl::PointCloud<pcl::PointXYZI>::Ptr source_cloud=getCloud2(pointcloud_files[i],skip);
		pcl::PointCloud<pcl::PointXYZI>::Ptr transformed_cloud (new pcl::PointCloud<pcl::PointXYZI> ());
		pcl::transformPointCloud (*source_cloud, *transformed_cloud,TI);
		ofstream cloud_out(out_name+"pcl_"+to_string(i)+".csv");
		for(int j=0;j<transformed_cloud->size();j++)
		{
			cloud_out<<transformed_cloud->points[j].x<<", ";
			cloud_out<<transformed_cloud->points[j].y<<", ";
			cloud_out<<transformed_cloud->points[j].z<<", ";
			cloud_out<<transformed_cloud->points[j].intensity<<endl;
		}
		cloud_out.close();
			
	}
}

