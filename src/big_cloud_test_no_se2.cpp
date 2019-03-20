#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <pcl/filters/crop_box.h>
#include <se_ndt/se_ndt.hpp>
#include <pcl/common/io.h>

using namespace std;
namespace po = boost::program_options;
pcl::PointCloud<pcl::PointXYZI>::Ptr  getCloud2(string filename, bool skip=false)
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
int main(int argc, char** argv)
{
	po::options_description desc("Allowed options");
	desc.add_options()
	("help", "produce help message")
	 ("pointclouds", po::value<std::vector<string> >()->multitoken(), "Point cloud files");


    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
    po::notify(vm);
	vector<string> pointcloud_files;

	if(vm.count("help"))
	{
		cout<<desc;
		return 0;
	}
	if (!vm["pointclouds"].empty() && (pointcloud_files= vm["pointclouds"].as<vector<string> >()).size() >= 2) {
		///cout<<"success pointcloud read";
	}else {cout<<"pointclouds read failure";};
	int num_files=pointcloud_files.size();


	Eigen::Affine3d Tt;
	Tt.setIdentity();
	NDTMatch_SE matcher ({1,2,0.5},{0,1,0,2},{200,200,200},{'e'},{1},0.01,50);// :-D
	pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_refI=getCloud2(pointcloud_files[0]);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_ref(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::copyPointCloud(*cloud_refI,*cloud_ref);
	std::vector<double> rsd_ref;
	for(int j=0;j<cloud_ref->size();j++)
		rsd_ref.push_back(1);
	for(int i=0;i<num_files;i++)
	{
		pcl::PointCloud<pcl::PointXYZI>::Ptr cloud3=getCloud2(pointcloud_files[i]);
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1(new pcl::PointCloud<pcl::PointXYZ>);
		pcl::copyPointCloud(*cloud3,*cloud1);

		std::vector<double> rsd_min;
		for(int j=0;j<cloud1->size();j++)
			rsd_min.push_back(1);
		Tt=matcher.match(cloud1,cloud_ref,{rsd_min},{rsd_ref});

		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				cout<<Tt(i,j)<<", ";
		cout<<endl;
	}

	return 0;
}

