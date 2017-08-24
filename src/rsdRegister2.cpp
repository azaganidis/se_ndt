#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <pcl/filters/crop_box.h>
#include <se_ndt/ndt_fuser_hmt_se.h>
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
pcl::PointCloud<pcl::PointXYZI>::Ptr cropIt(pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloud, std::vector<float>& b)
{
	pcl::PointCloud<pcl::PointXYZI>::Ptr outC(new pcl::PointCloud<pcl::PointXYZI>);
	pcl::CropBox<pcl::PointXYZI> cropBoxFilter(false);
	cropBoxFilter.setInputCloud(laserCloud);
	cropBoxFilter.setMin(Eigen::Vector4f (b[0],b[1],b[2], 1));
	cropBoxFilter.setMax(Eigen::Vector4f (b[3],b[4],b[5], 1));
	cropBoxFilter.setNegative(true);
	cropBoxFilter.filter(*outC);
	return outC;
}

int main(int argc, char** argv)
{
	po::options_description desc("Allowed options");
	desc.add_options()
	("help", "produce help message")
	("skip", "skip point cloud first line")
	("nc", "Do not concatenate transforms")
	 ("pointclouds", po::value<std::vector<string> >()->multitoken(), "Point cloud files")
	 ("b", po::value<std::vector<float> >()->multitoken(), "Bounding box--Atention! not working yet!")
	 ("sem1", po::value<std::vector<string> >()->multitoken(), "First semantic input files")
	 ("sem2", po::value<std::vector<string> >()->multitoken(), "Second semantic input files");


    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
    po::notify(vm);
	vector<string> rsd_min_files;
	vector<string> rsd_max_files;
	vector<string> pointcloud_files;
	vector<float> box;
	bool skip=false;
	if(vm.count("skip"))
		skip=true;
	bool conc=true;
	if(vm.count("nc"))
		conc=false;
	if(vm.count("help"))
	{
		cout<<desc;
		return 0;
	}
	if (!vm["pointclouds"].empty() && (pointcloud_files= vm["pointclouds"].as<vector<string> >()).size() >= 2) {
		///cout<<"success pointcloud read";
	}else {cout<<"pointclouds read failure";};
	if (!vm["sem1"].empty() && (rsd_min_files= vm["sem1"].as<vector<string> >()).size() >= 2) {
		///cout<<"success poles read";
	}else {cout<<"sem1 read failure";};
	if (!vm["sem2"].empty() && (rsd_max_files= vm["sem2"].as<vector<string> >()).size() >= 2) {
		///cout<<"success poles read";
	}else {cout<<"sem2 read failure";};
	bool h_box=false;
	if(vm.count("b")) {h_box=true;box=vm["b"].as<vector<float> >();if(box.size()!=6){cout<<"Wrong box size! must be 6."<<endl;return -1;}}
	int num_files=min(rsd_min_files.size(),pointcloud_files.size());


	Eigen::Affine3d T;
	Eigen::Affine3d Tt;
	T.setIdentity();
	Tt.setIdentity();
	//NDTMatch_SE matcher ({0.5,0.1,0.05},{0,1,0,1,2},{25,25,10},{3},{-1},0.60,25);
	NDTMatch_SE matcher ({1,2},{1,0},{80,80,20},{3,3},{-1,-1},0.50,25);
	//lslgeneric::NDTFuserHMT_SE matcher (the_initial_pose,{the_resolutions},{the_order_with which_the_resolutions_are_used},{the_size_of_the_map},{the_tail_segments},{ignore_values},reject_percentage,number_of_iterations);
	for(int i=0;i<num_files;i++)
	{
		pcl::PointCloud<pcl::PointXYZI>::Ptr cloud3=h_box?cropIt(getCloud2(pointcloud_files[i],skip),box):getCloud2(pointcloud_files[i],skip);
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1(new pcl::PointCloud<pcl::PointXYZ>);
		pcl::copyPointCloud(*cloud3,*cloud1);
		std::vector<double> rsd_min=getMeasure(rsd_min_files[i]);
		std::vector<double> rsd_max=getMeasure(rsd_max_files[i]);
		Tt=matcher.match(cloud1,{rsd_min,rsd_max});
		T=T*Tt;

		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				cout<<(conc?T(i,j):Tt(i,j))<<", ";
		cout<<endl;
	}

	return 0;
}

