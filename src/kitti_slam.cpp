#include <iostream>
#include <string>
#include <map>
#include <fstream>
#include <boost/program_options.hpp>
#include <se_ndt/se_ndt.hpp>
#include <ctime>
#define PI 3.14159265
std::map<int,int> get_label_map(){
    std::map<int,int> lm;
    lm[0 ]= 0;     // "unlabeled"
    lm[1 ]= 0;     // "outlier" mapped to "unlabeled" --------------------------mapped
    lm[10]= -1;     // "car"
    lm[11]= -1;     // "bicycle"
    lm[13]= -1;     // "bus" mapped to "other-vehicle" --------------------------mapped
    lm[15]= -1;     // "motorcycle"
    lm[16]= -1;     // "on-rails" mapped to "other-vehicle" ---------------------mapped
    lm[18]= -1;     // "truck"
    lm[20]= -1;     // "other-vehicle"
    lm[30]= -1;     // "person"
    lm[31]= -1;     // "bicyclist"
    lm[32]= -1;     // "motorcyclist"
    lm[40]= 1;     // "road"
    lm[44]= 1;     // "parking"
    lm[48]= 2;    // "sidewalk"
    lm[49]= 3;    // "other-ground"
    lm[50]= 4;    // "building"
    lm[51]= 5;    // "fence"
    lm[52]= 0;     // "other-structure" mapped to "unlabeled" ------------------mapped
    lm[60]= 6;     // "lane-marking" to "road" ---------------------------------mapped
    lm[70]= 7;    // "vegetation"
    lm[71]= 8;    // "trunk"
    lm[72]= 9;    // "terrain"
    lm[80]= 10;    // "pole"
    lm[81]= 11;    // "traffic-sign"
    lm[99]= -1;     // "other-object" to "unlabeled" ----------------------------mapped
    lm[252]= -1;    // "moving-car"
    lm[253]= -1;    // "moving-bicyclist"
    lm[254]= -1;    // "moving-person"
    lm[255]= -1;    // "moving-motorcyclist"
    lm[256]= -1;    // "moving-on-rails" mapped to "moving-other-vehicle" ------mapped
    lm[257]= -1;    // "moving-bus" mapped to "moving-other-vehicle" -----------mapped
    lm[258]= -1;    // "moving-truck"
    lm[259]= -1;    // "moving-other-vehicle"
    return lm;
}

using namespace std;
namespace po = boost::program_options;

pcl::PointCloud<pcl::PointXYZI>::Ptr getB(string filenameP, string filenameL, std::map<int,int> &label_map)
{
	pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloud(new pcl::PointCloud<pcl::PointXYZI>);
    std::ifstream inFileP(filenameP, ios::in|ios::binary);
    std::ifstream inFileL(filenameL, ios::in|ios::binary);
    if(inFileP.good() && inFileL.good())
    {
        float in_float[4];
        int label;
        while(!inFileP.eof() && !inFileL.eof())
        {
            inFileP.read((char *)in_float, 4*sizeof(float)); 
            pcl::PointXYZI point;
            point.x=in_float[0];
            point.y=in_float[1];
            point.z=in_float[2];
            point.intensity=in_float[3];
            inFileL.read((char *)&label, sizeof(int));
            point.intensity= label_map[label];
            if(point.intensity==-1)
                continue;
            (*laserCloud).points.push_back(point);
        }
    }
    return laserCloud;
}


int main(int argc, char** argv)
{
	string p_dir;
    auto label_map=get_label_map();

	po::options_description desc("Allowed options");
	desc.add_options()
	("help,h", "produce help message")
	("visualize,v", "Visualize with opengl")
	("labels,l", po::value<std::vector<string> >()->multitoken(), "Labels")
	("pointclouds,p", po::value<std::vector<string> >()->multitoken(), "Point cloud files");

	po::parsed_options parsed_options = po::command_line_parser(argc, argv)
      .options(desc).style(po::command_line_style::unix_style ).run();
    po::variables_map vm;
	po::store(parsed_options, vm);
    po::notify(vm);
	if(vm.count("help")||!vm.count("pointclouds")||!vm.count("labels"))
	{
        std::cout<<desc;
		return 0;
	}
	vector<string> pointcloud_files;
	pointcloud_files= vm["pointclouds"].as<vector<string> >();
	vector<string> label_files;
	label_files= vm["labels"].as<vector<string> >();
    NDTMatch_SE matcher ({4,0.8},{0,1},{80,80},12,50);
    if(vm.count("visualize"))
        matcher.visualize();

    int s_point=pointcloud_files.size();
	for(int t=0; t<s_point; t++)
	{
		pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_mv = getB(pointcloud_files[t], label_files[t], label_map);
		matcher.slam(cloud_mv);
    }
	return 0;
}

