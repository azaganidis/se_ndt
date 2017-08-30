#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <pcl/point_cloud.h>
#include <pcl/common/io.h>
#include <pcl/keypoints/uniform_sampling.h>

using namespace std;
namespace po = boost::program_options;
pcl::PointCloud<pcl::PointXYZI>::Ptr  getCloud(string filename,string se_filename, bool skip=false)
{
	ifstream infile(filename);
	ifstream se_infile(se_filename);
	string line = "";
	string se_line = "";
	pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloud(new pcl::PointCloud<pcl::PointXYZI>);
	if(skip)getline(infile, line);
	while (getline(infile, line)&&getline(se_infile, se_line)){
		stringstream strstr(line);
		string word = "";
		pcl::PointXYZI point;
		getline(strstr,word, ',');
		point.x=stof(word);
		getline(strstr,word, ',');
		point.y=stof(word);
		getline(strstr,word, ',');
		point.z=stof(word);
		point.intensity=stof(se_line);
		(*laserCloud).points.push_back(point);
	}
    return laserCloud;
}
int main(int argc, char** argv)
{
	string cloud_name, sem_name,outc,outs;
	float leaf;
	po::options_description desc("Allowed options");
	desc.add_options()
	("help", "produce help message")
	 ("cloud", po::value<string>(&cloud_name), "Point cloud file")
	 ("out_cld", po::value<string>(&outc), "Point cloud out file")
	 ("out_sem", po::value<string>(&outs), "Semantics out file")
	 ("leaf", po::value<float>(&leaf)->default_value(0.3), "Leaf size")
	 ("sem", po::value<string>(&sem_name), "Point cloud reference file");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
    po::notify(vm);
	if(vm.count("help")) { cout<<desc; return 0; }
	pcl::PointCloud<pcl::PointXYZI>::Ptr cloud=getCloud(cloud_name,sem_name);
	pcl::PointCloud<pcl::PointXYZI> cloud_out;

	pcl::UniformSampling<pcl::PointXYZI> us;
	us.setInputCloud (cloud);
	us.setRadiusSearch (leaf);
	us.filter (cloud_out);
	ofstream out_f(outc);
	ofstream out_s(outs);
	for(int i=0;i<cloud_out.size();i++)
	{
		out_f<<cloud_out.points[i].x<<", ";
		out_f<<cloud_out.points[i].y<<", ";
		out_f<<cloud_out.points[i].z<<", ";
		out_f<<cloud_out.points[i].intensity<<endl;
		out_s<<cloud_out.points[i].intensity<<endl;
	}
	return 0;
}

