#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <boost/program_options.hpp>
//#include <se_ndt/features_CI.hpp>
#include <se_ndt/se_ndt.hpp>
#include <se_ndt/ndt_matcher_d2d_se.h>

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
	string point_cloud_dir=string(argv[1]);
	string smoothness_dir=string(argv[2]);
	string poles_dir=string(argv[3]);
	string fname1=string(argv[4]);
	string fname2=string(argv[5]);
	pcl::PointCloud<pcl::PointXYZI>::Ptr cloud1=getCloud2(point_cloud_dir+'/'+fname1);
	pcl::PointCloud<pcl::PointXYZI>::Ptr cloud2=getCloud2(point_cloud_dir+'/'+fname2);
	std::vector<double> smoothness1=getMeasure(smoothness_dir+'/'+fname1);
	std::vector<double> poles1=getMeasure(poles_dir+'/'+fname1);
	std::vector<double> smoothness2=getMeasure(smoothness_dir+'/'+fname2);
	std::vector<double> poles2=getMeasure(poles_dir+'/'+fname2);
	size_t num_res=4;
	float size_x=50;
	float size_y=50;
	float size_z=50;
	double removeP=0.95;
		lslgeneric::NDTMatcherD2D_SE matcher;
		lslgeneric::NDTMap ***mapReading;
		lslgeneric::NDTMap ***mapReference;
		mapReading=new lslgeneric::NDTMap ** [num_res];
		mapReference=new lslgeneric::NDTMap ** [num_res];

		auto attrs={smoothness1,poles1};
		auto tails={3,3};
		int number_tails=4;

		matcher.NumInputs=number_tails;
		matcher.ITR_MAX =50;
		matcher.step_control=true;

		mapReading[0]=initMap(tails,{0.2},{size_x,size_y,size_z});
		mapReading[1]=initMap(tails,{0.5},{size_x,size_y,size_z});
		mapReading[2]=initMap(tails,{1.0},{size_x,size_y,size_z});
		mapReading[3]=initMap(tails,{2.0},{size_x,size_y,size_z});
		mapReference[0]=initMap(tails,{0.2},{size_x,size_y,size_z});
		mapReference[1]=initMap(tails,{0.5},{size_x,size_y,size_z});
		mapReference[2]=initMap(tails,{1.0},{size_x,size_y,size_z});
		mapReference[3]=initMap(tails,{2.0},{size_x,size_y,size_z});


		std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud1=getSegments(cloud1,attrs,tails,{-1,0},removeP);
		updateMap(mapReference[0],laserCloud1, number_tails);
		updateMap(mapReference[1],laserCloud1, number_tails);
		updateMap(mapReference[2],laserCloud1, number_tails);
		updateMap(mapReference[3],laserCloud1, number_tails);

		auto attrs2={smoothness2,poles2};
		std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud2=getSegments(cloud2,attrs2,tails,{-1,0},removeP);
		updateMap(mapReading[0],laserCloud2, number_tails);
		updateMap(mapReading[1],laserCloud2, number_tails);
		updateMap(mapReading[2],laserCloud2, number_tails);
		updateMap(mapReading[3],laserCloud2, number_tails);

		ET T;
		T.setIdentity();
		/*
		matcher.current_resolution=1;
		matcher.match(mapReference[1],mapReading[1],T,true);
		matcher.current_resolution=2;
		matcher.match(mapReference[2],mapReading[2],T,true);
		matcher.current_resolution=1;
		matcher.match(mapReference[1],mapReading[1],T,true);
		matcher.current_resolution=0.5;
		matcher.match(mapReference[0],mapReading[0],T,true);
		*/
		matcher.current_resolution=0.5;
		matcher.match(mapReference[1],mapReading[1],T,true);
		matcher.current_resolution=0.2;
		matcher.match(mapReference[0],mapReading[0],T,true);
		if(T.translation().transpose()*T.translation()>0.05)T.setIdentity();
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				cout<<T(i,j)<<", ";
		cout<<endl;

	return 0;
}

