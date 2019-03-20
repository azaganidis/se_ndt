#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <pcl/filters/crop_box.h>
#include <se_ndt/se_ndt.hpp>
#include <pcl/common/io.h>
#include <ndt_registration/ndt_matcher_d2d.h>

using namespace std;
namespace po = boost::program_options;
std::tuple<pcl::PointCloud<pcl::PointXYZ>::Ptr,std::vector<double> >  getCloud(string filename, int n_useless)
{
	pcl::PointCloud<pcl::PointXYZ>::Ptr laserCloud(new pcl::PointCloud<pcl::PointXYZ>);
	std::vector<double> labels;
	ifstream in(filename); // for example
	if(!in)std::cerr<<"could not open pointcloud file :"<<filename<<std::endl;
	string line;
	float useless;
	while (getline(in, line)){
		pcl::PointXYZ point;
		double label;
		stringstream sin(line);
		sin>>point.x>>point.y>>point.z;
		(*laserCloud).points.push_back(point);
		for(int i=0;i<n_useless;i++)
			sin>>useless;
		sin>>label;
		labels.push_back(label);

	}
    return std::make_tuple(laserCloud, labels);
}
std::vector<Eigen::Affine3d>  readTransform(string fname)
{
	ifstream infile(fname);
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
template<int M, template<typename> class F = std::less>
struct TupleCompare
{
	template<typename T>
	bool operator()(T const &t1, T const &t2)
	{
		return F<typename tuple_element<M,T>::type>()(std::get<M>(t1), std::get<M>(t2));
	}
};
double ang_diff(Eigen::Matrix4d A, Eigen::Matrix4d B)
{
		Eigen::Matrix4d T = A.colPivHouseholderQr().solve(B);
		Eigen::Vector3d translation = T.block<3,1>(0,3);
		return  acos((T.block<3,3>(0,0).trace()-1)/2);
}
double pos_diff(Eigen::Matrix4d A, Eigen::Matrix4d B)
{
		Eigen::Matrix4d T = A.colPivHouseholderQr().solve(B);
		Eigen::Vector3d translation = T.block<3,1>(0,3);
		return  translation.dot(translation);
}
std::vector<std::tuple<int,int,double,double,Eigen::Matrix4d> > getThat(string in_trans1, bool inverse_file=false)
{
	std::vector<Eigen::Affine3d> transforms1 = readTransform(in_trans1);
	std::vector<std::tuple<int,int,double,double,Eigen::Matrix4d> > results; 
	for(int i=0;i<transforms1.size();i++)
	{
		for(int j=i+1;j<transforms1.size();j++)
		{
			Eigen::Matrix4d A,B;
			if(inverse_file)
			{
				A = transforms1[i].inverse().matrix();
				B = transforms1[j].inverse().matrix();
			}
			else
			{
				A = transforms1[i].matrix();
				B = transforms1[j].matrix();
			}
			Eigen::Matrix4d T= A.colPivHouseholderQr().solve(B);
			Eigen::Vector3d translation = T.block<3,1>(0,3);
			double dT = translation.dot(translation);
			double dR = acos((T.block<3,3>(0,0).trace()-1)/2);
			results.push_back(std::make_tuple(i,j,sqrt(dT),dR,T));

		}
	}
	sort(results.begin(),results.end(), TupleCompare<2>());
	double max_distance=3;
	int i_max=50;
	std::vector<std::tuple<int,int,double,double,Eigen::Matrix4d> > results_distributed; 
	std::vector<std::tuple<int,int,double,double,Eigen::Matrix4d> >::iterator r_iterator=results.begin(); 
	results_distributed.push_back(*r_iterator);
	float distance_step=max_distance/i_max;
	while(std::get<2>(results_distributed.back())<max_distance)
	{
		while(std::get<2>(*r_iterator)<distance_step+std::get<2>(results_distributed.back())&&r_iterator!=std::prev(results.end()))
			r_iterator++;
		results_distributed.push_back(*r_iterator);
	}
	return results_distributed;
}

int main(int argc, char** argv)
{
	string transforms;
	float p;
	string p_dir;

	po::options_description desc("Allowed options");
	desc.add_options()
	("help", "produce help message")
	("pointclouds,p", po::value<std::vector<string> >()->multitoken(), "Point cloud files");

	po::parsed_options parsed_options = po::command_line_parser(argc, argv)
      .options(desc).style(po::command_line_style::unix_style ).run();
    po::variables_map vm;
	po::store(parsed_options, vm);
    po::notify(vm);
	if(vm.count("help")||!vm.count("pointclouds"))
	{
		cout<<desc;
		return 0;
	}
	vector<string> pointcloud_files;
	pointcloud_files= vm["pointclouds"].as<vector<string> >();
	NDTMatch_SE matcher ({100,20,100,4,1,2},{0,1,2,3,4,5,4},{100,100,100},{'=','=','=','=','=','=','=','=','='},{0,1,2,3,4,5,6,7,8},0.01,5);// :-D
	matcher.setNeighbours((int )2);

    Eigen::Affine3d T;
    T.setIdentity();
	for(std::vector<string>::iterator t=pointcloud_files.begin()+1;t<pointcloud_files.end();++t)
	{
		auto t_st = getCloud(*t,5);
		auto t_mv = getCloud(*(t+1),5);
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_st=std::get<0>(t_st);
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_mv=std::get<0>(t_mv);
		std::vector<double> label_st=std::get<1>(t_st);
		std::vector<double> label_mv=std::get<1>(t_mv);

		Eigen::Affine3d Td;
		Td=matcher.match(cloud_st,cloud_mv,{label_st,label_st,label_st,label_st,label_st,label_st,label_st,label_st,label_st},{label_mv,label_mv,label_mv,label_mv,label_mv,label_mv,label_mv,label_mv,label_mv});
        T=Td*T;
        cout<<T.translation().transpose()<<endl;
	}
	return 0;
}
