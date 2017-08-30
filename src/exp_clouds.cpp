#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <Eigen/Geometry>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/filters/random_sample.h>
#include <pcl/filters/conditional_removal.h>
#include <pcl/point_cloud.h>

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
std::list<int>  readSemantics(string filename)
{
	ifstream infile(filename); // for example
	std::list<int> semantics;
	string line = "";
	while(getline(infile,line))
	{
		semantics.push_back(stof(line));
	}
    return semantics;
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
	 ("semantics", po::value<std::vector<string> >()->multitoken(), "Semantics files")
	 ("transforms", po::value<std::string>(&transforms_name), "Transforms file")
	 ("out", po::value<std::string>(&out_name), "Output file prefix")
	("skip", "skip point cloud first line")
	("dev", po::value<float>(&dev)->default_value(1.0), "Standard deviation")
	("deva", po::value<float>(&deva)->default_value(1.0), "Standard deviation of angle")
	("num", po::value<int>(&num)->default_value(10), "Number of generated transforms");
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

	if (!vm["semantics"].empty() && (semantics_files= vm["semantics"].as<vector<string> >()).size() > 0) {
	}else {cout<<"semantics read failure";};

	ifstream in_trans;
	if(vm.count("transforms")) { in_trans.open(transforms_name, ifstream::in); }else{cout<<"no transform file"<<endl;return -1;};


	std::vector<Eigen::Affine3d> transforms;
	std::vector<Eigen::Matrix4d> inv_transforms;
	Eigen::Matrix4d I;
	Eigen::Matrix3f I_f;
	I.setIdentity();
	I_f.setIdentity();

	std::default_random_engine generator;

	std::set<std::pair<int,int> > angles;
	float max_range=50;
	int steps=5;
	int total_points=10000;
	pcl::ExtractIndices<pcl::PointXYZI> extract;
	transforms = readTransform(in_trans);
	//	inv_transforms.push_back(I.colPivHouseholderQr().solve(TM));
	for(int p_cl=0;p_cl<pointcloud_files.size();p_cl++)
	{
		pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudIn=getCloud2(pointcloud_files[p_cl],skip);
		std::list<int> semantics= readSemantics(semantics_files[p_cl]);
		for(int t_n=0;t_n<transforms.size();t_n++)
		{
			Eigen::Affine3d T=transforms.at(t_n);
			pcl::RandomSample<pcl::PointXYZI> sample(true);
			string foname="pcl_";
			foname+=(t_n<10?"0":"")+to_string(t_n)+".csv";
			string sename="sem_";
			sename+=(t_n<10?"0":"")+to_string(t_n)+".csv";
			ofstream file_out(foname,ofstream::out|ofstream::app);
			ofstream file_sem_out(sename,ofstream::out|ofstream::app);
			for(int k=1;k<=steps;k++)
			{
				///////////
				//
				//Conditional Removal
				float radius_k=k*max_range/steps+max_range/steps/pointcloud_files.size()*p_cl;
				double val=-pow(radius_k,2)+T.translation().transpose()*T.translation();
				pcl::TfQuadraticXYZComparison<pcl::PointXYZI>::Ptr quadr_comp(new pcl::TfQuadraticXYZComparison<pcl::PointXYZI>());

				Eigen::Vector3f mcvec=-T.translation().cast<float>();
				quadr_comp->setComparisonVector(mcvec);
				quadr_comp->setComparisonMatrix(I_f);
				quadr_comp->setComparisonScalar(val);
				quadr_comp->setComparisonOperator(pcl::ComparisonOps::GT);

				pcl::ConditionOr<pcl::PointXYZI>::Ptr rad_cond (new pcl::ConditionOr<pcl::PointXYZI> ());
				rad_cond->addComparison(quadr_comp);
				pcl::ConditionalRemoval<pcl::PointXYZI> condrem(true);
				condrem.setCondition (rad_cond);
				condrem.setInputCloud (laserCloudIn);
//				condrem.setKeepOrganized(true);
				pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_sampled_A(new pcl::PointCloud<pcl::PointXYZI>);
				pcl::PointIndices::Ptr inliers (new pcl::PointIndices ());	
				condrem.filter(*cloud_sampled_A);
				inliers->indices = (*condrem.getRemovedIndices());
				///////////////
				//
				//
				//Random sampling
				sample.setInputCloud(laserCloudIn);
				sample.setIndices(inliers);
				pcl::PointIndices::Ptr ind_rem(new pcl::PointIndices ());	
				sample.setSample(total_points/steps/pointcloud_files.size());
//				sample.setKeepOrganized(true);
				sample.filter(ind_rem->indices);


//				std::sort(ind_rem->indices.begin(),ind_rem->indices.end());
				std::list<int>::iterator it = semantics.begin();
				int last_one=-1;
				for(int i=0;i<ind_rem->indices.size();i++) 
				{
					int index_c= ind_rem->indices[i];
					if(index_c>=laserCloudIn->points.size())
					{
						cerr<<"out of bounds"<<endl;
						return -1;
					}
//					file_out<<index_c<<" ! "<<radius_k<<" ! "<<inliers->indices.size()<<" ! ";
					file_out<<laserCloudIn->points[index_c].x<<", ";
					file_out<<laserCloudIn->points[index_c].y<<", ";
					file_out<<laserCloudIn->points[index_c].z<<", ";
					file_out<<laserCloudIn->points[index_c].intensity<<endl;
					std::advance(it, index_c-last_one-1);
					last_one=index_c;
					file_sem_out<<(*it)<<endl;
					std::list<int>::iterator delete_this = it;
					std::advance(it,1);
					semantics.erase(delete_this);

				}
				//////
				//Remove from main point cloud
				extract.setInputCloud(laserCloudIn);
				extract.setIndices(ind_rem);
				extract.setNegative(true);
//				extract.setKeepOrganized(true);
				extract.filter(*laserCloudIn);


			}
			file_out.close();
		}


//		angles.push_back(std::pair<int,int>(114.59*a,114.59*b));///0.5 degrees rounding
	}
}

