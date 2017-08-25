#include <iostream>
#include <random>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <Eigen/Geometry>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/point_cloud.h>

using namespace std;
namespace po = boost::program_options;
using namespace Eigen;
int main(int argc, char** argv)
{
	int num;
	float dev;
	float deva;
	string pointcloud_in_name, transforms_name;
	po::options_description desc("Allowed options");
    desc.add_options()
	("help", "produce help message")
	 ("pointcloud", po::value<std::string>(&pointcloud_in_name), "Point cloud file")
	 ("transforms", po::value<std::string>(&transforms_name), "Transforms file")
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

	ifstream in_trans;
	if(vm.count("transforms")) { in_trans.open(transforms, ifstream::in); trans=true; }else{cout<<"no transform file"<<endl;return -1;};

	pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudIn=getCloud2(pointcloud_in_name,skip);

	std::vector<Eigen::Affine3d> transforms;
	std::vector<Eigen::Affine3d> inv_transforms;
	Eigen::Affine3d A;
	A.setIdentity();

	std::default_random_engine generator;

	std::set<std::pair<int,int> > angles;
	float radius_k=5;
	int points_per_radius=200;
	  pcl::ExtractIndices<pcl::PointXYZ> extract;
	while(in_trans.peek() != std::ifstream::traits_type::eof())
	{
		Eigen::Affine3d T= = readTransform(in_trans);
		transforms.push_back(T);
		inv_transforms.push_back(A.colPivHouseholderQr().solve(T));
		Eigen::Vector3d vec=T.translation();

		pcl::RandomSample<PointXYZI> sample(true);
		pcl::PointCloud<PointXYZI>::Ptr cloud_final(new pcl::PointCloud<PointXYZI>);
		for(int k=0;k<20;k++)
		{
			pcl::ConditionOr<pcl::PointXYZI>::Ptr rad_cond (new pcl::ConditionOr<pcl::PointXYZ> ());
			rad_cond->addComparison(pcl::TfQuadraticXYZComparison<pcl::PointXYZI>::ConstPtr (new 
					pcl::TfQuadraticXYZComparison<pcl::PointXYZI>
					(pcl::ComparisonOps::LT
					 ,A.rotation.matrix()
					 ,-T.translation()
					 ,-pow(radius_k*k,2)*T.translation().transpose*T.translation())));
			// build the filter
			pcl::ConditionalRemoval<pcl::PointXYZ> condrem;
			condrem.setCondition (rad_cond);
			condrem.setInputCloud (laserCloudIn);
			condrem.setKeepOrganized(true);
			// apply filter
			pcl::PointIndices::Ptr inliers (new pcl::PointIndices ());	
			pcl::PointCloud<PointXYZI>::Ptr cloud_filtered(new pcl::PointCloud<PointXYZI>);
			condrem.setNegative(true);
			condrem.filter (*cloud_filtered);
			inliers = condrem.getRemovedIndices();
			pcl::PointCloud<PointXYZI>::Ptr cloud_sampled(new pcl::PointCloud<PointXYZI>);
			sample.setInputCloud(laserCloudIn);
			sample.setIndices(inliers);
			std::vector<int> ind_rem;
			sample.setSample(points_per_radius);
			sample.filter(ind_rem);
			sample.filter(cloud_sampled);
			extract.setInputCloud(laserCloudIn);
			extract.setIndices(ind_rem);
			extract.filter(cloud_sampled);

		}


		angles.push_back(std::pair<int,int>(114.59*a,114.59*b));///0.5 degrees rounding
	}
	int cloudSize = laserCloudIn->points.size();
	for(int i=0;i<cloudSize;i++)
	{

	}
		pcl::transformPointCloud(*cloud1,*cloud1,T);

	Eigen::Affine3d T;
	T.setIdentity();
	  std::default_random_engine generator;
	std::normal_distribution<double> distr(0,dev);
	std::normal_distribution<double> distra(0,deva);
	for(int cr=0;cr<num;cr++)
	{
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				cout<<T(i,j)<<", ";
		cout<<endl;
		T.setIdentity();
		T=T*AngleAxisd(distra(generator),Vector3d(distr(generator),distr(generator),distr(generator)).normalized());
		T.translation()=Eigen::Vector3d(distr(generator),distr(generator),distr(generator));

	}
}

