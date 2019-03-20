#include <iostream>
#include <random>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <Eigen/Geometry>

using namespace std;
namespace po = boost::program_options;
using namespace Eigen;
int main(int argc, char** argv)
{
	int num;
	float dev;
	float deva;
	float devz;
	float seed_;
	float offset_x,offset_y, offset_z;
	po::options_description desc("Allowed options");
    desc.add_options()
	("help", "")
	("uniform", "Uniform distribution for X, Y. If defined, the option --dev is the limits (-,+)")
	("dev", po::value<float>(&dev)->default_value(1.0), "Standard deviation of X,Y.")
	("devz", po::value<float>(&devz)->default_value(0.1), "Standard deviation of z")
	("deva", po::value<float>(&deva)->default_value(0.1), "Standard deviation of angle")
	("offset-x", po::value<float>(&offset_x)->default_value(0.0), "Translation offset x")
	("offset-y", po::value<float>(&offset_y)->default_value(0.0), "Translation offset y")
	("seed", po::value<float>(&seed_)->default_value(0.0), "Seed")
	("offset-z", po::value<float>(&offset_z)->default_value(0.0), "Translation offset z")
	("num", po::value<int>(&num)->default_value(100), "Number of generated transforms");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
	if(vm.count("help"))
	{
		cout<<desc;
		return 0;
	}
	Eigen::Affine3d T;
	T.setIdentity();
	  std::default_random_engine generator;
	if(vm.count("seed"))
		generator.seed(seed_);
	std::normal_distribution<double> distr(0,dev);
	std::normal_distribution<double> distrz(0,devz);
	std::normal_distribution<double> distrax(0,deva);
	std::uniform_real_distribution<double> distraz(-M_PI,M_PI);
	std::uniform_real_distribution<double> distru(-dev,dev);
	T.translation()=Eigen::Vector3d(offset_x,offset_y,offset_z);
	for(int cr=0;cr<num;cr++)
	{
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				cout<<T(i,j)<<", ";
		cout<<endl;
		T.setIdentity();
		T=	AngleAxisd(distrax(generator),Eigen::Vector3d::UnitX())*
			AngleAxisd(distraz(generator),Eigen::Vector3d::UnitZ())*T;
		if(vm.count("uniform"))
		T.translation()=Eigen::Vector3d(offset_x+distru(generator),offset_y+distru(generator),offset_z+distrz(generator));
		else
		T.translation()=Eigen::Vector3d(offset_x+distr(generator),offset_y+distr(generator),offset_z+distrz(generator));

	}
}

