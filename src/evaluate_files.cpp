#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <Eigen/Geometry>

using namespace std;
namespace po = boost::program_options;
using namespace Eigen;
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
	bool inverse_1=false;
	bool inverse_2=false;
	string transforms_name1,  transforms_name2;
	po::options_description desc("Allowed options");
    desc.add_options()
	("help", "If the transforms are not inverse, but show the same direction, add -i.")
	 ("transforms1", po::value<std::string>(&transforms_name1), "Transforms file 1, the moving")
	 ("inv1,i", "Invert file 1")
	 ("inv2,j", "Invert file 2")
	 ("transforms2", po::value<std::string>(&transforms_name2), "Transforms file 2, the reference");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
	if(vm.count("inv1"))
		inverse_1=true;
	if(vm.count("inv2"))
		inverse_2=true;
	if(vm.count("help"))
	{
		cout<<desc;
		return 0;
	}

	ifstream in_trans1;
	ifstream in_trans2;
	if(vm.count("transforms1")) { in_trans1.open(transforms_name1, ifstream::in); }else{cout<<"no transform file"<<endl;return -1;};
	if(vm.count("transforms2")) { in_trans2.open(transforms_name2, ifstream::in); }else{cout<<"no transform file"<<endl;return -1;};


	Eigen::Matrix4d I;
	Eigen::Matrix3f I_f;
	I.setIdentity();
	I_f.setIdentity();

	std::vector<Eigen::Affine3d> transforms1 = readTransform(in_trans1);
	std::vector<Eigen::Affine3d> transforms2 = readTransform(in_trans2);
	double sse=0;
	double srse=0;
	double res_sum=0;
	for(int i=0;i<transforms1.size();i++)
	{
		if(inverse_1)
			transforms1[i]=transforms1[i].inverse();
		if(inverse_2)
			transforms2[i]=transforms2[i].inverse();
		Eigen::Vector3d vec1=transforms1[i].translation();
		Eigen::Vector3d vec2=transforms2[i].translation();

		Eigen::Affine3d dT_i= transforms2[i]*transforms2[0].inverse();
		double dT_it = sqrt(pow(dT_i.translation()(0),2)+pow(dT_i.translation()(1),2)+pow(dT_i.translation()(2),2));
		double dR_i= acos((dT_i.rotation().matrix().trace()-1)/2);

		Eigen::Affine3d dT= transforms1[i]*transforms2[i].inverse();
		double dT_t = pow(dT.translation()(0),2)+pow(dT.translation()(1),2)+pow(dT.translation()(2),2);
		double dR= acos((dT.rotation().matrix().trace()-1)/2);

		cerr<<dT_it<<", "<<dR_i<<", "<<i<<" :  SE : "<<dT_t<<" \t RSE : "<<sqrt(dT_t)<<"\t : "<<dR<<endl;
		cout<<dT_it<<", "<<dR_i<<", "<<sqrt(dT_t)<<", "<<dR<<endl;
		if(!isnan(dR))res_sum+=dR;
		sse+= dT_t;
		srse+=sqrt(dT_t);
	}
	sse/=transforms1.size();
	cerr<<sqrt(sse)<<"\t"<<srse/transforms1.size()<<"\t "<<res_sum/transforms1.size()<<endl;
}

