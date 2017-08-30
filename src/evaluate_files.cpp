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
	string transforms_name1,  transforms_name2;
	po::options_description desc("Allowed options");
    desc.add_options()
	("help", "produce help message")
	 ("transforms1", po::value<std::string>(&transforms_name1), "Transforms file 1")
	 ("transforms2", po::value<std::string>(&transforms_name2), "Transforms file 2");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
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
	for(int i=0;i<transforms1.size();i++)
	{
		Eigen::Vector3d vec1=transforms1[i].translation();
		Eigen::Vector3d vec2=transforms2[i].translation();
//		cout<<vec1.transpose()*vec1-2*vec1.transpose()*vec2+vec2.transpose()*vec2<<endl;
		double currr = (vec1.transpose()*vec1-2*vec1.transpose()*vec2+vec2.transpose()*vec2)(0,0);
		cerr<<i<<" :  SE : "<<currr<<" \t RSE : "<<sqrt(currr)<<endl;
		sse+= currr;
		srse+=sqrt(currr);
	}
	sse/=transforms1.size();
	cout<<sqrt(sse)<<"\t"<<srse/transforms1.size()<<endl;
}

