#include <iostream>
#include <iomanip>
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
template<int M, template<typename> class F = std::less>
struct TupleCompare
{
	template<typename T>
	bool operator()(T const &t1, T const &t2)
	{
		return F<typename tuple_element<M,T>::type>()(std::get<M>(t1), std::get<M>(t2));
	}
};
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
	 ("transforms,t", po::value<std::string>(&transforms_name1), "Transforms file 1, the moving");
	po::parsed_options parsed_options = po::command_line_parser(argc, argv)
		  .options(desc).style(po::command_line_style::unix_style ).run();
    po::variables_map vm;
	po::store(parsed_options, vm);
    po::notify(vm);
	if(vm.count("help"))
	{
		cout<<desc;
		return 0;
	}
	ifstream in_trans1;
	if(vm.count("transforms")) { 
		in_trans1.open(transforms_name1, ifstream::in); 
	}else{
		cout<<"no transform file"<<endl;return -1;
	};


	std::vector<Eigen::Affine3d> transforms1 = readTransform(in_trans1);
	std::vector<std::tuple<int,int,double,double,Eigen::Matrix4d> > results; 
	for(int i=0;i<transforms1.size();i++)
	{
		for(int j=i+1;j<transforms1.size();j++)
		{
			Eigen::Matrix4d A = transforms1[i].matrix();
			Eigen::Matrix4d B = transforms1[j].matrix();
			Eigen::Matrix4d T = A.colPivHouseholderQr().solve(B);
			Eigen::Vector3d translation = T.block<3,1>(0,3);
			double dT = translation.dot(translation);
			double dR = acos((T.block<3,3>(0,0).trace()-1)/2);
/*			Eigen::Affine3d dT= transforms1[i]*transforms1[j].inverse();
			double dT_t = pow(dT.translation()(0),2)+pow(dT.translation()(1),2)+pow(dT.translation()(2),2);
			double dR= acos((dT.rotation().matrix().trace()-1)/2);
			*/
			results.push_back(std::make_tuple(i,j,sqrt(dT),dR,T));
		}
	}
	sort(results.begin(),results.end(), TupleCompare<2>());
	Eigen::IOFormat SpaceSep(StreamPrecision,DontAlignCols, " ", " ", "", "", "", "");
	std::cout <<std::setprecision(5) <<std::fixed;
	for(int i=0;i<results.size();i++)
	{
		std::cout<<std::get<0>(results[i])<<" ";
		std::cout<<std::get<1>(results[i])<<" ";
		std::cout<<std::get<2>(results[i])<<" ";
		std::cout<<std::get<3>(results[i])<<" ";
		std::cout<<std::get<4>(results[i]).format(SpaceSep)<<std::endl;
	}

}

