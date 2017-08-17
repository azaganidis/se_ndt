#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <ndt_registration/ndt_matcher_d2d.h>

using namespace std;
namespace po = boost::program_options;
Eigen::Affine3d getTrans(string filename)
{
	ifstream infile(filename); // for example
	string line = "";
	Eigen::Affine3d T;
	T.setIdentity();
	for(int i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			string word="";
			char c=' ';
			while(c!='.'&&c!='-'&&(c<'0'||c>'9'))
			{
				c=infile.get();
			}
			while(c=='e'||c=='.'||c=='-'||(c>='0'&&c<='9'))
			{
				word+=c;
				c=infile.get();
			}
			T(i,j)=stof(word);
		}
	}
    return T;
}
int main(int argc, char** argv)
{
	po::options_description desc("Allowed options");
    desc.add_options()
	("help", "produce help message")
	 ("trans", po::value<std::vector<string> >()->multitoken(), "transform files");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
	vector<string> trans_files;
	if (!vm["trans"].empty() && (trans_files= vm["trans"].as<vector<string> >()).size() >= 2) {
		///cout<<"success pointcloud read";
	}else {cout<<"transform read failure";};
	int num_files=trans_files.size();
	Eigen::Affine3d T;
	T.setIdentity();
	for(int i=0;i<num_files;i++)
	{
		Eigen::Affine3d T_pred;
		T_pred=getTrans(trans_files[i]);
		T=T*T_pred;
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				cout<<T(i,j)<<", ";
		cout<<endl;
	}

	return 0;
}

