#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <ndt_registration/ndt_matcher_d2d.h>

using namespace std;
namespace po = boost::program_options;
Eigen::Affine3d getTrans(istream &infile)
{
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
	("invert", "apply Td*T instead of T*Td")
	 ("trans", po::value<std::vector<string> >()->multitoken(), "Transform files");
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);
	vector<string> trans_files;
	int is_pipe=1;
	bool invert=false;
	if(vm.count("invert"))
		invert=true;
	if(vm.count("help"))
	{
		cout<<desc;
		cout<<"use a pipe if only one file with transforms"<<endl;
		return 0;
	}
	if (!vm["trans"].empty() && (trans_files= vm["trans"].as<vector<string> >()).size() >= 2) {
		is_pipe=0;
		///cout<<"success pointcloud read";
	}
	int num_files=trans_files.size();
	Eigen::Affine3d T;
	T.setIdentity();
	for(int i=0;i<num_files+is_pipe;i++)
	{
		Eigen::Affine3d T_pred;
		if(!is_pipe)
		{
			ifstream in_file(trans_files[i]);
			T_pred=getTrans(in_file);
		}
		else
		{
			T_pred=getTrans(cin);
			if(cin.eof())
				return 0;
			i--;
		}
		if(invert)
			T=T_pred*T;
		else
			T=T*T_pred;
		for(int k=0;k<4;k++)
			for(int j=0;j<4;j++)
				cout<<T(k,j)<<", ";
		cout<<endl;
	}
	return 0;
}

