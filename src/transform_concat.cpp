#include <iostream>
#include <string>
#include <fstream>
#include <random>
#include <unistd.h>
#include <Eigen/Geometry>

using namespace std;
void printTransform(Eigen::Affine3d &T, ostream &outS)
{
	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
			outS<<T(i,j)<<", ";
	outS<<endl;
}
void readTransform(istream &infile)
{
	Eigen::Affine3d T0;
	T0.setIdentity();
	bool f_read=false;
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
		T0=T0*T;
		printTransform(T0,cout);
		
	}
}
int main(int argc, char** argv)
{
	std::istream* in_trans=&std::cin;
	std::ifstream in_file_trans;
	if(argc==2&&argv[1]!="-")
	{
		in_file_trans.open(argv[1], ifstream::in);
		in_trans=&in_file_trans;
	}
	readTransform(*in_trans);
}

