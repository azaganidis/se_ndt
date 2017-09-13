#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <Eigen/Geometry>

using namespace std;
Eigen::Affine3d getTrans(istream &infile)
{
	string line = "";
	Eigen::Affine3d T;
	T.setIdentity();
	while(getline(infile, line))
	{
		stringstream str(line);
		string word="";
		getline(str,word,',');
		T.translation()[0]=stof(word);
		getline(str,word,',');
		T.translation()[1]=stof(word);
		getline(str,word,',');
		T.translation()[2]=stof(word);
		getline(str,word,',');
		Eigen::AngleAxisd rollA(M_PI*stof(word) /180, Eigen::Vector3d::UnitY());
		getline(str,word,',');
		Eigen::AngleAxisd pitchA(M_PI*stof(word)/180 , Eigen::Vector3d::UnitX());
		getline(str,word,',');
		Eigen::AngleAxisd yawA(M_PI*stof(word)/180 , Eigen::Vector3d::UnitZ());
		Eigen::Quaternion<double> q = rollA* yawA* pitchA;
		T=T*q;
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				cout<<T(i,j)<<", ";
		cout<<endl;
	}
    return T;
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
	getTrans(*in_trans);
	return 0;
}

