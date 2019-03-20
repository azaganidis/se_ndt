#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <ndt_registration/ndt_matcher_d2d.h>

using namespace std;
Eigen::Affine3d getTrans(istream &infile, Eigen::Affine3d T, bool invert=false)
{
	Eigen::Affine3d Td;
	Td.setIdentity();
	for(int i=0;i<4;i++)
		for(int j=0;j<4;j++)
			infile>>Td(i,j);
    if(invert)
        T=Td*T;
    else
        T=T*Td;

    for(int k=0;k<4;k++)
        for(int j=0;j<4;j++)
            cout<<T(k,j)<<" ";
    cout<<endl;
    return T;
}
int main(int argc, char** argv)
{
	bool is_pipe=false;
	bool invert=false;
    if(argc==1)
        std::cerr<<"Usage: concat_transform -i file1 file2 file3.\n -i for invert, -p for pipe in"<<std::endl;
    for(int i=1;i<argc;i++)
        if(argv[i][0]=='-')
            switch (argv[i][1])
            {
                case 'i':
                    invert = true;
                    break;
                case 'p':
                    is_pipe= true;
                    break;
                default:
                    std::cerr<<"Usage: concat_transform -i file1 file2 file3.\n -i for invert, -p for pipe in"<<std::endl;
                    break;
            }
	Eigen::Affine3d T;
	T.setIdentity();
	for(int i=1;i<argc||is_pipe;i++)
	{
		if(!is_pipe&&argv[i][0]!='-')
		{
			ifstream in_file(argv[i]);
            while(!in_file.eof())
            {
                T=getTrans(in_file, T, invert);
            }
		}
		else
		{
            while(!cin.eof())
                T=getTrans(cin, T, invert);
		}
	}
	return 0;
}

