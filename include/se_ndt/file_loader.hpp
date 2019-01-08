#include <pcl/point_cloud.h>
#include <pcl/common/io.h>
using namespace std;
template <>
pcl::PointCloud<pcl::PointXYZI>::Ptr getCloud<pcl::PointXYZI>(string filename,char IFS, bool skip)
{
	ifstream infile(filename); // for example
	string line = "";
	pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloud(new pcl::PointCloud<pcl::PointXYZI>);
	if(skip)getline(infile, line);
	while (getline(infile, line)){
		stringstream strstr(line);
		string word = "";
		pcl::PointXYZI point;
		getline(strstr,word, IFS);
		point.x=stof(word);
		getline(strstr,word, IFS);
		point.y=stof(word);
		getline(strstr,word, IFS);
		point.z=stof(word);
		getline(strstr,word);
		point.intensity=stof(word);
		(*laserCloud).points.push_back(point);
	}
    return laserCloud;
}
template <>
pcl::PointCloud<pcl::PointXYZ>::Ptr getCloud<pcl::PointXYZ>(string filename,char IFS, bool skip)
{
	ifstream infile(filename); // for example
	string line = "";
	pcl::PointCloud<pcl::PointXYZ>::Ptr laserCloud(new pcl::PointCloud<pcl::PointXYZ>);
	if(skip)getline(infile, line);
	while (getline(infile, line)){
		stringstream strstr(line);
		string word = "";
		pcl::PointXYZ point;
		getline(strstr,word, IFS);
		point.x=stof(word);
		getline(strstr,word, IFS);
		point.y=stof(word);
		getline(strstr,word, IFS);
		point.z=stof(word);
		(*laserCloud).points.push_back(point);
	}
    return laserCloud;
}
Eigen::Affine3d readTransform(istream &infile)
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
