#include <iostream>
#include <sstream>
#include <string>
#include <map>
#include <fstream>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <pcl/filters/crop_box.h>
#include <se_ndt/se_ndt.hpp>
#include <pcl/common/io.h>
#include <pcl/common/transforms.h>
#include <pcl/filters/voxel_grid.h>
#include <ndt_registration/ndt_matcher_d2d.h>
#include <ctime>
#define PI 3.14159265
std::map<int,int> get_label_map(){
    std::map<int,int> lm;
    lm[0 ]= 0;     // "unlabeled"
    lm[1 ]= 0;     // "outlier" mapped to "unlabeled" --------------------------mapped
    lm[10]= -1;     // "car"
    lm[11]= -1;     // "bicycle"
    lm[13]= -1;     // "bus" mapped to "other-vehicle" --------------------------mapped
    lm[15]= -1;     // "motorcycle"
    lm[16]= -1;     // "on-rails" mapped to "other-vehicle" ---------------------mapped
    lm[18]= -1;     // "truck"
    lm[20]= -1;     // "other-vehicle"
    lm[30]= -1;     // "person"
    lm[31]= -1;     // "bicyclist"
    lm[32]= -1;     // "motorcyclist"
    lm[40]= 1;     // "road"
    lm[44]= 1;     // "parking"
    lm[48]= 2;    // "sidewalk"
    lm[49]= 3;    // "other-ground"
    lm[50]= 4;    // "building"
    lm[51]= 5;    // "fence"
    lm[52]= 0;     // "other-structure" mapped to "unlabeled" ------------------mapped
    lm[60]= 6;     // "lane-marking" to "road" ---------------------------------mapped
    lm[70]= 7;    // "vegetation"
    lm[71]= 8;    // "trunk"
    lm[72]= 9;    // "terrain"
    lm[80]= 10;    // "pole"
    lm[81]= 11;    // "traffic-sign"
    lm[99]= -1;     // "other-object" to "unlabeled" ----------------------------mapped
    lm[252]= -1;    // "moving-car"
    lm[253]= -1;    // "moving-bicyclist"
    lm[254]= -1;    // "moving-person"
    lm[255]= -1;    // "moving-motorcyclist"
    lm[256]= -1;    // "moving-on-rails" mapped to "moving-other-vehicle" ------mapped
    lm[257]= -1;    // "moving-bus" mapped to "moving-other-vehicle" -----------mapped
    lm[258]= -1;    // "moving-truck"
    lm[259]= -1;    // "moving-other-vehicle"
    return lm;
}


using namespace std;
namespace po = boost::program_options;
void filt_write (pcl::PointCloud<pcl::PointXYZI>::Ptr cloud, float max_dist)
{
	pcl::PointCloud<pcl::PointXYZI>::Ptr filt(new pcl::PointCloud<pcl::PointXYZI>);
    int npoints=cloud->points.size();
    for(int pn=0;pn<npoints;pn++)
    {
        auto t=cloud->points[pn];
        if(t.intensity>5)
            filt->points.push_back(t);
    }
    std::ofstream outfile("cloud_tmp.csv");
    for(int pn=0;pn<npoints;pn++)
    {
        bool flag=true;
        Eigen::Map<Eigen::Vector3f> v1(cloud->points[pn].data);
        for(auto f = filt->points.begin();f<filt->points.end();++f)
        {
            Eigen::Map<Eigen::Vector3f> v2(f->data);
            if((v1-v2).norm()<max_dist)
            {
                cloud->points[pn].intensity=7;
            }
        }
        outfile<<cloud->points[pn].x<<" ";
        outfile<<cloud->points[pn].y<<" ";
        outfile<<cloud->points[pn].z<<" ";
        outfile<<cloud->points[pn].intensity<<endl;
 
    }
}


pcl::PointCloud<pcl::PointXYZI>::Ptr filter_class_omp(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud, float max_dist)
{
#define n_threads 12
	pcl::PointCloud<pcl::PointXYZI>::Ptr out(new pcl::PointCloud<pcl::PointXYZI>);
	pcl::PointCloud<pcl::PointXYZI>::Ptr filt(new pcl::PointCloud<pcl::PointXYZI>);
    std::vector<pcl::PointCloud<pcl::PointXYZI>::Ptr > out_omp;
    std::vector<pcl::PointCloud<pcl::PointXYZI>::Ptr > filt_omp;
    for(int i=0;i<n_threads;i++)
    {
        out_omp.push_back(pcl::PointCloud<pcl::PointXYZI>::Ptr(new pcl::PointCloud<pcl::PointXYZI>));
        filt_omp.push_back(pcl::PointCloud<pcl::PointXYZI>::Ptr(new pcl::PointCloud<pcl::PointXYZI>));
    }
    int npoints=cloud->points.size();
    #pragma omp parallel num_threads(n_threads)
    {
    int thread_id = omp_get_thread_num();
    #pragma omp for
    for(int pn=0;pn<npoints;pn++)
    {
        auto t=cloud->points[pn];
        if(t.intensity>5)
            filt_omp[thread_id]->points.push_back(t);
    }
    }
    for(int i=0;i<n_threads;i++)
        (*filt)+=(*filt_omp[i]);
    #pragma omp parallel num_threads(n_threads)
    {
    int thread_id = omp_get_thread_num();
    #pragma omp for
    for(int pn=0;pn<npoints;pn++)
    {
        bool flag=true;
        Eigen::Map<Eigen::Vector3f> v1(cloud->points[pn].data);
        for(auto f = filt->points.begin();f<filt->points.end();++f)
        {
            Eigen::Map<Eigen::Vector3f> v2(f->data);
            if((v1-v2).norm()<max_dist)
            {
                flag=false;
                break;
            }
        }
        if(cloud->points[pn].intensity<2)
            cloud->points[pn].intensity=0;

        if(flag)
            out_omp[thread_id]->points.push_back(cloud->points[pn]);
    }
    }
    for(int i=0;i<n_threads;i++)
        (*out)+=(*out_omp[i]);
    return out;
}

std::tuple<pcl::PointCloud<pcl::PointXYZ>::Ptr,std::vector<double> >  getCloud(string filename, int n_useless, bool binary=false)
{
	pcl::PointCloud<pcl::PointXYZ>::Ptr laserCloud(new pcl::PointCloud<pcl::PointXYZ>);
	std::vector<double> labels;
    if(!binary)
    {
        ifstream in(filename); // for example
        if(!in)std::cerr<<"could not open pointcloud file :"<<filename<<std::endl;
        string line;
        float useless;
        while (getline(in, line)){
            pcl::PointXYZ point;
            double label;
            stringstream sin(line);
            sin>>point.x>>point.y>>point.z;
            (*laserCloud).points.push_back(point);
            for(int i=0;i<n_useless;i++)
                sin>>useless;
            sin>>label;
            labels.push_back(label);

        }
    }
    else
    {
        ifstream inFile(filename, ios::in|ios::binary);
        if(inFile.good())
        {
            float in_float[4];
            while(!inFile.eof())
            {
                inFile.read((char *)in_float, 4*sizeof(float)); 
                pcl::PointXYZ point;
                point.x=in_float[0];
                point.y=in_float[1];
                point.z=in_float[2];
                double label=in_float[3];
                (*laserCloud).points.push_back(point);
                labels.push_back(label);
            }
        }
    }
    return std::make_tuple(laserCloud, labels);
}
pcl::PointCloud<pcl::PointXYZI>::Ptr getB(string filenameP, string filenameL, std::map<int,int> &label_map)
{
	pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloud(new pcl::PointCloud<pcl::PointXYZI>);
    ifstream inFileP(filenameP, ios::in|ios::binary);
    ifstream inFileL(filenameL, ios::in|ios::binary);
    if(inFileP.good() && inFileL.good())
    {
        float in_float[4];
        int label;
        while(!inFileP.eof() && !inFileL.eof())
        {
            inFileP.read((char *)in_float, 4*sizeof(float)); 
            pcl::PointXYZI point;
            point.x=in_float[0];
            point.y=in_float[1];
            point.z=in_float[2];
            point.intensity=in_float[3];
            inFileL.read((char *)&label, sizeof(int));
            point.intensity= label_map[label];
            if(point.intensity==-1)
                continue;
            (*laserCloud).points.push_back(point);
        }
    }
    return laserCloud;
}
void transform_interpolated(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud, Eigen::Affine3d& T)
{
    Eigen::AngleAxisd rotation(T.rotation());
    Eigen::Vector3d translation(T.translation());
    for(int i=0;i<cloud->points.size();i++)
    {
        float theta=atan2(cloud->points[i].y,cloud->points[i].x);
        theta+=PI/2;
        if(theta>2*PI) theta=theta-2*PI;
        if(theta<0) theta=theta+2*PI;
        float scale=1-theta/(2*PI);//If oposit direction then not 1-
        Eigen::Vector3d v(cloud->points[i].x, cloud->points[i].y, cloud->points[i].z);
        Eigen::AngleAxisd r= Eigen::AngleAxisd(-rotation.angle()*scale, rotation.axis());
        v=-scale * translation + r*v;
        cloud->points[i].x=v.x();
        cloud->points[i].y=v.y();
        cloud->points[i].z=v.z();
    }
}
std::vector<Eigen::Affine3d>  readTransform(string fname)
{
	ifstream infile(fname);
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



void print_transform(Eigen::Affine3d &T)
{
        fprintf(stdout, "%e %e %e %e %e %e %e %e %e %e %e %e\n", 
                T(0,0),T(0,1),T(0,2),T(0,3),
                T(1,0),T(1,1),T(1,2),T(1,3),
                T(2,0),T(2,1),T(2,2),T(2,3));
        fflush(stdout);
}
int main(int argc, char** argv)
{
	string transforms;
	string p_dir;
    auto label_map=get_label_map();

	po::options_description desc("Allowed options");
	desc.add_options()
	("help,h", "produce help message")
	("visualize,v", "Visualize with opengl")
	("bin,b", "Input is in binary format?")
	("labels,l", po::value<std::vector<string> >()->multitoken(), "Labels")
	("pointclouds,p", po::value<std::vector<string> >()->multitoken(), "Point cloud files");

	po::parsed_options parsed_options = po::command_line_parser(argc, argv)
      .options(desc).style(po::command_line_style::unix_style ).run();
    po::variables_map vm;
	po::store(parsed_options, vm);
    po::notify(vm);
	if(vm.count("help")||!vm.count("pointclouds"))
	{
		cout<<desc;
		return 0;
	}
	vector<string> pointcloud_files;
	pointcloud_files= vm["pointclouds"].as<vector<string> >();
	vector<string> label_files;
	label_files= vm["labels"].as<vector<string> >();
    Eigen::Affine3d calib, calib_inv;
    calib.matrix() << 4.276802385584e-04,-9.999672484946e-01,-8.084491683471e-03,-1.198459927713e-02,-7.210626507497e-03,8.081198471645e-03,-9.999413164504e-01,-5.403984729748e-02,9.999738645903e-01,4.859485810390e-04,-7.206933692422e-03,-2.921968648686e-01,0,0,0,1;
    calib_inv=calib.inverse();
    NDTMatch_SE matcher ({4,0.8},{0,1},{80,80},12,50);
    if(vm.count("visualize"))
        matcher.visualize();

    Eigen::Affine3d T, prnt;
    //matcher.T=calib;
    int s_point=pointcloud_files.size();
	for(int t=0; t<s_point; t++)
	{
		pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_mv = getB(pointcloud_files[t], label_files[t], label_map);
        //transform_interpolated(cloud_mv, Td); // Could be usefull in other datasets
		T=matcher.slam(cloud_mv);
        //prnt=T*calib_inv;
        //print_transform(prnt);
    }
	return 0;
}

