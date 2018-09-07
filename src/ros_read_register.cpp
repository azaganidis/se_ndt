#include <iostream>
#include <sstream>
#include <string>
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

#include <ros/ros.h>
#include <pcl_ros/point_cloud.h>
#include <pcl/point_types.h>
#include <tf2_ros/transform_broadcaster.h>
#include <geometry_msgs/TransformStamped.h>
#include <tf2_eigen/tf2_eigen.h>
#include <pcl_conversions/pcl_conversions.h>
using namespace std;
namespace po = boost::program_options;
void send_transform(Eigen::Affine3d &T)
{
    static tf2_ros::TransformBroadcaster br;
    geometry_msgs::TransformStamped tS=tf2::eigenToTransform(T);
    tS.header.stamp=ros::Time::now();
    tS.header.frame_id="world";
    tS.child_frame_id="velodyne";
    br.sendTransform(tS);
}
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
pcl::PointCloud<pcl::PointXYZI>::Ptr getB(string filename)
{
	pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloud(new pcl::PointCloud<pcl::PointXYZI>);
    ifstream inFile(filename, ios::in|ios::binary);
    if(inFile.good())
    {
        float in_float[4];
        while(!inFile.eof())
        {
            inFile.read((char *)in_float, 4*sizeof(float)); 
            pcl::PointXYZI point;
            point.x=in_float[0];
            point.y=in_float[1];
            point.z=in_float[2];
            point.intensity=in_float[3];
            (*laserCloud).points.push_back(point);
        }
    }
    return laserCloud;
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
double ang_diff(Eigen::Matrix4d A, Eigen::Matrix4d B)
{
		Eigen::Matrix4d T = A.colPivHouseholderQr().solve(B);
		return  acos((T.block<3,3>(0,0).trace()-1)/2);
}
double pos_diff(Eigen::Matrix4d A, Eigen::Matrix4d B)
{
		Eigen::Matrix4d T = A.colPivHouseholderQr().solve(B);
		Eigen::Vector3d translation = T.block<3,1>(0,3);
		return  translation.dot(translation);
}
std::vector<std::tuple<int,int,double,double,Eigen::Matrix4d> > getThat(string in_trans1, bool inverse_file=false)
{
	std::vector<Eigen::Affine3d> transforms1 = readTransform(in_trans1);
	std::vector<std::tuple<int,int,double,double,Eigen::Matrix4d> > results; 
	for(unsigned int i=0;i<transforms1.size();i++)
	{
		for(unsigned int j=i+1;j<transforms1.size();j++)
		{
			Eigen::Matrix4d A,B;
			if(inverse_file)
			{
				A = transforms1[i].inverse().matrix();
				B = transforms1[j].inverse().matrix();
			}
			else
			{
				A = transforms1[i].matrix();
				B = transforms1[j].matrix();
			}
			Eigen::Matrix4d T= A.colPivHouseholderQr().solve(B);
			Eigen::Vector3d translation = T.block<3,1>(0,3);
			double dT = translation.dot(translation);
			double dR = acos((T.block<3,3>(0,0).trace()-1)/2);
			results.push_back(std::make_tuple(i,j,sqrt(dT),dR,T));

		}
	}
	sort(results.begin(),results.end(), TupleCompare<2>());
	double max_distance=3;
	int i_max=50;
	std::vector<std::tuple<int,int,double,double,Eigen::Matrix4d> > results_distributed; 
	std::vector<std::tuple<int,int,double,double,Eigen::Matrix4d> >::iterator r_iterator=results.begin(); 
	results_distributed.push_back(*r_iterator);
	float distance_step=max_distance/i_max;
	while(std::get<2>(results_distributed.back())<max_distance)
	{
		while(std::get<2>(*r_iterator)<distance_step+std::get<2>(results_distributed.back())&&r_iterator!=std::prev(results.end()))
			r_iterator++;
		results_distributed.push_back(*r_iterator);
	}
	return results_distributed;
}
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
    float val=0;
    //bool bin_in=false;
    string topic="velodyne_points", posesF;

	po::options_description desc("Allowed options");
	desc.add_options()
	("help,h", "produce help message")
	//("bin,b", "Input is in binary format?")
    ("topic,t", po::value<string>(&topic),"Topic to publish")
	("value,v", po::value<float >(&val), "Point cloud files")
	("poses,g", po::value<string >(&posesF), "Ground truth poses.")
	("pointclouds,p", po::value<std::vector<string> >()->multitoken(), "Point cloud files");

	po::parsed_options parsed_options = po::command_line_parser(argc, argv)
      .options(desc).style(po::command_line_style::unix_style ).run();
    po::variables_map vm;
	po::store(parsed_options, vm);
    po::notify(vm);
	//if(vm.count("bin"))
     //   bin_in = true;
	if(vm.count("help")||!vm.count("pointclouds"))
	{
		cout<<desc;
		return 0;
	}
	vector<string> pointcloud_files;
	pointcloud_files= vm["pointclouds"].as<vector<string> >();
    Eigen::Affine3d calib;
    ifstream poses;
    if(vm.count("poses"))
        poses.open(posesF);

    calib.matrix() << 4.276802385584e-04,-9.999672484946e-01,-8.084491683471e-03,-1.198459927713e-02,-7.210626507497e-03,8.081198471645e-03,-9.999413164504e-01,-5.403984729748e-02,9.999738645903e-01,4.859485810390e-04,-7.206933692422e-03,-2.921968648686e-01,0,0,0,1;
    cerr<<calib.matrix()<<endl;

    Eigen::Affine3d T;
    T.setIdentity();
//    T=calib.inverse();
    int s_point=pointcloud_files.size();
		Eigen::Affine3d Td;
        Td.setIdentity();

    pcl::VoxelGrid<pcl::PointXYZI > sor;
    float filter =0.8f;
//    filter=val;
    sor.setLeafSize(filter,filter,filter);



	ros::init (argc,argv,"pub_sendt");
	ros::NodeHandle nh;
	ros::Publisher pub = nh.advertise<pcl::PointCloud<pcl::PointXYZI> >(topic, 10);

	NDTMatch_SE matcher ({4},{0},{100,100,100},{'=','=','=','=','=','=','=','='},{0,1,2,3,4,5,6,7},0.01,500);// :-D
	NDTMatch_SE matcher2 ({0.8},{0},{20,20,20},{'=','=','=','=','=','=','=','='},{0,1,2,3,4,5,6,7},0.01,1);// :-D
	NDTMatch_SE matcherNDT ({0.8},{0},{20,20,20},{'*'},{0},0.01,1);// :-D
	matcher.setNeighbours((int )2);
	matcherNDT.setNeighbours((int )1);
	matcher2.setNeighbours((int )1);
    ros::Rate loop_rate(val);
	for(int t=0; t<s_point&&ros::ok(); t++)
	{
		pcl::PointCloud<pcl::PointXYZI>::Ptr cloud_mv = getB(pointcloud_files[t]);
        auto r_time = ros::Time::now();
        clock_t begin_time = clock();
		Td=matcher.matchFaster(Td,cloud_mv);
//		Td=matcher2.matchFaster(Td,cloud_mv);
        //pcl::PointCloud<pcl::PointXYZ>::Ptr cl(new pcl::PointCloud<pcl::PointXYZ>);
        //pcl::copyPointCloud(*cloud_mv, *cl);
		//Td=matcherNDT.match(Td,cl, {std::vector<double>()});
        cout<<"1 "<<float( clock() -begin_time ) / CLOCKS_PER_SEC<<endl;begin_time=clock();
        T=T*Td;
        pcl_conversions::toPCL(r_time,cloud_mv->header.stamp);
        cloud_mv->header.frame_id="velodyne";
        Eigen::Affine3d ts=T;
        if(vm.count("poses"))
        {
            Eigen::Matrix<double,4,4> t_Temp;
            t_Temp.setIdentity();
            for(int i=0;i<3;i++)
                for(int j=0;j<4;j++)
                    poses>>t_Temp(i,j);
            ts.matrix()=calib.matrix().inverse()*t_Temp;
            pcl::transformPointCloud(*cloud_mv, *cloud_mv, calib);
        }
        pub.publish(cloud_mv);
        print_transform(ts);
        send_transform(ts);
        ros::spinOnce();
        loop_rate.sleep();
    }
	return 0;
}
