#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <pcl/filters/crop_box.h>
#include <se_ndt/se_ndt.hpp>
#include <pcl/common/io.h>
#include <ndt_registration/ndt_matcher_d2d.h>

using namespace std;
namespace po = boost::program_options;
std::tuple<pcl::PointCloud<pcl::PointXYZ>::Ptr,std::vector<double> >  getCloud(string filename, int n_useless)
{
	pcl::PointCloud<pcl::PointXYZ>::Ptr laserCloud(new pcl::PointCloud<pcl::PointXYZ>);
	std::vector<double> labels;
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
    return std::make_tuple(laserCloud, labels);
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
		Eigen::Vector3d translation = T.block<3,1>(0,3);
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
	for(int i=0;i<transforms1.size();i++)
	{
		for(int j=i+1;j<transforms1.size();j++)
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

int main(int argc, char** argv)
{
	string transforms;
	float p;
	string p_dir;

	po::options_description desc("Allowed options");
	desc.add_options()
	("help", "produce help message")
	 ("transforms,t", po::value<std::string >(&transforms), "File with the transforms")
	 ("parameter,p", po::value<float >(&p), "var")
	 ("ground,g", "Use ground truth labels")
	 ("pointcloud_dir,d", po::value<string >(&p_dir), "Pointcloud file directory. Point X Y Z Class");

	po::parsed_options parsed_options = po::command_line_parser(argc, argv)
      .options(desc).style(po::command_line_style::unix_style ).run();
    po::variables_map vm;
	po::store(parsed_options, vm);
    po::notify(vm);
	if(vm.count("help")||!vm.count("pointcloud_dir")||!vm.count("transforms"))
	{
		cout<<desc;
		return 0;
	}
	bool ground_labels=false;
	if(vm.count("ground"))
			ground_labels=true;
	std::vector<std::tuple<int,int,double,double,Eigen::Matrix4d> >  test_on = getThat(transforms,true);

	NDTMatch_SE matcher ({40,16,2,0.8},{0,1,2,3},{100,100,100},{'=','=','=','=','=','=','=','=','='},{0,1,2,3,4,5,6,7,8},0.01,1);// :-D
//	NDTMatch_SE matcher ({1,2},{0,1,0},{100,100,100},{'=','=','=','=','=','=','=','=','='},{0,1,2,3,4,5,6,7,8},0.01,5);// :-D

	//NDTMatch_SE matcher ({1,6,10,20,30,60},{5,4,3,2,0,1,0},{100,100,30},{'*'},{0},0.01,5);// :-D
	//NDTMatch_SE matcher ({1,6,10,20,30,60},{5,4,3,2,0,1,0},{100,100,30},{'*'},{0},0.01,5);// :-D
	matcher.setNeighbours((int )1);

	std::vector<std::tuple<int,int,double,double,Eigen::Matrix4d> >::iterator it; 
	double trans_tot=0;
	double ang_tot=0;
	double trans_tot_r=0;
	double ang_tot_r=0;
	int num=0;

	for(it = test_on.begin();it<test_on.end();it++)
	{
		std::string st_cloud = p_dir+std::to_string(std::get<0>(*it))+"_out.csv";
		std::string mv_cloud = p_dir+std::to_string(std::get<1>(*it))+"_out.csv";
		Eigen::Matrix4d A = std::get<4>(*it);
		Eigen::Matrix4d m, Iden;
		m.setIdentity();
		Iden.setIdentity();
		auto t_st = getCloud(st_cloud,ground_labels?4:5);
		auto t_mv = getCloud(mv_cloud,ground_labels?4:5);
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_st=std::get<0>(t_st);
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_mv=std::get<0>(t_mv);
		std::vector<double> label_st=std::get<1>(t_st);
		std::vector<double> label_mv=std::get<1>(t_mv);

		/*
		cout<<st_cloud <<" "<<mv_cloud<<std::endl;
		int a=0;
		while(ang_diff(Iden,m*A)>0.3)
		{
			a++;
			m.block<3,3>(0,0) = Eigen::AngleAxisd(a*0.3,Eigen::Vector3d::UnitZ()).matrix();
		}
		A=m*A;
		pcl::transformPointCloud(*cloud_mv, *cloud_mv, m);
		ofstream mco("cloud_mv.csv");
		for( int i=0;i<cloud_mv->points.size();i++)
		{
			mco<<cloud_mv->points[i].x<<" ";
			mco<<cloud_mv->points[i].y<<" ";
			mco<<cloud_mv->points[i].z<<" ";
			mco<<label_mv[i]<<endl;
		}
		ofstream sco("cloud_st.csv");
		for( int i=0;i<cloud_st->points.size();i++)
		{
			sco<<cloud_st->points[i].x<<" ";
			sco<<cloud_st->points[i].y<<" ";
			sco<<cloud_st->points[i].z<<" ";
			sco<<label_st[i]<<endl;
		}
		mco.close();
		sco.close();
		std::cout<<A<<std::endl;
		
		Eigen::Affine3d Tndt;
		Tndt.setIdentity();
		m_ndt.match(*cloud_mv, *cloud_st,Tndt,false); 
		Eigen::Matrix4d B=Tndt.matrix();
		*/

		double initTE=pos_diff(Iden,A);
		double initRE=ang_diff(Iden,A);
		Eigen::Matrix4d B;

		B=matcher.match(cloud_st,cloud_mv,{label_st,label_st,label_st,label_st,label_st,label_st,label_st,label_st,label_st},{label_mv,label_mv,label_mv,label_mv,label_mv,label_mv,label_mv,label_mv,label_mv}).matrix();
	//	Eigen::Matrix4d B=matcher.match(cloud_st,cloud_mv,{label_st},{label_mv}).matrix();

		Eigen::Matrix4d T = A.colPivHouseholderQr().solve(B);
		Eigen::Vector3d translation = T.block<3,1>(0,3);
		double dT = translation.dot(translation);
		double t_m_d=(T.block<3,3>(0,0).trace()-1.0)/2.0;
		if(t_m_d>1)t_m_d=1;
		double dR = acos(t_m_d);
		std::cout<<sqrt(initTE)<<" "<<initRE<<" "<<sqrt(dT)<<" "<<dR<<std::endl;
		num++;
		trans_tot+=dT;
		ang_tot+=dR;
		trans_tot_r+=initTE;
		ang_tot_r+=initRE;
	}
	std::cerr<<sqrt(trans_tot_r/num)<<" "<<ang_tot_r/num<<" ";
	std::cerr<<sqrt(trans_tot/num)<<" "<<ang_tot/num<<std::endl;
	return 0;
}

		/*
	}
	Eigen::Affine3d T;
	Eigen::Affine3d Tt;
	T.setIdentity();
	Tt.setIdentity();
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_ref=getCloud(pointcloud_files[0]);
	std::vector<double> rsd_ref=getMeasure(rsd_min_files[0]);
	for(int i=0;i<num_files;i++)
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1=getCloud(pointcloud_files[i]);
		if(trans){
			Eigen::Affine3d T=readTransform(in_trans);
			pcl::transformPointCloud(*cloud1,*cloud1,T);
		}
		std::vector<double> rsd_min=getMeasure(rsd_min_files[i]);
					this simulates dynamic classes
		for(unsigned int j=0;j<cloud1->size();j++)
			if(rsd_min[j]==3)
			{
				Eigen::Vector3f p;
				p[0]=cloud1->points[j].x;
				p[1]=cloud1->points[j].y;
				p[2]=cloud1->points[j].z;
				p=Eigen::AngleAxisf(0.25*M_PI, Eigen::Vector3f::UnitX())*p+Eigen::Vector3f::UnitX();
				cloud1->points[j].x=p[0];
				cloud1->points[j].y=p[1];
				cloud1->points[j].z=p[2];
			}
			
		Tt=matcher.match(cloud1,cloud_ref,{rsd_min,rsd_min,rsd_min,rsd_min,rsd_min,rsd_min,rsd_min,rsd_min},{rsd_ref,rsd_ref,rsd_ref,rsd_ref,rsd_ref,rsd_ref,rsd_ref,rsd_ref});
		T=T*Tt;

		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				cout<<(conc?T(i,j):Tt(i,j))<<", ";
		cout<<endl;
	}

	return 0;
}
*/
