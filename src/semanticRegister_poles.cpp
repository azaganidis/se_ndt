#define BULK_ADJUSTMENT 1

#ifdef BULK_ADJUSTMENT
#include "g2o/core/sparse_optimizer.h"
#include "g2o/types/slam3d/types_slam3d.h"
#include "g2o/core/optimization_algorithm_factory.h"
using namespace g2o;
#endif

#include <iostream>
#include <sstream>
#include <string>
#include <fstream>
#include <unistd.h>
#include <boost/program_options.hpp>
#include <pcl/filters/crop_box.h>
#include <se_ndt/ndt_fuser_hmt_se.h>

#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
using namespace std;
namespace po = boost::program_options;
#define G2O_OPT 1

Eigen::VectorXd get7d(Eigen::Affine3d p)
{
	Eigen::VectorXd v(7);
	Eigen::Quaternion<double> q(p.rotation());
	v<<p.translation(),q.coeffs();
	v(6)=sqrt(1-v.segment(3,3).norm());
	return v;
}
vector<double>  getMeasure(string filename)
{
	ifstream infile(filename); // for example
	string line = "";
	vector<double> measure;
	while (getline(infile, line)){
		measure.push_back(stod(line));
	}
    return measure;
}
pcl::PointCloud<pcl::PointXYZI>::Ptr cropIt(pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloud, std::vector<float>& b)
{
	pcl::PointCloud<pcl::PointXYZI>::Ptr outC(new pcl::PointCloud<pcl::PointXYZI>);
	pcl::CropBox<pcl::PointXYZI> cropBoxFilter(false);
	cropBoxFilter.setInputCloud(laserCloud);
	cropBoxFilter.setMin(Eigen::Vector4f (b[0],b[1],b[2], 1));
	cropBoxFilter.setMax(Eigen::Vector4f (b[3],b[4],b[5], 1));
	cropBoxFilter.setNegative(true);
	cropBoxFilter.filter(*outC);
	return outC;
}
void writeCSV(pcl::PointCloud<pcl::PointXYZI>::Ptr p_, std::string fname)
{
	std::ofstream fout(fname);
	for(int i=0;i<p_->points.size();i++)
	{
		fout<<p_->points[i].x<<", ";
		fout<<p_->points[i].y<<", ";
		fout<<p_->points[i].z<<", ";
		fout<<p_->points[i].intensity<<std::endl;
	}
	fout.close();
}
int main(int argc, char** argv)
{
	string transforms,out_folder;
	string ndt_folder;
	char IFS=',';
	float parameter;
	bool Inv=false;
	po::options_description desc("Allowed options");
	desc.add_options()
	("help,h", "produce help message")
	("skip", "skip point cloud first line")
	("nc,n", "Do not concatenate transforms")
	 ("inv,i",  "Inverse apply transforms")
	 ("transforms,t", po::value<std::string >(&transforms), "File with initial transforms")
	 ("pointclouds,p", po::value<std::vector<string> >()->multitoken(), "Point cloud files")
	 ("b", po::value<std::vector<float> >()->multitoken(), "Bounding box--Atention! not working yet!")
	 ("parameter", po::value<float>(&parameter), "Bounding box--Atention! not working yet!")
	 ("ifs,f", po::value<char>(&IFS), "Pointcloud IFS")
	 ("ndt_folder", po::value<std::string >(&ndt_folder)->default_value("/tmp/semantic/"), "Folder to load/save computed NDTs")
	 ("pcl_out_folder", po::value<std::string >(&out_folder)->default_value("/tmp/"), "Folder to save registered clouds")
	 ("sem,s", po::value<std::vector<string> >()->multitoken(), "First semantic input files");

	po::parsed_options parsed_options = po::command_line_parser(argc, argv)
      .options(desc).style(po::command_line_style::unix_style ).run();

   std::vector<std::vector<std::string>> sem_files;
   for (const po::option& o : parsed_options.options) {
      if (o.string_key == "sem")
         sem_files.push_back(o.value);
   }

    po::variables_map vm;
//    po::store(po::parse_command_line(argc, argv, desc, po::command_line_style::unix_style ^ po::command_line_style::allow_short), vm);
	po::store(parsed_options, vm);
    po::notify(vm);
	vector<string> pointcloud_files;
	vector<float> box;

	bool trans=false;
	ifstream in_trans;
	if(vm.count("transforms")) { in_trans.open(transforms, ifstream::in); trans=true; }

	bool skip=false;
	if(vm.count("skip")) skip=true;
	if(vm.count("inv"))Inv=true;
	bool conc=true;
	if(vm.count("nc")) conc=false;
	if(vm.count("help"))
	{
		cout<<desc;
		return 0;
	}
	if (!vm["pointclouds"].empty() && (pointcloud_files= vm["pointclouds"].as<vector<string> >()).size() >= 2) {
		///cout<<"success pointcloud read";
	}else {cout<<"pointclouds read failure";};
	bool h_box=false;
	if(vm.count("b")) {h_box=true;box=vm["b"].as<vector<float> >();if(box.size()!=6){cout<<"Wrong box size! must be 6."<<endl;return -1;}}
	int num_files=pointcloud_files.size();
//	cout<<"Number of semantic labels, only if file per label: "<<sem_files.size()<<std::endl;
	for(int i=0;i<sem_files.size();i++)
	{
		num_files=min(num_files,(int ) sem_files.at(i).size());
//		cout<<"Number of input files for semantic label: "<<num_files<<endl;
}


	Eigen::Affine3d T;
	Eigen::Affine3d Tt,IdentityM;
	T.setIdentity();
	Tt.setIdentity();
	IdentityM.setIdentity();
	//NDTMatch_SE matcher ({0.5,0.1,0.05},{0,1,0,1,2},{25,25,10},{3},{-1},0.60,25);
	//NDTMatch_SE matcher ({200,60,200,70,50,5},{0,1,2,3,4,5},{200,200,200},{'*'},{0},0.01,50);// :-D
	//lslgeneric::NDTFuserHMT_SE matcher (the_initial_pose,{the_resolutions},{the_order_with which_the_resolutions_are_used},{the_size_of_the_map},{the_tail_segments},{ignore_values},reject_percentage,number_of_iterations);
#ifndef BULK_ADJUSTMENT
	NDTMatch_SE matcher ({1,2,0.5},{0,1,0,2},{200,200,200},{'=','=','=','*'},{1,1,1,1},0.01,50);// :-D
//	lslgeneric::NDTFuserHMT_SE matcher (T,{1,2,0.5},{1},{200,200,200},{'*'},{0},0.01,50);
	for(int i=0;i<num_files;i++)
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1=getCloud<pcl::PointXYZ>(pointcloud_files[i],IFS,skip);
		if(trans){
			Eigen::Affine3d T=readTransform(in_trans);
			pcl::transformPointCloud(*cloud1,*cloud1,T);
		}
		//Tt=matcher.match(IdentityM,cloud1,{std::vector<double>()});
		Tt=matcher.match(cloud1,{getMeasure(sem_files[0][i]),getMeasure(sem_files[1][i]),getMeasure(sem_files[2][i]),getMeasure(sem_files[2][i])});
		//matcher.updateMap();
		if(i!=0)
		{
			if(Inv) T=Tt*T;
			else T=T*Tt;
			for(int i=0;i<4;i++)
				for(int j=0;j<4;j++)
					cout<<(conc?T(i,j):Tt(i,j))<<", ";
			cout<<endl;
			pcl::PointCloud<pcl::PointXYZI>::Ptr cloud2=getCloud<pcl::PointXYZI>(pointcloud_files[i],IFS,skip);
			pcl::transformPointCloud(*cloud2,*cloud2,conc?T:Tt);
			writeCSV(cloud2,out_folder+pointcloud_files[i].substr(pointcloud_files[i].find_last_of("/\\")+1));
			Eigen::MatrixXd tmM(6,6);
			tmM=matcher.getPoseCovariance(Tt);
			cerr<<Eigen::Map<Eigen::RowVectorXd>(tmM.data(),tmM.size())<<std::endl;
		}
	}
#else
	struct stat st = {0};

	if (stat(ndt_folder.c_str(), &st) == -1) {
		mkdir(ndt_folder.c_str(), 0700);
	}
	//NDTMatch_SE matcher ({200,60,200,70,50,5},{0,1,2,3,4,5},{200,200,200},{'*'},{0},0.01,50);// :-D
	//NDTMatch_SE matcher ({100,20,5,2,1,0.5},{0,1,2,3,4,3,4,5},{200,200,200},{'*'},{0},0.01,5,true);// :-D
	//NDTMatch_SE matcher ({64,9,2,1},{0,1,3,2,2},{200,200,200},{'=','=','=','=','=','=','=','='},{1,2,3,4,5,6,7,8},0.01,5);// :-D
	//NDTMatch_SE matcher ({170,11,0.5},{0,1,2},{200,200,200},{'=','=','=','=','=','=','=','='},{1,2,300,4,5,6,7,8},0.01,5);// :-D the only good is that
	NDTMatch_SE matcher ({0.3},{0},{100,100,100},{'=','=','=','*'},{1,1,1,1},0.01,50);// :-D
	matcher.setNeighbours(2);
	matcher.IFS=IFS;
	matcher.skip=skip;
	matcher.precomputed_ndt_folder=ndt_folder;
	matcher.useSaved=false;
#ifdef G2O_OPT
	HyperGraph::VertexSet vertices;
	HyperGraph::EdgeSet edges;
	SparseOptimizer optimizer;
	OptimizationAlgorithmProperty solverProperty;
	OptimizationAlgorithmFactory* solverFactory = OptimizationAlgorithmFactory::instance();
	optimizer.setAlgorithm(solverFactory->construct("lm_fix6_3", solverProperty));
	if (! optimizer.solver()) {
		cerr << "Error allocating solver. Allocating \"" << "strSolver" << "\" failed!" << endl;
		return 0;
	}
	for(int i=0;i<num_files;i++)
	{
		VertexSE3* vert=new VertexSE3();
		vert->setId(i);
		if(i==0)
		{
			vert->setFixed(true);
			Eigen::Isometry3d Tr;
			Tr.setIdentity();
			vert->setEstimate(Tr);
		}
		optimizer.addVertex(vert);
	}
	for(int nIter=0;nIter<3;nIter++)
	{
		int i=0;
#else
	for(int k=0;k<4;k++)
		for(int l=0;l<4;l++)
			cout<<T(k,l)<<", ";
	cout<<endl;
	std::vector<Eigen::Affine3d > allTrans;
	std::vector<Eigen::Matrix<double,6,6> > allInf;
	Eigen::Matrix<double,6,6> InTot;
	Eigen::Matrix<double,6,1> VTot;
	InTot.setIdentity();
	allInf.push_back(InTot);
	T.setIdentity();
	allTrans.push_back(T);
	int i=1;
#endif
	for(i;i<num_files;i++)
	{
		int s_=i-5;
		if(s_<0)s_=0;
#ifndef G2O_OPT
			InTot.setZero();
#endif
		for(int j=s_;j<i;j+=1)
		{
//			Tt=matcher.match(pointcloud_files[i],pointcloud_files[j],{std::vector<double>()},{std::vector<double>()});
			Tt.setIdentity();
#ifdef G2O_OPT
			Tt=static_cast<VertexSE3*>(optimizer.vertex(i))->estimate()*static_cast<VertexSE3*>(optimizer.vertex(j))->estimate().inverse();
#endif
			Tt=matcher.match(Tt,pointcloud_files[j],pointcloud_files[i],{getMeasure(sem_files[0][j]),getMeasure(sem_files[1][j]),getMeasure(sem_files[2][j]),getMeasure(sem_files[2][j])},{getMeasure(sem_files[0][i]),getMeasure(sem_files[1][i]),getMeasure(sem_files[2][i]),getMeasure(sem_files[2][i])});
			std::cerr<<"["<<i<<","<<j<<"]\t"<<Tt.translation().transpose()<<std::endl;
			Eigen::Matrix<double,6,6> InM;
			InM.setIdentity();
			InM.block<6,6>(0,0)=matcher.getPoseCovariance(Tt).block<6,6>(0,0);
#ifndef G2O_OPT
			Tt=allTrans[j]*Tt;
			InM=(InM.inverse()+allInf[j].inverse()).inverse();
			Eigen::MatrixXd ea(6,1);
			ea.block<3,1>(0,0)=Tt.translation();
			ea.block<3,1>(3,0)=Tt.rotation().eulerAngles(2,0,2);
			VTot=VTot+InM*ea;
			InTot=InTot+InM;
		}
		VTot=InTot.inverse()*VTot;
			cerr<<VTot<<endl;
		Tt.translation()=VTot.block<3,1>(0,0);
		Tt.matrix().block<3,3>(0,0)=(Eigen::AngleAxisd(VTot[3], Eigen::Vector3d::UnitZ())
			 * Eigen::AngleAxisd(VTot[4], Eigen::Vector3d::UnitX())
			 * Eigen::AngleAxisd(VTot[5], Eigen::Vector3d::UnitZ())).matrix(); 
		allTrans.push_back(Tt);
		allInf.push_back(InTot);
		pcl::PointCloud<pcl::PointXYZI>::Ptr cloud2=getCloud<pcl::PointXYZI>(pointcloud_files[i],IFS,skip);
		pcl::transformPointCloud(*cloud2,*cloud2,Tt);
		writeCSV(cloud2,out_folder+pointcloud_files[i].substr(pointcloud_files[i].find_last_of("/\\")+1));
		for(int k=0;k<4;k++)
			for(int l=0;l<4;l++)
				cout<<Tt(k,l)<<", ";
		cout<<endl;
	}
#else
			if(Tt.translation().norm()>(i-j)*1)
				continue;
			EdgeSE3 *edge=NULL;
			for(std::set<HyperGraph::Edge*>::iterator it=optimizer.vertex(i)->edges().begin();
					it !=optimizer.vertex(i)->edges().end();++it)
			{
				EdgeSE3 *e=static_cast<EdgeSE3*>(*it);
				if(e->vertices()[0]==optimizer.vertex(j))edge=e;
			}
			bool incremental=true;
			if(edge==NULL||incremental)
			{
				edge=new EdgeSE3();
				edge->vertices()[0]=optimizer.vertex(j);
				edge->vertices()[1]=optimizer.vertex(i);
			}
			Eigen::Isometry3d b;
			b.translation() = Tt.translation();
			b.linear() = Tt.rotation();
			edge->setMeasurement(b);
			edge->setInformation(InM);
			cerr<<InM<<endl;
//			edge->setInformation(matcher.getPoseCovariance(Tt));
			optimizer.addEdge(edge);
		}
	}
	optimizer.initializeOptimization();
	optimizer.optimize(2000);
	}
	for(int k=0;k<num_files;k++)
	{
		Eigen::Affine3d Tp=static_cast<VertexSE3*>(optimizer.vertex(k))->estimate();
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				cout<<Tp(i,j)<<", ";
		cout<<endl;
			pcl::PointCloud<pcl::PointXYZI>::Ptr cloud2=getCloud<pcl::PointXYZI>(pointcloud_files[k],IFS,skip);
			pcl::transformPointCloud(*cloud2,*cloud2,Tp);
			writeCSV(cloud2,out_folder+pointcloud_files[k].substr(pointcloud_files[k].find_last_of("/\\")+1));
	}
	optimizer.save("tutorial_before.g2o");
	optimizer.clear();
	// destroy all the singletons
	//Factory::destroy();
	OptimizationAlgorithmFactory::destroy();
	HyperGraphActionLibrary::destroy();
#endif

#endif

	return 0;
}

