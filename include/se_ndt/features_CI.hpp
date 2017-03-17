#include <cmath>
#include <iostream>
#include <vector>
#include <chrono>
#include <ctime>

#include <opencv/cv.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/filters/voxel_grid.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/common/transforms.h>
#include <pcl/filters/approximate_voxel_grid.h>
#include <pcl/filters/crop_box.h>
#include <pcl/range_image/range_image.h>
#include <pcl/io/pcd_io.h>
#include <pcl/features/range_image_border_extractor.h>
#include <pcl/keypoints/narf_keypoint.h>
#include <pcl/filters/random_sample.h>
#include <pcl/console/parse.h>
#include <ndt_registration/ndt_matcher_d2d.h>
#include <se_ndt/se_ndt.hpp>
#include <se_ndt/ndt_matcher_d2d_se.h>
#include <omp.h>

#include <pcl/filters/extract_indices.h>

#include <boost/thread/thread.hpp>

pcl::PointCloud<pcl::PointXYZ> getNARF(pcl::PointCloud<pcl::PointXYZI>::Ptr point_cloud_ptr)
{
	  pcl::PointCloud<pcl::PointXYZ> point_cloud ;
	  pcl::copyPointCloud(*point_cloud_ptr,point_cloud);
	float support_size = 0.01f;
	float noise_level = 0.0;
	float min_range = 0.0f;
	int border_size = 0;
	float angular_resolution = 0.05f;
	  angular_resolution = pcl::deg2rad (angular_resolution);
	    Eigen::Affine3f scene_sensor_pose (Eigen::Affine3f::Identity ());
		pcl::RangeImage::CoordinateFrame coordinate_frame = pcl::RangeImage::CAMERA_FRAME;
		  pcl::PointCloud<pcl::PointWithViewpoint> far_ranges;
	boost::shared_ptr<pcl::RangeImage> range_image_ptr (new pcl::RangeImage);
	pcl::RangeImage& range_image = *range_image_ptr;   
	range_image.createFromPointCloud (point_cloud, angular_resolution, pcl::deg2rad (360.0f), pcl::deg2rad (180.0f),
								   scene_sensor_pose, coordinate_frame, noise_level, min_range, border_size);
	range_image.integrateFarRanges (far_ranges);
	range_image.setUnseenToMaxRange ();

  pcl::RangeImageBorderExtractor range_image_border_extractor;
  pcl::NarfKeypoint narf_keypoint_detector (&range_image_border_extractor);
  narf_keypoint_detector.setRangeImage (&range_image);
  narf_keypoint_detector.getParameters ().support_size = support_size;
  
  pcl::PointCloud<int> keypoint_indices;
  narf_keypoint_detector.compute (keypoint_indices);
  std::cout << "Found "<<keypoint_indices.points.size ()<<" key points.\n";

  pcl::PointCloud<pcl::PointXYZ>::Ptr keypoints_ptr (new pcl::PointCloud<pcl::PointXYZ>);
  pcl::PointCloud<pcl::PointXYZ>& keypoints = *keypoints_ptr;
  keypoints.points.resize (keypoint_indices.points.size ());
  for (size_t i=0; i<keypoint_indices.points.size (); ++i)
    keypoints.points[i].getVector3fMap () = range_image.points[keypoint_indices.points[i]].getVector3fMap ();
 return keypoints;
}

using namespace std;

int occluded(pcl::PointXYZI a, pcl::PointXYZI b, float d)
{
	float d1 = sqrt(a.x*a.x+a.y*a.y+a.z*a.z);
	float d2 = sqrt(b.x*b.x+b.y*b.y+b.z*b.z);
	if(d1>d2)
	{
		float dx=b.x-a.x*d2/d1;
		float dy=b.y-a.y*d2/d1;
		float dz=b.z-a.z*d2/d1;
		if(sqrt(dx*dx+dy*dy+dz*dz)/d2<d)
			return 1;
	}
	return 0;
}
std::vector<double> getCornerness2(pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudIn,int K)
{
	int cloudSize = laserCloudIn->points.size();
	pcl::KdTreeFLANN<pcl::PointXYZI> kdtree;
	kdtree.setInputCloud (laserCloudIn);
	std::vector<double> cornerness(cloudSize,-1);
	int rem=0;
#define n_threads 8
    #pragma omp parallel num_threads(n_threads)
	{
        #pragma omp for
		for(int i=0;i<cloudSize;i++)
		{
			std::vector<int> pointIdxKNNSearch(K);
			std::vector<float> pointDistance(K);
			////////////
//			std::vector<int> pointIdxKNNSearchOO(2);
//			std::vector<float> pointDistanceOO(2);
			///////////
			std::vector<float> diffXYZ (3,0);
			if(kdtree.nearestKSearch(laserCloudIn->points[i],K,pointIdxKNNSearch,pointDistance)>0)
			{
				int minCornIndex=i;
				int maxCornIndex=i;
				bool cf=1;
				int numKNN=pointIdxKNNSearch.size();
				for(int j=0;j<numKNN;++j)
				{
					diffXYZ[0]+=laserCloudIn->points[pointIdxKNNSearch[j]].x-laserCloudIn->points[i].x;
					diffXYZ[1]+=laserCloudIn->points[pointIdxKNNSearch[j]].y-laserCloudIn->points[i].y;
					diffXYZ[2]+=laserCloudIn->points[pointIdxKNNSearch[j]].z-laserCloudIn->points[i].z;
					float disOcc=0.2;
					if(pointDistance[j]>disOcc&&
							occluded(laserCloudIn->points[i],laserCloudIn->points[pointIdxKNNSearch[j]],disOcc))
					{
						cf=0;
						break;
					}
					if(cornerness[j]<cornerness[minCornIndex]||cornerness[minCornIndex]==-1)
						minCornIndex = j;
					else if(cornerness[j]>cornerness[maxCornIndex])
						maxCornIndex = j;
					else cornerness[j]=-1;
				}
				float distanceOfPoint =sqrt( pow(laserCloudIn->points[i].x,2)+pow(laserCloudIn->points[i].y,2)+pow(laserCloudIn->points[i].z,2));
				if(cf)cornerness[i]=numKNN>0?sqrt(pow(diffXYZ[0],2)+pow(diffXYZ[1],2)+pow(diffXYZ[2],2))/distanceOfPoint/numKNN:-1;
				if(cornerness[i]!=-1)
				{
					if(cornerness[i]<cornerness[minCornIndex]&&minCornIndex!=i)
						cornerness[minCornIndex]=-1;
					else if(cornerness[i]>cornerness[maxCornIndex]&&maxCornIndex!=i)
						cornerness[maxCornIndex]=-1;
					else if (maxCornIndex!=i&&minCornIndex!=i)cornerness[i]=-1;
				}
			}
		}
	}
	return cornerness;
}
std::vector<double> getCornerness(pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudIn,int K)
{
	int cloudSize = laserCloudIn->points.size();
	pcl::KdTreeFLANN<pcl::PointXYZI> kdtree;
	kdtree.setInputCloud (laserCloudIn);
	std::vector<double> cornerness(cloudSize,-1);
	int rem=0;
#define n_threads 8
    #pragma omp parallel num_threads(n_threads)
	{
        #pragma omp for
		for(int i=0;i<cloudSize;i++)
		{
			std::vector<int> pointIdxKNNSearch(K);
			std::vector<float> pointDistance(K);
			////////////
//			std::vector<int> pointIdxKNNSearchOO(2);
//			std::vector<float> pointDistanceOO(2);
			///////////
			std::vector<float> diffXYZ (3,0);
			if(kdtree.radiusSearch(laserCloudIn->points[i],0.2,pointIdxKNNSearch,pointDistance)>0)
			{
				//////////
				/*
				kdtree.nearestKSearch(laserCloudIn->points[i],2,pointIdxKNNSearchOO,pointDistanceOO);
				Eigen::Vector3f a=laserCloudIn->points[i].getVector3fMap();
				Eigen::Vector3f b=laserCloudIn->points[pointIdxKNNSearchOO[1]].getVector3fMap();
				Eigen::Vector3f atb;
				float a1=a.dot(b/b.norm());
				float angle_AB=(a/a.norm()).dot(b/b.norm());
				atb=b/b.norm()*a1;
				atb=atb-b;
				a=a-b;
				float angle_AATB=(a/a.norm()).dot(atb/atb.norm());

				cout<<acos(angle_AATB)*180/M_PI<<endl;
				if(acos(angle_AATB)*180/M_PI<5)
				{
					cout<<"rem"<<rem<<endl;
					rem++;
					cornerness[i]=-1;
					continue;
				}

				*/
				/////////
				int minCornIndex=i;
				int maxCornIndex=i;
				bool cf=1;
				int numKNN=pointIdxKNNSearch.size();
				for(int j=0;j<numKNN;++j)
				{
					diffXYZ[0]+=laserCloudIn->points[pointIdxKNNSearch[j]].x-laserCloudIn->points[i].x;
					diffXYZ[1]+=laserCloudIn->points[pointIdxKNNSearch[j]].y-laserCloudIn->points[i].y;
					diffXYZ[2]+=laserCloudIn->points[pointIdxKNNSearch[j]].z-laserCloudIn->points[i].z;
					float disOcc=0.2;
					if(pointDistance[j]>disOcc&&
							occluded(laserCloudIn->points[i],laserCloudIn->points[pointIdxKNNSearch[j]],disOcc))
					{
						cf=0;
						break;
					}
					if(cornerness[j]<cornerness[minCornIndex]||cornerness[minCornIndex]==-1)
						minCornIndex = j;
					else if(cornerness[j]>cornerness[maxCornIndex])
						maxCornIndex = j;
					else cornerness[j]=-1;
				}
				float distanceOfPoint =sqrt( pow(laserCloudIn->points[i].x,2)+pow(laserCloudIn->points[i].y,2)+pow(laserCloudIn->points[i].z,2));
				if(cf)cornerness[i]=numKNN>0?/*sqrt(distanceOfPoint)*/sqrt(pow(diffXYZ[0],2)+pow(diffXYZ[1],2)+pow(diffXYZ[2],2))/numKNN:-1;
				if(cornerness[i]!=-1)
				{
					if(cornerness[i]<cornerness[minCornIndex]&&minCornIndex!=i)
						cornerness[minCornIndex]=-1;
					else if(cornerness[i]>cornerness[maxCornIndex]&&maxCornIndex!=i)
						cornerness[maxCornIndex]=-1;
					else if (maxCornIndex!=i&&minCornIndex!=i)cornerness[i]=-1;
				}
			}
		}
	}
	return cornerness;
}
pcl::PointCloud<pcl::PointXYZI>::Ptr Voxel_rm_NaN(pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudIn)
{
    std::vector<int> indices;
    pcl::removeNaNFromPointCloud(*laserCloudIn, *laserCloudIn, indices);
	float voxel_size=0;
	if(voxel_size!=0)
	{
			pcl::PointCloud<pcl::PointXYZ>::Ptr non_filtered_cloud (new pcl::PointCloud<pcl::PointXYZ>);
			pcl::PointCloud<pcl::PointXYZ>::Ptr input_cloud (new pcl::PointCloud<pcl::PointXYZ>);
			copyPointCloud(*laserCloudIn,*non_filtered_cloud);
			pcl::VoxelGrid<pcl::PointXYZ> approximate_voxel_filter;
			approximate_voxel_filter.setLeafSize (voxel_size, voxel_size, voxel_size);
			approximate_voxel_filter.setInputCloud (non_filtered_cloud);
			approximate_voxel_filter.filter (*input_cloud);
			copyPointCloud(*input_cloud,*laserCloudIn);
	}
	//std::vector<float> meanDist =getMeanDist(laserCloudIn);
	//std::vector<float> intensity (cloudSize);
	//std::vector<float> corRate =getCorRate(laserCloudIn,10,50,100);
	//	for(int i=0;i<cloudSize;i++)
	//		intensity[i]=laserCloudIn->points[i].intensity;
	//auto Idst=sort_pointcloud(meanDist);
	//auto Iint=sort_pointcloud(intensity);
	//auto IcR=sort_pointcloud(corRate);
	//laserCloud[1]=getNARF(laserCloudIn);
	return laserCloudIn;
}
pcl::PointCloud<pcl::PointXYZ>  getCloud(string filename)
{
	ifstream infile(filename); // for example
	string line = "";
	pcl::PointCloud<pcl::PointXYZ> laserCloud;
	getline(infile, line);
	while (getline(infile, line)){
	//cin>>line;
	//while (cin>>line){
		stringstream strstr(line);
		string word = "";
		getline(strstr,word, ',');
		pcl::PointXYZ point;
		getline(strstr,word, ',');
		point.x=stof(word);
		getline(strstr,word, ',');
		point.y=stof(word);
		getline(strstr,word, ',');
		point.z=stof(word);
		getline(strstr,word, ',');
		//point.intensity=stof(word);
		laserCloud.points.push_back(point);
	}
	return laserCloud;
}


pcl::PointCloud<pcl::PointXYZI>::Ptr cropIt(pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloud)
{
	pcl::PointCloud<pcl::PointXYZI>::Ptr outC(new pcl::PointCloud<pcl::PointXYZI>);
	pcl::CropBox<pcl::PointXYZI> cropBoxFilter(false);
	cropBoxFilter.setInputCloud(laserCloud);
	cropBoxFilter.setMin(Eigen::Vector4f (-0.35, -0.35, 0, 1));
	cropBoxFilter.setMax(Eigen::Vector4f ( 0.2,  0.3, 0.65, 1));
	cropBoxFilter.setNegative(true);
	cropBoxFilter.filter(*outC);
	return outC;
}

pcl::PointCloud<pcl::PointXYZI>::Ptr  getCloudI(string filename)
{
	ifstream infile(filename); // for example
	string line = "";
	pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloud(new pcl::PointCloud<pcl::PointXYZI>);
	getline(infile, line);
	while (getline(infile, line)){
	//cin>>line;
	//while (cin>>line){
		stringstream strstr(line);
		string word = "";
		getline(strstr,word, ',');
		pcl::PointXYZI point;
		getline(strstr,word, ',');
		point.x=stof(word);
		getline(strstr,word, ',');
		point.y=stof(word);
		getline(strstr,word, ',');
		point.z=stof(word);
		getline(strstr,word, ',');
		point.intensity=stof(word);
		(*laserCloud).points.push_back(point);
	}
    return cropIt(laserCloud);
}


std::vector<Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> > getPose(string filename)
{
	ifstream infile(filename); // for example
	string line = "";
	std::vector< Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> >poses;
	getline(infile, line);
	while (getline(infile, line)){
	//cin>>line;
	//while (cin>>line){
		stringstream strstr(line);
		string word = "";
		Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> poseI;
		getline(strstr,word, ',');
		getline(strstr,word, ',');
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
			{
				getline(strstr,word, ',');
				poseI(i,j)=stof(word);
			}
		poses.push_back(poseI);
	}
	return poses;
}
typedef Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> ET;
class timer
{
public:
	double reading_time;
	double reference_time;
	std::chrono::time_point<std::chrono::system_clock> startT, end;
    // default constructor that stores the start time
    void start()
    {
		startT= std::chrono::system_clock::now();
    }
	void stop_read(){
		end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-startT;
		reading_time= (elapsed_seconds).count();
	}
	void stop_ref(){
		end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-startT;
		reference_time= (elapsed_seconds).count();
	}
	float stop_reg(){
		end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-startT;
		double reg_time = (elapsed_seconds).count();
		return (reg_time+reading_time+reference_time);
	}

};
#define RECURS true
class matchCI{
	public:
		float size_x=20,size_y=20,size_z=20,sensor_range=50;
		//double resolution[20] = {3.0,1.0,3,1,0.1};
		double resolution[20] = {1.0,2,1,0.5};
		int num_res=4;
		int num_inputs=2;
		float removeP=0.5;
		int nK=5;
		lslgeneric::NDTMatcherD2D_SE matcher;
		lslgeneric::NDTMap ***mapReading;
		lslgeneric::NDTMap ***mapReference;
		string base_name;
		ofstream out_file;

		matchCI(string directory,bool activate,float remove_percent=0.5)
		{
			mapReading=new lslgeneric::NDTMap ** [num_res];
			mapReference=new lslgeneric::NDTMap ** [num_res];
			if(!activate)return;
			removeP=remove_percent;
			int argc=1;
			char **argv;
			string s1="reference",s2="reading";
			base_name=directory;
			out_file.open(directory+"results.csv");
			out_file<<endl;

			matcher.NumInputs=num_inputs;
			matcher.ITR_MAX =5;
			matcher.step_control=true;

			mapReading[0]=initMap({3},{0.5},{size_x,size_y,size_z});
			mapReading[1]=initMap({3},{1.0},{size_x,size_y,size_z});
			mapReading[2]=initMap({3},{2.0},{size_x,size_y,size_z});
			mapReference[0]=initMap({3},{0.5},{size_x,size_y,size_z});
			mapReference[1]=initMap({3},{1.0},{size_x,size_y,size_z});
			mapReference[2]=initMap({3},{2.0},{size_x,size_y,size_z});
		};
		~matchCI(){out_file.close();};


		string prev_ref,prev_read;
	ET estimateTransformNDT_CI(string protocol_line)
	{
		stringstream strstr(protocol_line);
		string word="";
		getline(strstr,word, ',');
		ET T;
		timer reg_timer;

		if(word.compare(prev_ref)!=0)
		{
			pcl::PointCloud<pcl::PointXYZI>::Ptr r_c = getCloudI(base_name+word);
			pcl::PointCloud<pcl::PointXYZI>::Ptr reference =Voxel_rm_NaN(r_c);
			std::vector<double> cornerness =getCornerness(reference,nK);
			auto attrs={cornerness};
			std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloudIn=getSegments(reference,attrs,{3},{-1},removeP);
			prev_ref=word;

			reg_timer.start();
			updateMap(mapReference[0],laserCloudIn,2);
			updateMap(mapReference[1],laserCloudIn,2);
			updateMap(mapReference[2],laserCloudIn,2);
			reg_timer.stop_ref();

		}
		getline(strstr,word, ',');
		word.erase(0,1);
		if(word.compare(prev_read)!=0)
		{
			pcl::PointCloud<pcl::PointXYZI>::Ptr r_c = getCloudI(base_name+word);
			pcl::PointCloud<pcl::PointXYZI>::Ptr reading = Voxel_rm_NaN(r_c);
			std::vector<double> cornerness =getCornerness(reading,nK);
			auto attrs={cornerness};
			std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> laserCloudIn=getSegments(reading,attrs,{3},{-1},removeP);
			prev_read=word;

			reg_timer.start();
			updateMap(mapReading[0],laserCloudIn,2);
			updateMap(mapReading[1],laserCloudIn,2);
			updateMap(mapReading[2],laserCloudIn,2);
			reg_timer.stop_read();

		}
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
			{
				getline(strstr,word, ',');
				T(i,j)=stof(word);
			}
		ET Tn=T;
		reg_timer.start();

		matcher.current_resolution=1;
		matcher.match(mapReference[1],mapReading[1],T,true);
		matcher.current_resolution=2;
		matcher.match(mapReference[2],mapReading[2],T,true);
		matcher.current_resolution=1;
		matcher.match(mapReference[1],mapReading[1],T,true);
		matcher.current_resolution=0.5;
		matcher.match(mapReference[0],mapReading[0],T,true);

		float final_time =	reg_timer.stop_reg();
		stringstream kstr;
		kstr<<final_time<<", ";
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				kstr<<T(i,j)<<", ";
		string buf=kstr.str();
		buf.resize(buf.size()-2);
		out_file<<buf<<endl;
		return T;
	}
};
class matchNDT{
	public:
		float size_x=10,size_y=10,size_z=5,sensor_range=50;
		double __res[4] = {1,2,1,0.5};
		int num_res=4;
		lslgeneric::NDTMatcherD2D *matcher;
		string base_name;
		ofstream out_file;

		matchNDT(string directory,bool activate)
		{
			if(!activate)return;
			std::vector<double> resolutions (__res, __res+sizeof(__res)/sizeof(double));
			matcher = new lslgeneric::NDTMatcherD2D(false,false,resolutions);
			base_name=directory;
			out_file.open(directory+"results_ndt.csv");
			matcher->ITR_MAX =100;
			matcher->step_control=true;
		};
		~matchNDT(){out_file.close();};


		string prev_ref,prev_read;
		pcl::PointCloud<pcl::PointXYZ> reference,reading;
	ET estimateTransformNDT(string protocol_line)
	{
		stringstream strstr(protocol_line);
		string word="";
		getline(strstr,word, ',');
		ET T;
		if(word.compare(prev_ref)!=0)
		{
			reference = getCloud(base_name+word);
			prev_ref=word;
		}
		getline(strstr,word, ',');
		word.erase(0,1);
		if(word.compare(prev_read)!=0)
		{
			reading = getCloud(base_name+word);
			prev_read=word;
		}
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
			{
				getline(strstr,word, ',');
				T(i,j)=stof(word);
			}
		time_t begin_time =time(NULL);
		matcher->match(reference,reading,T,true);
		double final_time = difftime(time(NULL),begin_time);
		out_file<<final_time<<",";
		for(int i=0;i<4;i++)
			for(int j=0;j<4;j++)
				out_file<<T(i,j)<<",";
		out_file<<endl;
		return T;
	}
};
