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
#include <pcl/io/pcd_io.h>
#include <pcl/console/parse.h>
#include <se_ndt/se_ndt.hpp>
#include <se_ndt/ndt_matcher_d2d_se.h>
#include <omp.h>

#include <pcl/filters/extract_indices.h>


#include <pcl/io/pcd_io.h>
#include <pcl/features/normal_3d.h>
#include <pcl/features/rsd.h>

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
std::vector<double> getRSD(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud, float a,float b, float c,float d)
{
	int cloudSize =cloud->points.size();
	std::vector<double> rsd_v(cloudSize,-1);
	// Object for storing the normals.
	pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
	// Object for storing the RSD descriptors for each point.
	pcl::PointCloud<pcl::PrincipalRadiiRSD>::Ptr descriptors(new pcl::PointCloud<pcl::PrincipalRadiiRSD>());

	pcl::NormalEstimation<pcl::PointXYZI, pcl::Normal> normalEstimation;
	normalEstimation.setInputCloud(cloud);
	normalEstimation.setKSearch(a);
	pcl::search::KdTree<pcl::PointXYZI>::Ptr kdtree(new pcl::search::KdTree<pcl::PointXYZI>);
	normalEstimation.setSearchMethod(kdtree);
	normalEstimation.compute(*normals);

	// RSD estimation object.
	pcl::RSDEstimation<pcl::PointXYZI, pcl::Normal, pcl::PrincipalRadiiRSD> rsd;
	rsd.setInputCloud(cloud);
	rsd.setInputNormals(normals);
	rsd.setSearchMethod(kdtree);
	// Search radius, to look for neighbors. Note: the value given here has to be
	// larger than the radius used to estimate the normals.
	rsd.setRadiusSearch(b);
	// Plane radius. Any radius larger than this is considered infinite (a plane).
	rsd.setPlaneRadius(c);
	// Do we want to save the full distance-angle histograms?
	rsd.setSaveHistograms(false);

	rsd.compute(*descriptors);
	if(d==0)
		for(int i=0;i<cloudSize;i++)
			rsd_v[i]=descriptors->points[i].r_min;//Recomended 6 0.1 4 
	if(d==1)
		for(int i=0;i<cloudSize;i++)
			rsd_v[i]=descriptors->points[i].r_max;//Recomended 6 0.1 50
	if(d==2)
		for(int i=0;i<cloudSize;i++)
			rsd_v[i]=descriptors->points[i].r_max/descriptors->points[i].r_min;
	if(d==3)
		for(int i=0;i<cloudSize;i++)
			rsd_v[i]=descriptors->points[i].r_min/descriptors->points[i].r_max;

	return rsd_v;
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
std::vector<double> getCornerness(pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudIn,float R)
{
	int cloudSize = laserCloudIn->points.size();
	pcl::KdTreeFLANN<pcl::PointXYZI> kdtree;
	kdtree.setInputCloud (laserCloudIn);
	std::vector<double> cornerness(cloudSize,-1);
	int rem=0;
#define n_threads 12
    #pragma omp parallel num_threads(n_threads)
	{
        #pragma omp for
		for(int i=0;i<cloudSize;i++)
		{
			std::vector<int> pointIdxKNNSearch(10);
			std::vector<float> pointDistance(10);
			////////////
//			std::vector<int> pointIdxKNNSearchOO(2);
//			std::vector<float> pointDistanceOO(2);
			///////////
			std::vector<float> diffXYZ (3,0);
			if(kdtree.radiusSearch(laserCloudIn->points[i],R,pointIdxKNNSearch,pointDistance)>0)
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



