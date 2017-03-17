#include <vector>
#include <pcl/point_cloud.h>
#include <ndt_map/ndt_map.h>
#include <ndt_map/lazy_grid.h>
#include <ndt_map/ndt_cell.h>
using namespace std;
size_t* sort_pointcloud(vector<double> &in,float disregard)
{
	vector<size_t> idx(in.size());
	size_t *idx2 = new size_t[in.size()];
	iota(idx.begin(),idx.end(),0);
	sort(idx.begin(), idx.end(),[&in](size_t i1,size_t i2){return in[i1]<in[i2];});
	std::vector<size_t>::iterator up= std::upper_bound (idx.begin(), idx.end(), disregard,[&in](double v,size_t i1){return in[i1]>v;}); //                   
	size_t mid = (idx.end()-up)/2;
	  std::swap_ranges(up, up+mid, idx.begin());
	  int j=0;
	  for(auto i:idx)
		  idx2[i]=j++;
	return idx2;
}
inline size_t index_selector(size_t **I,int p,int num,std::vector<size_t> Tails,size_t number_points)
{
	size_t minI=0,minA=I[0][p];
	size_t maxI=0,maxA=I[0][p];
	for(auto i=0;i<num;i++)
	{
		if(*(Tails.begin()+i)&1)if(I[i][p]<=minA) { minA=I[i][p];minI=i; }
		if(*(Tails.begin()+i)&2)if(I[i][p]>=maxA) { maxA=I[i][p];maxI=i; }
	}
	return number_points-maxA<minA?maxI+num:minI;
}
inline bool checkInLimits(size_t **in,int p,int num,int cu,int cl)
{
	for(int i=0;i<num;i++)
		if(in[i][p]>cu||in[i][p]<cl)
			return true;
	return false;
}
vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> getSegments(pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudIn,initializer_list<vector<double> >& attributes_,initializer_list<size_t > distribution_tails_,initializer_list<float> disregard_, float rejectPerc)
{
	vector<vector<double> > attributes(attributes_);
	vector<float> disregard(disregard_);
	vector<size_t> distribution_tails(distribution_tails_);
	int cloudSize = laserCloudIn->points.size();
	int num_attributes=attributes.size();
	size_t *sorted_indexes[num_attributes];
	for(auto i=0;i< num_attributes;i++)
		sorted_indexes[i]=sort_pointcloud(attributes.at(i),disregard.at(i));

	int cut_off_l=(1-rejectPerc)/2*cloudSize;
	int cut_off_u=(1+rejectPerc)/2*cloudSize;
	////////
	//Create look-up
	size_t look_up[2*num_attributes];
	size_t num_tails=0;
	for(int i=0;i<num_attributes;i++)
	{
		if(distribution_tails[i]&1)
		{
			look_up[2*i]= num_tails;
			num_tails++;
		}
		if(distribution_tails[i]&2)
		{
			look_up[2*i+1]=num_tails;
			num_tails++;
		}
	}
	///////
	vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud;
	for(int i=0;i<num_tails;i++)
		laserCloud.emplace_back(new pcl::PointCloud<pcl::PointXYZ>);
	for(int i=0;i<cloudSize;i++)
	{
		if(!checkInLimits(sorted_indexes,i,num_attributes,cut_off_u,cut_off_l)) continue;
		pcl::PointXYZ point;
		point.x=laserCloudIn->points[i].x;
		point.y=laserCloudIn->points[i].y;
		point.z=laserCloudIn->points[i].z;
		//point.intensity=laserCloudIn->points[i].intensity;
		int index=look_up[index_selector(sorted_indexes, i,num_attributes,distribution_tails,cloudSize)];
		laserCloud[index]->points.push_back(point);
	}
	for(auto i=0;i<num_attributes;i++)
		delete sorted_indexes[i];


	return laserCloud;
}
inline size_t count_tails(vector<size_t>& distribution_tails)
{
	size_t number_tails1=count(distribution_tails.begin(),distribution_tails.end(),1);
	size_t number_tails2=count(distribution_tails.begin(),distribution_tails.end(),2);
	size_t number_tails3=count(distribution_tails.begin(),distribution_tails.end(),3);
	return number_tails1+number_tails2+2*number_tails3;
}

lslgeneric::NDTMap **initMap(initializer_list<size_t > distribution_tails_,initializer_list<float> resolutions_, initializer_list<float>size_)
{
	vector<size_t> distribution_tails(distribution_tails_);
	vector<float> size(size_),resolutions(resolutions_);
	size_t number_tails=count_tails(distribution_tails);
	if(number_tails!=resolutions.size()&&resolutions.size()!=1)
		cerr<<"Number of resolutions different than number of segments. Taking resolution index modulus."<<endl;
	if(size.size()!=3)
	{
		float max_size=*max_element(size.begin(),size.end());
		size[0]=max_size;
		size[1]=max_size;
		size[2]=max_size;
		cerr<<"Wrong size parameter. Using size x=y=z="<<max_size<<endl;
	}
	lslgeneric::NDTMap **map = new lslgeneric::NDTMap * [number_tails];
	for(size_t i=0;i<number_tails;i++)
	{
		lslgeneric::LazyGrid *grid = new lslgeneric::LazyGrid(resolutions[i%resolutions.size()]);
		lslgeneric::NDTMap *mapTMP=new lslgeneric::NDTMap(grid);
		map[i]=mapTMP;
		map[i]->guessSize(0,0,0,size[0],size[1],size[2]);
		map[i]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE);
	}
	return map;
}
void updateMap(lslgeneric::NDTMap **map,std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> input_clouds,size_t number_tails,float sensor_range=100)
{
	for(size_t i=0;i<number_tails;i++)
	{
		map[i]->loadPointCloud(*input_clouds[i],sensor_range);
		map[i]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE);
	}
}

