#ifndef SE_NDT
#define SE_NDT
#include <vector>
#include <pcl/point_cloud.h>
#include <ndt_map/ndt_map.h>
#include <ndt_map/lazy_grid.h>
#include <ndt_map/ndt_cell.h>
using namespace std;
typedef Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> ET;
size_t count_tails(vector<int>& distribution_tails);
size_t* sort_pointcloud(vector<double> &in,float disregard);
inline size_t index_selector(size_t **I,int p,int num,std::vector<int> Tails,size_t number_points);
inline bool checkInLimits(size_t **in,int p,int num,int cu,int cl);
vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> getSegments(pcl::PointCloud<pcl::PointXYZ>::Ptr laserCloudIn,initializer_list<vector<double> >& attributes_,initializer_list<int > distribution_tails_,initializer_list<float> disregard_, float rejectPerc);
lslgeneric::NDTMap **initMap(initializer_list<int> distribution_tails_,initializer_list<float> resolutions_, initializer_list<float>size_);
void loadMap(lslgeneric::NDTMap **map,std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> input_clouds,size_t number_tails,float sensor_range=100);
#endif/*SE_NDT*/
