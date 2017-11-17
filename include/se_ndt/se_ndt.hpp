#ifndef SE_NDT
#define SE_NDT
#include <vector>
#include <pcl/point_cloud.h>
#include <pcl/common/io.h>
#include <ndt_map/ndt_map.h>
#include <ndt_map/lazy_grid.h>
#include <ndt_map/ndt_cell.h>
#include "se_ndt/ndt_matcher_d2d_se.h"
using namespace std;
Eigen::Matrix<double,6,6> getHes(Eigen::Matrix<double,6,6> Hessian,Eigen::Matrix<double,6,1> score_gradient);
class NDTMatch_SE{
 public:
		vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> getSegments(pcl::PointCloud<pcl::PointXYZ>::Ptr laserCloudIn,initializer_list<vector<double> >& attributes_,initializer_list<int > distribution_tails_,initializer_list<float> disregard_, float rejectPerc);
		unsigned int NumInputs;
		lslgeneric::NDTMap ***map;		 ///< da map
		lslgeneric::NDTMap ***mapLocal;		 ///< da map
		vector<float> resolutions;
		initializer_list<float> ignore,size;
		initializer_list<int> resolutions_order,tails;
		float removeProbability;
		NDTMatch_SE(initializer_list<float> b,initializer_list<int> c,initializer_list<float> d,initializer_list<int> e,initializer_list<float> ig,float removeP,int max_iter,bool load_from_file=false);
		NDTMatch_SE(){};
		~NDTMatch_SE()
		{
			for(unsigned int i=0;i<resolutions.size();i++)
			{
				for(auto j=0;j<NumInputs;j++)
				{
					delete map[i][j];
					delete mapLocal[i][j];
				}
				delete[] map[i];
				delete[] mapLocal[i];
			}
			delete[] map;
		}
		Eigen::Affine3d match(Eigen::Affine3d Tinit, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2, initializer_list<vector<double> > attributes1, initializer_list<vector<double> > attributes2);
		Eigen::Affine3d match(Eigen::Affine3d Tinit, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, initializer_list<vector<double> > attributes);
		Eigen::Affine3d match(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, initializer_list<vector<double> > attributes);
		Eigen::Affine3d match(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1,pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2, initializer_list<vector<double> > attributes1, initializer_list<vector<double> > attributes2);
		Eigen::Affine3d match(std::string cloudF1, std::string cloudF2,initializer_list<vector<double> > attributes1,initializer_list<vector<double> > attributes2);
		lslgeneric::NDTMatcherD2D_SE matcher;
		Eigen::Matrix<double,6,6> getPoseCovariance(Eigen::Affine3d T);
		void setNeighbours(short int i){matcher.n_neighbours=i;};
		char IFS=',';
		bool skip=false;
		float sensor_range=100;
    private:
		bool firstRun;
		bool load_from_file;
		std::vector<int> semantic_labels;

		Eigen::Vector3d localMapSize;
};
typedef Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> ET;
size_t count_tails(vector<int>& distribution_tails);
size_t* sort_pointcloud(vector<double> &in,float disregard);
Eigen::Affine3d readTransform(istream &infile);
inline size_t index_selector(size_t **I,int p,int num,std::vector<int> Tails,size_t number_points);
inline bool checkInLimits(size_t **in,int p,int num,int cu,int cl);
lslgeneric::NDTMap **initMap(int number_tails,initializer_list<float> resolutions_, initializer_list<float>size_,bool load_from_file=false);
void loadMap(lslgeneric::NDTMap **map,std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> input_clouds,float sensor_range=100);
template <typename T> typename pcl::PointCloud<T>::Ptr getCloud(string filename,char IFS, bool skip);

#endif/*SE_NDT*/
