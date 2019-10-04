#ifndef SE_NDT
#define SE_NDT
#include <vector>
#include <csignal>
#include <pcl/point_cloud.h>
#include <pcl/common/io.h>
#include <pcl/common/geometry.h>
#include <ndt_map/ndt_map.h>
#include <ndt_map/lazy_grid.h>
#include <ndt_map/ndt_cell.h>
#include <ndt_map/pointcloud_utils.h>
#include "se_ndt/ndt_matcher_d2d_se.h"
#include "se_ndt/ndt_histogram.h"
#include <profiler.hpp>
#ifdef GL_VISUALIZE
    #include "ndt_visualisation/ndt_viz.h"
#endif
#include <Eigen/StdVector>
#include <thread>
//#define GROUND_TRUTH_POSE 
using namespace std;
Eigen::Matrix<double,6,6> getHes(Eigen::Matrix<double,6,6> Hessian,Eigen::Matrix<double,6,1> score_gradient);
class NDTMatch_SE{
    public:
        void loadMap(int start_index, int stop_index, pcl::PointXYZL &pose);
        std::vector<std::vector<perception_oru::NDTCell*> > toRVIZ;
        perception_oru::NDTMap** loadSavedMap(int index);
        bool useSaved=false;
        unsigned int NumInputs;
		perception_oru::NDTMap ***map;
		perception_oru::NDTMap ***mapLocal_prev;
		perception_oru::NDTMap ***mapLocal;
		vector<float> resolutions,size;
		vector<int> resolutions_order;
        Eigen::Affine3d T, Td;
        std::vector<perception_oru::NDTHistogram> hists;
        std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d> > Ti;
//        std::vector<Eigen::Matrix<double, 7,7>, Eigen::aligned_allocator<Eigen::Matrix<double,7,7> > > covS;
        pcl::KdTreeFLANN<pcl::PointXYZL> pose_kdtree;
        pcl::PointCloud<pcl::PointXYZL>::Ptr poses;
        unsigned int num_clouds=0;
        float max_size;
        //std::vector<Eigen::Matrix<double, 1024,1> > gfeat;

        Eigen::Affine3d matchFaster(Eigen::Affine3d Tinit, pcl::PointCloud<pcl::PointXYZI>::Ptr cloud);
        Eigen::Affine3d matchFaster(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud);
        vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> getSegmentsFast(pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudIn);
        Eigen::Affine3d matchFaster_OM(Eigen::Affine3d Tinit, pcl::PointCloud<pcl::PointXYZI>::Ptr cloud);
        Eigen::Affine3d mapUpdate(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud, bool constant_velocity=true);

		NDTMatch_SE(initializer_list<float> b,initializer_list<int> c,initializer_list<float> d,int nIn,int max_iter);
        void loadMap(perception_oru::NDTMap **map,std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> input_clouds,float sensor_range=100);
		NDTMatch_SE(int nIn): NDTMatch_SE({4,0.8},{0,1},{80,80,80},nIn,50){ };
		NDTMatch_SE(): NDTMatch_SE({4,0.8},{0,1},{80,80,80},8,50){ };
		~NDTMatch_SE();
		perception_oru::NDTMatcherD2D_SE matcher;
		Eigen::Matrix<double,7,7> getPoseInformation(Eigen::Affine3d T, perception_oru::NDTMap ***m1, perception_oru::NDTMap*** m2, bool inverse);
//        Eigen::Matrix<double,7,7> CovSum;
		void setNeighbours(short int i){matcher.n_neighbours=i;};
		float sensor_range=100;
        Eigen::Affine3d simple_match(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud1, pcl::PointCloud<pcl::PointXYZI>::Ptr cloud2);
        void print_vals()
        {
            std::ofstream fout("results.txt");
            int n_poses=hists.size();
            for(int i=0;i<n_poses;i++)
                for(int j=i+1;j<n_poses;j++)
                {
                    fout<<j-i<<" ";
                    fout<<pcl::geometry::distance(poses->at(j), poses->at(i))<<" ";
                    fout<<hists[i].getSimilarity(hists[j])<<" ";
                    //fout<<(gfeat[poses->at(j).label]-gfeat[poses->at(i).label]).norm()/1024<<" ";
                    fout<<std::endl;
                }
            fout.close();
        };
    private:
        int last_loop_close_id=0;
        double last_loop_close_sim=0;
        bool matchToSaved(Eigen::Affine3d &Td_, pcl::PointXYZL &pose_end, pcl::PointXYZL &pose_current, int start_index, int iP, Eigen::Matrix<double,7,7> &Ccl);
        int find_start(pcl::PointCloud<pcl::PointXYZL>::iterator pi, float max_size);
		bool firstRun;
#ifdef GL_VISUALIZE
    public:
        NDTViz *viewer = NULL;
        void visualize();
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    private:
        void visualize_thread();
#endif
};
typedef Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> ET;
Eigen::Matrix<double,7,6> getJacobian(Eigen::VectorXd v);
#endif/*SE_NDT*/
