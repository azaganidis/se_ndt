#pragma once
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
#include "pose_optimizer.h"
#include <profiler.hpp>
#ifdef GL_VISUALIZE
    #include "ndt_visualisation/ndt_viz.h"
#endif
#include <Eigen/StdVector>
#include <thread>
using namespace std;

Eigen::Matrix<double,6,6> getHes(Eigen::Matrix<double,6,6> Hessian,Eigen::Matrix<double,6,1> score_gradient);
Eigen::Matrix<double,7,6> getJacobian(Eigen::VectorXd v);

class NDTMatch_SE{
    public:
		NDTMatch_SE(initializer_list<float> b,initializer_list<int> c,
                initializer_list<float> d,int nIn,int max_iter);
		NDTMatch_SE(int nIn): NDTMatch_SE({4,0.8},{0,1},{80,80,80},nIn,50){ };
		~NDTMatch_SE();


        vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> getSegmentsFast(
                pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudIn);
        Eigen::Affine3d slam(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud);
        bool matchToSaved(Eigen::Affine3d &Td_,
                Eigen::Vector3d &pose_ref, int start_index, int iP,
                Eigen::Matrix<double,7,7> &Ccl, int stop_index, int target_index);
        int find_start(std::map<int,perception_oru::NDTHistogram*>::iterator pi, float max_size);
        void loadMap(perception_oru::NDTMap **map,
                std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> input_clouds,
                float sensor_range=100);
		Eigen::Matrix<double,7,7> getPoseInformation(Eigen::Affine3d T,
                perception_oru::NDTMap ***m1, perception_oru::NDTMap*** m2, bool inverse);
        void print_vals();

        bool useSaved=false;
        unsigned int NumInputs;
		perception_oru::NDTMap ***map;
		perception_oru::NDTMap ***mapLocal_prev;
		perception_oru::NDTMap ***mapLocal;
		vector<float> resolutions,size;
		vector<int> resolutions_order;
        Eigen::Affine3d T, Td;
        std::map<int,perception_oru::NDTHistogram*> key_hists;
        unsigned int num_clouds=0;
        float max_size;
		perception_oru::NDTMatcherD2D_SE matcher;


		float sensor_range=100;
    private:
        PoseOptimizer pose_graph;
        int last_loop_close_id=0;
		bool firstRun;
#ifdef GL_VISUALIZE
    public:
        NDTViz *viewer = NULL;
        void visualize();
    private:
        void visualize_thread();
#endif
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
#endif/*SE_NDT*/
