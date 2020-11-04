#pragma once
#ifndef SE_NDT
#define SE_NDT
#include <pcl/point_cloud.h>
#include <pcl/common/io.h>
#include <pcl/common/geometry.h>
#include <ndt_map/ndt_map.h>
#include <ndt_map/lazy_grid.h>
#include <ndt_map/ndt_cell.h>
#include <ndt_map/pointcloud_utils.h>
#include "se_ndt/ndt_matcher_d2d_se.h"
#include "se_ndt/ndt_histogram.h"
#include <Eigen/StdVector>
#include "pose_optimizer.h"
#ifdef GL_VISUALIZE
    #include "ndt_visualisation/ndt_viz.h"
#endif
using namespace std;

class NDTMatch_SE{
    public:
		NDTMatch_SE(initializer_list<float> b,initializer_list<int> c,
                initializer_list<float> d,int nIn,int max_iter);
		NDTMatch_SE(int nIn): NDTMatch_SE({4,0.8},{0,1},{80,80,80},nIn,50){ };
		~NDTMatch_SE();

        Eigen::Affine3d slam(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud);
        bool matchToSaved(Eigen::Affine3d &Td_,
                Eigen::Vector3d &pose_ref, int start_index, int iP,
                Eigen::Matrix<double,7,7> &Ccl, int stop_index, int target_index);
        int find_start(std::map<int,perception_oru::NDTHistogram*>::iterator pi, float max_size);


		Eigen::Matrix<double,7,7> getPoseInformation(Eigen::Affine3d T,
                perception_oru::NDTMap ***m1, perception_oru::NDTMap*** m2, bool inverse);
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
