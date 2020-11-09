#pragma once
#ifndef SE_NDT
#define SE_NDT
#include <ndt_map/ndt_map.h>
#include "se_ndt/ndt_matcher_d2d_se.h"
#include "se_ndt/ndt_histogram.h"
#include "pose_optimizer.h"
#ifdef GL_VISUALIZE
    #include "ndt_visualisation/ndt_viz.h"
#endif
#include <boost/archive/binary_iarchive.hpp>
using namespace std;

class NDTMatch_SE{
    public:
		NDTMatch_SE(initializer_list<float> b,initializer_list<int> c,
                initializer_list<float> d,int nIn,int max_iter);
		NDTMatch_SE(int nIn): NDTMatch_SE({4,0.8},{0,1},{80,80,80},nIn,50){ };
		~NDTMatch_SE();

        void slamSimple(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud);
        Eigen::Affine3d slam(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud);
        bool matchToSaved(Eigen::Affine3d &Td_,
                Eigen::Vector3d &pose_ref, int start_index,
                Eigen::Matrix<double,7,7> &Ccl, int stop_index, int current);
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
        PoseOptimizer pose_graph;
    private:
        int last_loop_close_id=0;
		bool firstRun;
#ifdef GL_VISUALIZE
    public:
        NDTViz *viewer = NULL;
        std::thread* visualize();
    private:
        void visualize_thread();
#endif

void save_map(int index, perception_oru::NDTMap ***map_)
{
        std::ofstream ofs("/tmp/maps/"+std::to_string(index),std::ios::binary);
        boost::archive::binary_oarchive oa(ofs);
        for(unsigned int i=0;i<resolutions.size();i++)
            for(size_t j=0;j<NumInputs;j++)
            {
                oa.template register_type<perception_oru::LazyGrid>();
                oa&(*map_[i][j]);
            }
        ofs.close();
}
void load_map(int index, perception_oru::NDTMap ***map_)
{
        std::ifstream ifs("/tmp/maps/"+std::to_string(index),std::ios::binary);
        boost::archive::binary_iarchive ia(ifs);
        for(unsigned int i=0;i<resolutions.size();i++)
            for(size_t j=0;j<NumInputs;j++)
            {
                ia.template register_type<perception_oru::LazyGrid>();
                ia&(*map_[i][j]);
            }
        ifs.close();
}
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
#endif/*SE_NDT*/
