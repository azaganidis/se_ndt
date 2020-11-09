#pragma once
#ifndef POSE_OPTIMIZER_HH
#define POSE_OPTIMIZER_HH
#include "ceres/ceres.h"
#include "types.h"
class PoseOptimizer{
    public:
        PoseOptimizer():quaternion_local_parameterization(new ceres::EigenQuaternionParameterization){};
        ~PoseOptimizer(){write("pose_graph_out.txt");};
        void addPose(Eigen::Affine3d T, int id);
        void addConstraint(Eigen::Affine3d &Td, int id0,int id1,
                Eigen::Matrix<double,7,7> infoM);
        bool solve();
        bool write(const std::string& filename);
        Eigen::Vector3d get();
        Eigen::Affine3d getT();
        template<typename T>
        std::vector<int> get_in_range(std::map<int,T> &key_hists,
                float meters,int idiff) {
            std::vector<int> result;
            int idlast=poses.rbegin()->first;
            for(auto it=key_hists.begin(); it!=key_hists.end();++it)
            {
                int i=it->first;
                if(idlast-i<idiff)
                    break;
                float distance= poses[idlast].distance(poses[i]);
                if(distance<meters)
                    result.push_back(i);
            }
            return result;
        }

        float distance(int i);
        Eigen::Vector3d get(int i);
        std::vector<std::pair<double*,double*> > forGL;
        Eigen::Affine3d getT(int i);
    private:
        bool isFirst=true;
        ceres::LocalParameterization* quaternion_local_parameterization;
        MapOfPoses poses;
        VectorOfConstraints constraints;
        ceres::Problem problem;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};
#endif
