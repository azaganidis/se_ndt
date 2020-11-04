// Ceres Solver - A fast non-linear least squares minimizer
// Copyright 2016 Google Inc. All rights reserved.
// http://ceres-solver.org/
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
// * Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the following disclaimer.
// * Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the following disclaimer in the documentation
//   and/or other materials provided with the distribution.
// * Neither the name of Google Inc. nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.
//
// Author: vitus@google.com (Michael Vitus)

#include <iostream>
#include <fstream>
#include <string>

#include "pose_graph_3d_error_term.h"
#include "pose_optimizer.h"


void PoseOptimizer::addPose(Eigen::Affine3d &T,int id){
    if (poses.find(id) != poses.end()) {
        std::cerr << "Duplicate vertex with ID: " << id;
        return;
    }
    Pose3d ps;
    ps(T);
    poses[id]=ps;
}
void PoseOptimizer::addConstraint(Eigen::Affine3d &Td, int id0, int id1,
        Eigen::Matrix<double,7,7> infoM){
    Constraint3d c;
    c(id0,id1,Td,infoM);
    constraints.push_back(c);
    MapOfPoses::iterator pose_begin_iter = poses.find(id0);
    MapOfPoses::iterator pose_end_iter = poses.find(id1);
    const Eigen::Matrix<double, 6, 6> sqrt_information =c.information.llt().matrixL();
    // Ceres will take ownership of the pointer.
    ceres::CostFunction* cost_function =
        PoseGraph3dErrorTerm::Create(c.t_be, sqrt_information);

    problem.AddResidualBlock(cost_function, NULL,
                              pose_begin_iter->second.p.data(),
                              pose_begin_iter->second.q.coeffs().data(),
                              pose_end_iter->second.p.data(),
                              pose_end_iter->second.q.coeffs().data());

    problem.SetParameterization(pose_begin_iter->second.q.coeffs().data(),
                                 quaternion_local_parameterization);
    problem.SetParameterization(pose_end_iter->second.q.coeffs().data(),
                                 quaternion_local_parameterization);
    forGL.push_back(std::make_pair(pose_begin_iter->second.p.data(),
                pose_end_iter->second.p.data()));
}

// Returns true if the solve was successful.
bool PoseOptimizer::solve() {
    if(isFirst){
        isFirst=false;
        MapOfPoses::iterator pose_start_iter = poses.begin();
        problem.SetParameterBlockConstant(pose_start_iter->second.p.data());
        problem.SetParameterBlockConstant(pose_start_iter->second.q.coeffs().data());
    }
    ceres::Solver::Options options;
    options.max_num_iterations = 200;
    options.linear_solver_type = ceres::SPARSE_NORMAL_CHOLESKY;

    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);

    std::fstream optLog;
    optLog.open("opt_log.txt",std::istream::out);
    optLog << summary.FullReport() << '\n';
    optLog.close();

    return summary.IsSolutionUsable();
}

bool PoseOptimizer::write(const std::string& filename) 
{
    solve();
    std::fstream outfile;
    outfile.open(filename.c_str(), std::istream::out);
    if (!outfile) {
        std::cerr<< "Error opening the file: " << filename;
        return false;
    }
    Eigen::Affine3d calib, calib_inv;
    calib.matrix() << 4.276802385584e-04,-9.999672484946e-01,-8.084491683471e-03,-1.198459927713e-02,-7.210626507497e-03,8.081198471645e-03,-9.999413164504e-01,-5.403984729748e-02,9.999738645903e-01,4.859485810390e-04,-7.206933692422e-03,-2.921968648686e-01,0,0,0,1;
    calib_inv=calib.inverse();
    Eigen::IOFormat ofrmt(Eigen::StreamPrecision, Eigen::DontAlignCols, " "," ","","","","");
    for (MapOfPoses::iterator poses_iter = poses.begin();poses_iter != poses.end();
            ++poses_iter) {
        MapOfPoses::value_type& pair = *poses_iter;
        //outfile<<(pair.second.getT()*calib_inv).matrix().topRows(3).format(ofrmt)<<std::endl;
        //outfile<<pair.second.getT().matrix().topRows(3).format(ofrmt)<<std::endl;
        outfile<<(calib*pair.second.getT()*calib_inv).matrix().topRows(3).format(ofrmt)<<std::endl;
    }
    return true;
}
float PoseOptimizer::distance(int i){
    Pose3d *last=&poses.rbegin()->second;
    return last->distance(poses[i]);
}
std::vector<int> PoseOptimizer::get_in_range(std::vector<int> &key_poses, float meters,int idiff)
{
    std::vector<int> result;
    int idlast=poses.rbegin()->first;
    for(int i: key_poses)
    {
        if(idlast-i<idiff)
            break;
        float distance= poses[idlast].distance(poses[i]);
        if(distance<meters)
            result.push_back(i);
    }
    return result;
}
Eigen::Vector3d PoseOptimizer::get(int i){
    return poses[i].p;
}
Eigen::Vector3d PoseOptimizer::get(){
    return poses.rbegin()->second.p;
}
Eigen::Affine3d PoseOptimizer::getT(int i){
    return poses[i].getT();
}
Eigen::Affine3d PoseOptimizer::getT(){
    return poses.rbegin()->second.getT();
}
