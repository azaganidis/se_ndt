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

#ifndef EXAMPLES_CERES_TYPES_H_
#define EXAMPLES_CERES_TYPES_H_

#include <istream>
#include <map>
#include <string>
#include <vector>

#include "Eigen/Core"
#include "Eigen/Geometry"

struct Pose3d {
  Eigen::Vector3d p;
  Eigen::Quaterniond q;
  void operator()(Eigen::Affine3d& input){
    p=input.translation();
    q=input.rotation();
    q.normalize();
  }
  Eigen::Affine3d getT(){
      Eigen::Affine3d out;
      out.linear()=q.matrix();
      out.translation()=p;
      return out;
  }
  float distance(Pose3d &other){
      return (p-other.p).norm();
  }
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

struct Constraint3d {
    int id_begin;
    int id_end;
    Pose3d t_be;
    Eigen::Matrix<double, 6, 6> information;
    void operator()(int id0, int id1, Eigen::Affine3d& pose, Eigen::Matrix<double,7,7>& info)
    {
        id_begin=id0;
        id_end=id1;
        t_be(pose);
        information=info.block<6,6>(0,0);
    }
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

typedef std::vector<Constraint3d, Eigen::aligned_allocator<Constraint3d> >
    VectorOfConstraints;

typedef std::map<int, Pose3d, std::less<int>,
                 Eigen::aligned_allocator<std::pair<const int, Pose3d> > >
    MapOfPoses;
#endif  // EXAMPLES_CERES_TYPES_H_
