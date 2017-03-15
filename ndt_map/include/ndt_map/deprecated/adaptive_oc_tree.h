/*
 * Software License Agreement (BSD License)
 *
 *  Copyright (c) 2010, AASS Research Center, Orebro University.
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above
 *     copyright notice, this list of conditions and the following
 *     disclaimer in the documentation and/or other materials provided
 *     with the distribution.
 *   * Neither the name of AASS Research Center nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 *  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 *  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 *  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 *  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *  CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 *  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 *  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 *
 */
#ifndef ADAPTIVE_OCT_TREE_HH
#define ADAPTIVE_OCT_TREE_HH

#include <ndt_map/oc_tree.h>
#include <ndt_map/ndt_cell.h>
#include <vector>

namespace lslgeneric
{

/** \brief Implements an OctTree with adaptive leaf splitting as a post process step
    \details The OctTree is split unitl a conservative initial size. After this, a post
processing step is used to determine if the gaussians in each leaf are well-fitting.
At the moment one heuristic option is available - using the residual squares sum and
the residual variance. A second option using the omnibus normality test will be added
soon.
 */
template<typename PointT>
class AdaptiveOctTree : public OctTree<PointT>
{
protected:
    std::vector<OctTree<PointT>*> splitTree(OctTree<PointT> *leaf);
    std::vector<OctTree<PointT>*> myTreeLeafs;
    virtual void computeTreeLeafs();

    double computeResidualSquare(NDTCell<PointT> *cell);
    double computeDornikHansen(NDTCell<PointT> *cell);

    bool useDornikHansen;
    bool useFlatness;
    double RSS_THRESHOLD, DH_SIGNIFICANCE_LVL, FLAT_FACTOR;
    bool parametersSet;

public:

    ///dummy default constructor
    AdaptiveOctTree();
    ///creates an oct tree node with known center and size
    AdaptiveOctTree(pcl::PointXYZ center, double xsize, double ysize,
                    double zsize, NDTCell<PointT>* type, OctTree<PointT> *_parent=NULL, unsigned int _depth=0);
    virtual ~AdaptiveOctTree();

    virtual void postProcessPoints();
    virtual SpatialIndex<PointT>* clone();

    ///use this to set the parameters for the OctTree. If not called before creating the
    ///first leaf, default parameters will be used. \note be careful, remember that the parameters are static, thus global
    void setParameters(bool _useDornikHansen = false,
                       bool _useFlatness = true,
                       double _RSS_THRESHOLD = 1000,
                       double _DH_SIGNIFICANCE_LVL = 0.5,
                       double _MIN_CELL_SIZE = 1,
                       double _FLAT_FACTOR = 10,
                       double _BIG_CELL_SIZE = 5,
                       double _SMALL_CELL_SIZE = 5
                      );
    double MIN_CELL_SIZE;

};

} //end namespace

#include<ndt_map/impl/adaptive_oc_tree.hpp>
#endif
