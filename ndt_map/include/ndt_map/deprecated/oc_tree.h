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
#ifndef OCT_TREE_HH
#define OCT_TREE_HH

#include <ndt_map/spatial_index.h>
#include <ndt_map/ndt_cell.h>

#include <vector>
#include <cstdio>
#include <Eigen/Core>

namespace lslgeneric
{

template <typename PointT>
class AdaptiveOctTree;
/** \brief An Oct Tree data structure for storing 3D points
  * \details This is an implementation of a \ref SpatialIndex that splits space
  * using a tree structure. Each node has 8 children that do not necessarily have
  * associated points. Points are stored using \ref NDTCell cells. When inserting
  * points the tree is split until the size BIG_CELL is reached. Further splits occur
  * when the number of points in a leaf goes above the MAX_POINTS threshold
  */
template <typename PointT>
class OctTree : public SpatialIndex<PointT>
{
protected:
    OctTree* parent_;
    OctTree *children_[8];
    NDTCell<PointT> *myCell_;

    std::vector< Cell<PointT>* > myLeafs_;
    unsigned int depth_;
    bool leaf_;
    bool leafsCached_;

    ///checks in which child node a point would belong
    virtual size_t getIndexForPoint(const PointT &pt) const;

    ///fills in the leafs vector when needed
    virtual void computeLeafCells();


public:
    //--- OctTree Parameters ---//
    /// @param maximum depth of the tree, after which no more splits
    int MAX_DEPTH;
    /// @param number of points after which to split cell
    int MAX_POINTS;
    /// @param at this level do not split any more
    double SMALL_CELL_SIZE;
    /// @param split tree up to this size before allocating a cell
    double BIG_CELL_SIZE;
    bool parametersSet_;

    ///dummy default constructor
    OctTree();
    ///creates an oct tree node with known center and size
    OctTree(PointT center, double xsize, double ysize,
            double zsize, NDTCell<PointT>* type, OctTree<PointT> *_parent=NULL, unsigned int _depth=0);
    virtual ~OctTree();

    ///use this to set the parameters for the OctTree - will *only* apply to leafs of the current node.
    void setParameters(double _BIG_CELL_SIZE	=4,
                       double _SMALL_CELL_SIZE   =0.5 ,
                       int _MAX_POINTS		=10,
                       int _MAX_DEPTH		=20
                      );

    ///add a point to the index
    virtual Cell<PointT>* addPoint(const PointT &point);

    ///returns a pointer to the cell containing the point or NULL if not found
    virtual Cell<PointT>* getCellForPoint(const PointT &point);
    inline virtual Cell<PointT>* getMyCell()
    {
        return myCell_;
    }
    virtual OctTree<PointT>* getLeafForPoint(const PointT &point);

    ///sets the prototype for a cell
    virtual void setCellType(Cell<PointT> *type);

    ///iterator through all cells in index, points at the begining
    virtual typename SpatialIndex<PointT>::CellVectorItr begin();
    ///iterator through all cells in index, points at the end
    virtual typename SpatialIndex<PointT>::CellVectorItr end();

    ///recursively print the tree
    void print();

    ///returns a child at the specified index
    inline OctTree<PointT>* getChild(int idx)
    {
        if(idx<8 && idx>=0)
        {
            return children_[idx];
        }
        return NULL;
    }
    inline bool isLeaf()
    {
        return leaf_;
    }

    /// cloning method inherited from spatial index
    virtual SpatialIndex<PointT>* clone() const;
    virtual SpatialIndex<PointT>* copy() const;

    virtual void setCenter(const double &cx, const double &cy, const double &cz);
    virtual void setSize(const double &sx, const double &sy, const double &sz);

    virtual void getNeighbors(const PointT &point, const double &radius, std::vector<Cell<PointT>*> &cells);
    virtual bool intersectSphere(const PointT center, const double &radius) const;

    virtual Cell<PointT>* getClosestLeafCell(const PointT &pt);
    virtual NDTCell<PointT>* getClosestNDTCell(const PointT &pt);

    friend class AdaptiveOctTree<PointT>;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

} //end namespace

#include<ndt_map/impl/oc_tree.hpp>

#endif
