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
#ifndef CELL_HH
#define CELL_HH

#include <vector>
#include <cstdio>
#include <cmath>

namespace lslgeneric
{

/** \brief Base class for a rectangular 3D cell
  * \details The Cell class provides a base for all types of derived cells.
  * The logic behind this is to abstract some of the functionalities that
  * are expected from a spatial cell and eventually be able to interchange
  * cell classes in \ref SpatialIndex implementations.
  */
template<typename PointT>
class Cell
{

protected:
    PointT center_;
    double xsize_, ysize_, zsize_;

public:
    Cell() { }
    inline Cell(PointT &center, double &xsize, double &ysize, double &zsize)
    {
        center_ = center;
        xsize_ = xsize;
        ysize_ = ysize;
        zsize_ = zsize;
    }
    inline Cell(const Cell& other)
    {
        center_ = other.center_;
        xsize_ = other.xsize_;
        ysize_ = other.ysize_;
        zsize_ = other.zsize_;
    }
    virtual ~Cell() { }

    inline void setCenter(const PointT &cn)
    {
        center_ = cn;
    }
    inline void setDimensions(const double &xs, const double &ys, const double &zs)
    {
        xsize_ = xs;
        ysize_ = ys;
        zsize_ = zs;
    }

    inline PointT getCenter() const
    {
        return center_;
    }
    inline void getDimensions(double &xs, double &ys, double &zs) const
    {
        xs = xsize_;
        ys = ysize_;
        zs = zsize_;
    }
    inline bool isInside(const PointT pt) const
    {
        if(pt.x < center_.x-xsize_/2 || pt.x > center_.x+xsize_/2)
        {
            return false;
        }
        if(pt.y < center_.y-ysize_/2 || pt.y > center_.y+ysize_/2)
        {
            return false;
        }
        if(pt.z < center_.z-zsize_/2 || pt.z > center_.z+zsize_/2)
        {
            return false;
        }
        return true;
    }
    virtual double getDiagonal() const
    {
        return std::sqrt(xsize_*xsize_+ysize_*ysize_+zsize_*zsize_);
    }

    ///clone - create an empty object with same type
    inline virtual Cell<PointT>* clone() const
    {
        Cell<PointT> *ret = new Cell<PointT>();
        return(ret);
    }
    ///copy - create the same object as a new instance
    virtual Cell<PointT>* copy() const
    {
        Cell<PointT> *ret = new Cell<PointT>();
        ret->setDimensions(xsize_,ysize_,zsize_);
        ret->setCenter(center_);
        return(ret);
    }

    virtual void addPoint(PointT &pt) { }
};


} //end namespace

//#include<impl/cell.hpp>

#endif
