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
#ifndef LSL_LAZZY_GRID_HH
#define LSL_LAZZY_GRID_HH

#include <ndt_map/spatial_index.h>
#include <ndt_map/ndt_cell.h>
#include <boost/serialization/base_object.hpp>
#include "boost/serialization/serialization.hpp"
#include <boost/serialization/export.hpp>
#include "boost/serialization/array.hpp"
#include "boost/serialization/list.hpp"
#include <boost/serialization/vector.hpp>
#include "stdio.h"
#include <boost/serialization/string.hpp>
#include <boost/serialization/version.hpp>
#include <boost/serialization/split_member.hpp>
#include <boost/filesystem.hpp>
#include <set>
namespace perception_oru
{
class CellVector3d {
public:
  CellVector3d(size_t x=0, size_t y=0, size_t z=0,  NDTCell* init_val=NULL) :
    d1(x), d2(y), d3(z), data(x*y*z, init_val)
  {}
  NDTCell* GetVal(size_t i, size_t j, size_t k, bool &successfull)  {
    if(withinRange(i,j,k)){
      successfull=true;
      return data[i*d2*d3 + j*d3 + k];}
    else{
      successfull=false;
      return(NDTCell*) NULL;
    }
  }
  NDTCell* GetVal(size_t i, size_t j, size_t k)  {
      bool not_used;
      return GetVal( i,  j,  k, not_used);
  }
  void SetVal(size_t i, size_t j, size_t k,  NDTCell* cell){
    if(withinRange(i,j,k))
      data[i*d2*d3 + j*d3 + k]=cell;
  }
private:
  bool withinRange(size_t i, size_t j, size_t k){
    if(i<d1 && i>=0 && j<d2&& j>=0 &&k<d3 &&k>=0)
      return true;
    else return false;
  }
  friend class boost::serialization::access;
  template<class Archive>
  void serialize(Archive & ar, const unsigned int version)//In order to clal this you need to register it to boost using "ar.template register_type<LazyGrid>();"
  {
    ar & d1 & d2 & d3;
    ar & data;
  }
  size_t d1,d2,d3;
  std::vector<NDTCell*> data;
};

/** \brief A spatial index represented as a grid map
    \details A grid map with delayed allocation of cells.
*/
class LazyGrid : public SpatialIndex
{
public:
    unsigned int cloud_index=0;
    int semantic_index;
    std::string finalpath;
    void loadNDTs(int index);
    LazyGrid(double cellSize);
    LazyGrid(double cellSizeX_, double cellSizeY_, double cellSizeZ_);
    LazyGrid(LazyGrid *prot);
    LazyGrid(double sizeXmeters, double sizeYmeters, double sizeZmeters,
             double cellSizeX, double cellSizeY, double cellSizeZ,
             double _centerX, double _centerY, double _centerZ,
             NDTCell *cellPrototype );
    virtual ~LazyGrid();

    bool isValid(const double* const p, NDTCell* cell);
    virtual NDTCell* getCellForPoint(const double* const point);
    virtual NDTCell* addPoint(const double* const point);
    void clearCells();

    struct Darray
    {
        int size[3];
        NDTCell** array;
        Darray(): array(NULL), size{}{}
        void initialize(int x,int y,int z){
            if(x*y*z<1)
                x=y=z=1;
            size[0]=x;
            size[1]=y;
            size[2]=z;
            array=(NDTCell**) calloc(x*y*z, sizeof(NDTCell*));
        };
        void xyzBound(int &x,int &y,int &z)
        {
            x=x%size[0];
            y=y%size[1];
            z=z%size[2];
            while(x<0)x+=size[0];
            while(y<0)y+=size[1];
            while(z<0)z+=size[2];
        }
        NDTCell* operator()(int x, int y, int z)
        {
            xyzBound(x,y,z);
            return array[z+y*size[2]+x*size[2]*size[1]];
        }
        void set(int x, int y, int z, NDTCell* c)
        {
            xyzBound(x,y,z);
            array[z+y*size[2]+x*size[2]*size[1]] = c;
        }
    };
    void setSensorPose(const double *pose);
    bool inRange(const double* const p);
 
    //these two don't make much sense...
    ///iterator through all cells in index, points at the begining
    virtual typename SpatialIndex::CellVectorItr begin();
    virtual typename SpatialIndex::CellVectorConstItr begin() const;
    ///iterator through all cells in index, points at the end
    virtual typename SpatialIndex::CellVectorItr end();
    virtual typename SpatialIndex::CellVectorConstItr end() const;
    virtual int size();
    void deallocateCell(int i,int j, int k);

    ///clone - create an empty object with same type
    virtual SpatialIndex* clone() const;
    virtual SpatialIndex* copy() const;

    ///method to return all cells within a certain radius from a point
    virtual void getNeighbors(const double* const point, const double &radius, std::vector<NDTCell*> &cells);
	///method to return all cells within a certain radius from a point. Clone the cell once and the return shared pointers.
    virtual void getNeighborsShared(const double* const point, const double &radius, std::vector<boost::shared_ptr< NDTCell > > &cells);

    ///sets the cell factory type
    virtual void setCellType(NDTCell *type);

    virtual void setCenter(const double &cx, const double &cy, const double &cz);
    virtual void setSize(const double &sx, const double &sy, const double &sz);
    bool insertCell(NDTCell cell);

    virtual NDTCell* getClosestNDTCell(const double * const pt, bool checkForGaussian=true);
    virtual std::vector<NDTCell*> getClosestNDTCells(const double* const pt, int n_neigh, bool checkForGaussian=true);
    virtual std::vector<boost::shared_ptr< NDTCell > > getClosestNDTCellsShared(const double* const pt, int n_neigh, bool checkForGaussian=true);
    virtual std::vector<NDTCell*> getClosestCells(const double* const pt);

    virtual inline void getCellAt(int indX, int indY, int indZ, NDTCell* &cell){
      cell = dataArray(indX,indY,indZ);
    }
  
    virtual inline void getCellAt(const double* const pt, NDTCell* &cell){
       int indX,indY,indZ;
       this->getIndexForPoint(pt,indX,indY,indZ);
       this->getCellAt(indX,indY,indZ,cell);
    }


  ///automatically allocate the cell if needed (checks that the indexes etc. are correct).
  virtual bool getCellAtAllocate(const double* const t, NDTCell* &cell);

    //FIXME: these two are now not needed anymore
    virtual inline void getNDTCellAt(int indX, int indY, int indZ, NDTCell* &cell){
        cell = (dataArray(indX,indY,indZ));
    }
    virtual inline void getNDTCellAt(const double* const pt, NDTCell* &cell){
			int indX,indY,indZ;
			this->getIndexForPoint(pt,indX,indY,indZ);
			this->getNDTCellAt(indX,indY,indZ,cell);
    }

    void getCellSize(double &cx, double &cy, double &cz);
    void getGridSize(int &cx, int &cy, int &cz);
    void getGridSizeInMeters(double &cx, double &cy, double &cz);
    void getCenter(double &cx, double &cy, double &cz);
    virtual void getIndexForPoint(const double* const t, int &idx, int &idy, int &idz);
    NDTCell* getProtoType()
    {
        return protoType;
    }

    virtual void initialize();

//    NDTCell ****getDataArrayPtr()
//    {
//        return dataArray;
//    }

    ///reads map contents from .jff file
    virtual int loadFromJFF(FILE * jffin);
    bool traceLine(const Eigen::Vector3d &origin, const Eigen::Vector3d &endpoint, const Eigen::Vector3d &diff, const double& maxz, std::vector<NDTCell*> &cells);
    
    void loadCells(int index_start, int index_end);
    void InitializeDefaultValues();
    std::string GetDataString();
    std::string ToString();
    volatile int sensor_pose[3]={0,0,0};
    Darray dataArray;
    double cellSizeX, cellSizeY, cellSizeZ;
    Eigen::Vector3d translation;
    void addNDTCell(NDTCell* cell);
protected:
    bool initialized=false;
    //bool ***linkedCells;
    NDTCell *protoType=NULL;
    std::set<NDTCell*> activeCells;
    bool centerIsSet, sizeIsSet;

    double sizeXmeters, sizeYmeters, sizeZmeters;
    double centerX, centerY, centerZ;
    int sizeX,sizeY,sizeZ;

    virtual bool checkCellforNDT(int indX, int indY, int indZ, bool checkForGaussian=true);
private:
    void dealocateCells(int i, int d, boost::archive::binary_oarchive& oa, unsigned int &min_index);
    LazyGrid(){InitializeDefaultValues();}

    friend class boost::serialization::access;
    template<class Archive>
    void save(Archive& ar,  const unsigned int version) const {
      ar & boost::serialization::base_object<SpatialIndex>(*this);
      ar & protoType;
      ar & sizeX & sizeY & sizeZ;
      for(int i=0; i<sizeX*sizeY*sizeZ; i++)
           ar & dataArray.array[i];
      ar & activeCells;
      ar & centerIsSet & sizeIsSet   & initialized; //Serialize all primitive types
      ar & sizeXmeters & sizeYmeters & sizeZmeters;
      ar & cellSizeX   & cellSizeY   & cellSizeZ;
      ar & centerX     & centerY     & centerZ;
    }
    template<class Archive>
    void load(Archive & ar, const unsigned int version) {
      ar & boost::serialization::base_object<SpatialIndex>(*this);
      ar & protoType;
      ar & sizeX & sizeY & sizeZ;

      dataArray.initialize(sizeX,sizeY,sizeZ);
      for(int i=0; i<sizeX*sizeY*sizeZ; i++)
           ar & dataArray.array[i];
      ar & activeCells;
      ar & centerIsSet & sizeIsSet   & initialized; //Serialize all primitive types
      ar & sizeXmeters & sizeYmeters & sizeZmeters;
      ar & cellSizeX   & cellSizeY   & cellSizeZ;
      ar & centerX     & centerY     & centerZ;
    }

    BOOST_SERIALIZATION_SPLIT_MEMBER()

public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};


}; //end namespace
//#include<ndt_map/impl/lazy_grid.hpp>

#endif
