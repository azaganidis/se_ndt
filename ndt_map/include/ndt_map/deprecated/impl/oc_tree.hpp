#include<pointcloud_vrml/pointcloud_utils.h>

namespace lslgeneric
{

template <typename PointT>
void OctTree<PointT>::setParameters(double _BIG_CELL_SIZE	,
                                    double _SMALL_CELL_SIZE     ,
                                    int _MAX_POINTS		,
                                    int _MAX_DEPTH
                                   )
{

    //defaults
    this->MAX_POINTS		= _MAX_POINTS		;
    this->MAX_DEPTH		= _MAX_DEPTH		;
    this->BIG_CELL_SIZE	= _BIG_CELL_SIZE	;
    this->SMALL_CELL_SIZE	= _SMALL_CELL_SIZE	;
    this->parametersSet_ = true;

    if(!this->leaf_)
    {
        for(unsigned int it=0; it<8; it++)
        {
            if(this->children_[it]!=NULL)
            {
                this->children_[it]->setParameters(_BIG_CELL_SIZE,_SMALL_CELL_SIZE,
                                                   _MAX_POINTS,_MAX_DEPTH);
            }
        }
    }
}
/** returns the child index that should hold the point.
    \note the point is not necessarily within the child's boundaries,
    but is guaranteed to be in the same quadrant (octant?)
  */
template <typename PointT>
size_t OctTree<PointT>::getIndexForPoint(const PointT &pt) const
{
    /** index table: (convenient computations)
      *	id  x	y   z
      * 0   +	+   +
      * 1   +	+   -
      * 2   +	-   -
      * 3   +	-   +
      * 4   -	+   +
      * 5   -	+   -
      * 6   -	-   -
      * 7   -	-   +
      */

    size_t ret = 0;
    if(pt.x < myCell_->getCenter().x)
    {
        ret += 4;
    }
    if(pt.y < myCell_->getCenter().y)
    {
        ret += 2;
    }
    if(ret%4==0)
    {
        if(pt.z < myCell_->getCenter().z)
        {
            ret += 1;
        }
    }
    else
    {
        if(pt.z > myCell_->getCenter().z)
        {
            ret += 1;
        }
    }
    return ret;
}

/** empty constructor
  */
template <typename PointT>
OctTree<PointT>::OctTree()
{

    this->parent_=NULL;
    this->leaf_=true;
    this->depth_=0;
    for(unsigned int it=0; it<8; it++)
    {
        this->children_[it]=NULL;
    }
    myCell_ = NULL;
    leafsCached_ = false;
    this->setParameters();
}

/** parametrized constructor
  */
template <typename PointT>
OctTree<PointT>::OctTree(PointT center, double xsize, double ysize,
                         double zsize, NDTCell<PointT>* type, OctTree<PointT> *_parent, unsigned int _depth)
{

    parent_=_parent;
    leaf_=true;
    depth_=_depth;
    for(unsigned int it=0; it<8; it++)
    {
        children_[it]=NULL;
    }
    NDTCell<PointT>* tmp = dynamic_cast<NDTCell<PointT>*>(type->clone());
    if(tmp==NULL)
    {
        //ERR("dynamic cast of cell failed!!\n");
        return;
    }

    tmp->setCenter(center);
    tmp->setDimensions(xsize,ysize,zsize);
    tmp->points_.clear();

    myCell_ = tmp;
    leafsCached_ = false;
    this->setParameters();
}

/** destructor, deletes all pointers starting from parent and going down
  */
template <typename PointT>
OctTree<PointT>::~OctTree()
{

    if(!leaf_)
    {
        for(unsigned int it=0; it<8; it++)
        {
            if(children_[it]!=NULL)
            {
                delete children_[it]; //calls destructor of child, so we are ok
                children_[it]=NULL;
            }
        }
    }
    delete myCell_;
}

/** adds a point to the index. Iterates down the tree and if necessary
  creates new leafs and splits current ones.
  \note at the moment root is not grown in case of points outside!
  */
template <typename PointT>
Cell<PointT>* OctTree<PointT>::addPoint(const PointT &point_c)
{

    PointT point = point_c;
    if(std::isnan(point.x) ||std::isnan(point.y) ||std::isnan(point.z))
    {
        return NULL;
    }
    leafsCached_ = false;
    if(this->leaf_)
    {
        double xs,ys,zs;
        myCell_->getDimensions(xs,ys,zs);

        double cellSize = (xs+ys+zs)/3.; //average for now

        if(myCell_->points_.size()<(unsigned int)MAX_POINTS && cellSize <= BIG_CELL_SIZE)
        {
            if(!myCell_->isInside(point))
            {
                //DBG(1,"OctTree: addPoint (%lf,%lf,%lf) not in boundary!\n",point.x,point.y,point.z);
                return NULL;
            }
            myCell_->points_.push_back(point);
        }
        else
        {
            if(depth_>(unsigned int) MAX_DEPTH || cellSize <= 2*SMALL_CELL_SIZE )
            {
                //TSV: we have to be sure we won't violate the space constraint if we split
                //just store point, we can't split any more
                if(!myCell_->isInside(point))
                {
                    //DBG(1,"OctTree: addPoint (%lf,%lf,%lf) not in boundary!\n",point.x,point.y,point.z);
                    return NULL;
                }
                myCell_->points_.push_back(point);
                return myCell_;
            }

            PointT myCenter = myCell_->getCenter();

            //branch leaf
            for(int it=0; it<8; it++)
            {
                PointT newCenter;

                //computes the center_ of the it'th child
                newCenter.x = (myCenter.x + pow(-1.,it/4)*xs/4.);
                newCenter.y = (myCenter.y + pow(-1.,it/2)*ys/4.);
                newCenter.z = (myCenter.z + pow(-1.,(it+1)/2)*zs/4.);

                children_[it] = new OctTree<PointT>(newCenter,xs/2,ys/2,
                                                    zs/2, myCell_, this, depth_+1);
                children_[it]->setParameters(BIG_CELL_SIZE,SMALL_CELL_SIZE,
                                             MAX_POINTS,MAX_DEPTH);
            }
            //add current points
            for(unsigned int jt=0; jt<myCell_->points_.size(); jt++)
            {
                size_t ind = getIndexForPoint(myCell_->points_[jt]);
                children_[ind]->addPoint(myCell_->points_[jt]);
            }
            //finally add the new point
            size_t ind = getIndexForPoint(point);
	    Cell<PointT>* ptcell = children_[ind]->addPoint(point);
            this->leaf_=false;
            this->myCell_->points_.clear();
	    return ptcell;
        }
    }
    else
    {
        //pass down to correct child
        size_t ind = getIndexForPoint(point);
        return children_[ind]->addPoint(point);
    }
}

/** returns the cell that should hold the point
  */
template <typename PointT>
Cell<PointT>* OctTree<PointT>::getCellForPoint(const PointT &point)
{

    OctTree<PointT>* pointLeaf = this->getLeafForPoint(point);
    return (pointLeaf==NULL) ? NULL : pointLeaf->myCell_;

}

/** finds the leaf that should hold the point
  */
template <typename PointT>
OctTree<PointT>* OctTree<PointT>::getLeafForPoint(const PointT &point)
{

    if(this->leaf_ && myCell_!= NULL)
    {
        if(myCell_->isInside(point))
        {
            return this;
        }
    }
    else
    {
        size_t ind = getIndexForPoint(point);
        if(children_[ind]!=NULL)
        {
            return children_[ind]->getLeafForPoint(point);
        }
    }
    return NULL;

}

/** finds all children cells that are located at current leafs
  */
template <typename PointT>
void OctTree<PointT>::computeLeafCells()
{
    if(this->isLeaf())
    {
        myLeafs_.push_back(this->myCell_);
        return;
    }

    myLeafs_.clear();
    std::vector<OctTree<PointT>*> next;
    next.push_back(this);

    while(next.size()>0)
    {
        OctTree<PointT> *cur = next.front();
        if(cur!=NULL)
        {
            if(cur->isLeaf())
            {
                myLeafs_.push_back(cur->myCell_);
            }
            else
            {
                for(int i=0; i<8; i++)
                {
                    OctTree<PointT>* tmp = cur->getChild(i);
                    if(tmp!=NULL)
                    {
                        next.push_back(tmp);
                    }
                }
            }
        }
        next.erase(next.begin());
    }
}

/** sets the cell factory type
  */
template <typename PointT>
void OctTree<PointT>::setCellType(Cell<PointT> *type)
{

    myCell_ = dynamic_cast<NDTCell<PointT>*>(type->clone());
    if(myCell_ == NULL)
    {
        //cast failed, it's not a derivative of oct cell
        myCell_ = new NDTCell<PointT>();
    }

}

/** iterator poining to the first leaf cell.
  * \note recomputes the vector when changes have occured!
  */
template <typename PointT>
typename SpatialIndex<PointT>::CellVectorItr OctTree<PointT>::begin()
{
    if(!leafsCached_)
    {
        myLeafs_.clear();
        computeLeafCells();
    }
    leafsCached_ = true;
    return myLeafs_.begin();
}

/** iterator poining at last leaf cell
  */
template <typename PointT>
typename SpatialIndex<PointT>::CellVectorItr OctTree<PointT>::end()
{
    if(!leafsCached_)
    {
        myLeafs_.clear();
        computeLeafCells();
    }
    leafsCached_ = true;
    return myLeafs_.end();
}

/** returns a new OctTree spatial index
  */
template <typename PointT>
SpatialIndex<PointT>* OctTree<PointT>::clone() const
{
    OctTree<PointT> *tr = new OctTree<PointT>();
    tr->setParameters(BIG_CELL_SIZE,SMALL_CELL_SIZE,MAX_POINTS,MAX_DEPTH);
    if(myCell_ != NULL)
    {
        tr->setCellType(myCell_);
    }
    return tr;
}

/** copies the spatial index
  */
template <typename PointT>
SpatialIndex<PointT>* OctTree<PointT>::copy() const
{
    OctTree<PointT> *tr;
    if(myCell_ != NULL)
    {
        PointT center = myCell_->getCenter();
        double sx,sy,sz;
        myCell_->getDimensions(sx,sy,sz);
        tr = new OctTree<PointT>(center,sx,sy,sz,myCell_);
    }
    else
    {
        tr = new OctTree<PointT>();
    }
    tr->setParameters(BIG_CELL_SIZE,SMALL_CELL_SIZE,MAX_POINTS,MAX_DEPTH);
    return tr;
}

/** sets the center_.
  \note this is not going to re-compute cells that are already in the tree!
  */
template <typename PointT>
void OctTree<PointT>::setCenter(const double &cx, const double &cy, const double &cz)
{
    if(myCell_ == NULL)
    {
        return;
    }
    PointT center_;
    center_.x = cx;
    center_.y = cy;
    center_.z = cz;

    myCell_->setCenter(center_);
}

/** sets the size
  \note this is not going to re-compute cells that are already in the tree!
  */
template <typename PointT>
void OctTree<PointT>::setSize(const double &sx, const double &sy, const double &sz)
{
    if(myCell_ == NULL)
    {
        return;
    }
    myCell_->setDimensions(sx,sy,sz);
}

/** returns all the neighboring cells within a radius
   \param point the location around which we are looking for neighbors. The point must be inside the boundaries of a current leaf!
   \param radius the neighbor tolerance
   \param cells output
  */
template <typename PointT>
void OctTree<PointT>::getNeighbors(const PointT &point, const double &radius, std::vector<Cell<PointT>*> &cells)
{

    cells.clear();
    //first find the leaf that contains the point
    OctTree<PointT> *target = this->getLeafForPoint(point);
    if(target==NULL) return;

    OctTree<PointT> *mparent = target->parent_;
    OctTree<PointT> *mthis = target;
    std::vector<OctTree<PointT>*> toExpand;
    //check if any of the siblings intersect the sphere (key,constraint)

    while(mparent!=NULL)
    {
        for(unsigned int it=0; it<8; it++)
        {
            if(mparent->children_[it] == NULL) continue;
            if(mparent->children_[it]->intersectSphere(point,radius)
                    && mparent->children_[it]!=mthis )
            {
                //if yes, add them to the list to expand
                toExpand.push_back(mparent->children_[it]);
            }
        }
        //go up to parent
        mthis=mparent;
        mparent=mparent->parent_;
        //for all nodes in list, go down to leafs that intersect
        for(unsigned int nt=0; nt<toExpand.size(); nt++)
        {
            if(toExpand[nt] == NULL )
            {
                //ERR("ERROR in nearest neighbor!!\n");
                continue;
            }

            PointT center_ = (toExpand[nt]->myCell_->getCenter());
            Eigen::Vector3d d;
            d<<center_.x-point.x, center_.y-point.y, center_.z-point.z;
            if(toExpand[nt]->isLeaf() &&
                    d.norm() < radius)
            {
                cells.push_back(toExpand[nt]->myCell_);
            }
            else
            {
                for(unsigned int it=0; it<8; it++)
                {
                    if(toExpand[nt]->children_[it]==NULL) continue;
                    if(toExpand[nt]->children_[it]->intersectSphere(point,radius))
                    {
                        toExpand.push_back(toExpand[nt]->children_[it]);
                    }
                }
            }
        }

        toExpand.clear();
    }

}

/** checks if the tree node intersects the sphere located at center_ and of size radius
  */
template <typename PointT>
bool OctTree<PointT>::intersectSphere(PointT center_, const double &radius) const
{

    PointT mcenter_ = myCell_->getCenter();
    Eigen::Vector3d d;
    d<<center_.x-mcenter_.x, center_.y-mcenter_.y, center_.z-mcenter_.z;
    double dist = d.norm();
    Eigen::Vector3d localRadius;
    myCell_->getDimensions(localRadius(0),localRadius(1),localRadius(2));
    double lRad = localRadius.norm()/2;
    double interDist = lRad+radius;
    return (interDist>dist);

}

template <typename PointT>
Cell<PointT>* OctTree<PointT>::getClosestLeafCell(const PointT &point)
{
    if(this->leaf_ && myCell_!= NULL)
    {
        if(myCell_->isInside(point))
        {
            return myCell_;
        }
    }
    else
    {
        size_t ind = getIndexForPoint(point);
        if(children_[ind]!=NULL)
        {
            return children_[ind]->getClosestLeafCell(point);
        }
        else
        {
            //the leaf we should be in is empty
            //start from here and find closest neighbor
            double minDist = INT_MAX;
            int index = -1;
            for(int i=0; i<8; i++)
            {
                if(children_[i]==NULL) continue;
                PointT center = children_[i]->myCell_->getCenter();
                Eigen::Vector3d d;
                d <<center.x-point.x, center.y-point.y, center.z-point.z;
                double dist = d.norm();

                if(dist<minDist)
                {
                    index = i;
                    minDist = dist;
                }
            }
            //std::cout<<"minDist"<<minDist<<std::endl;
            if(index>=0 && index<8)
            {
                return children_[index]->getClosestLeafCell(point);
            }
        }
    }
    return myCell_;
}

template <typename PointT>
NDTCell<PointT>* OctTree<PointT>::getClosestNDTCell(const PointT &point)
{
    if(this->leaf_)
    {
        //we have reached the bottom of the tree.
        //if we have a valid ndt cell, return it.
        if(myCell_->isInside(point))
        {
            NDTCell<PointT>* nd = dynamic_cast<NDTCell<PointT>*>(myCell_);
            if(nd!=NULL)
            {
                if(nd->hasGaussian_)
                {
                    return nd;
                }
            }
        }
    }

    //we go down the tree recursively
    size_t ind = getIndexForPoint(point);
    if(children_[ind]!=NULL)
    {
        return children_[ind]->getClosestNDTCell(point);
    }


    //the leaf we should be in is empty
    //iterate through all leafs connected to current node, find closest ndt cell
    OctTree<PointT> *my_parent = this->parent_;
    OctTree<PointT> *me = this;
    typename SpatialIndex<PointT>::CellVectorItr it;
    NDTCell<PointT> *closest = NULL, *temp = NULL;
    double minDist = INT_MAX, dist = INT_MAX;

    while(true)
    {
        it = me->begin();
        while(it!=me->end())
        {
            temp = dynamic_cast<NDTCell<PointT>*> (*it);
            if(temp!=NULL)
            {
                if(temp->hasGaussian_)
                {
                    dist = lslgeneric::geomDist<PointT>(temp->getCenter(),point);
                    if(dist < minDist)
                    {
                        minDist = dist;
                        closest = temp;
                    }
                }
            }
            it++;
        }
        if(closest!=NULL)
        {
//	    cout<<"got it!\n";
            break;
        }
        if(my_parent != NULL)
        {
            me = my_parent;
            my_parent = me->parent_;
        }
        else
        {
            //nothing more can be done...
            break;
        }
    }

    return closest;
}
}
