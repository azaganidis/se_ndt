#include <climits>
#include <pcl/common/distances.h>
#include <ndt_map/cell_vector.h>

namespace perception_oru
{

CellVector::CellVector():mp(new pcl::PointCloud<pcl::PointXYZ>())
{
    pcl::console::setVerbosityLevel(pcl::console::L_ERROR);
    protoType= new NDTCell();
    treeUpdated = false;
}

CellVector::CellVector(NDTCell* cellPrototype):mp(new pcl::PointCloud<pcl::PointXYZ>())
{
    protoType = cellPrototype->clone();
    treeUpdated = false;
}

CellVector::CellVector(const CellVector& other):protoType(NULL)
{
    for(unsigned int i =0; i< other.activeCells.size(); i++)
    {
        NDTCell* r = *std::next(other.activeCells.begin(), i);
        if(r == NULL) continue;
        for(size_t i=0; i<r->points_.size(); i++)
        {
            this->activeCells.insert(r->copy());
        }
    }
}

CellVector::~CellVector()
{
    //go through all cells and delete the non-NULL ones
    for(SpatialIndex::CellVectorItr it=activeCells.begin(); it!=activeCells.end(); ++it)
    {
        if((*it)!=NULL)
        {
            delete (*it);
        }
    }
}

NDTCell* CellVector::getCellForPoint(const double* const &point)
{
    NDTCell* ret = NULL;
    if (treeUpdated)
    {
        std::vector<int> id;
        std::vector<float> dist;
        int NCELLS = 1;
        id.reserve(NCELLS);
        dist.reserve(NCELLS);
        const pcl::PointXYZ pt(point);
        if(!meansTree.nearestKSearch(pt,NCELLS,id,dist)) return ret;

        ret = *std::next(activeCells.begin(), id[0]);
    }
    else
    {
        float min_dist = std::numeric_limits<float>::max( );
        typename SpatialIndex::CellVectorItr it = this->begin();
        while(it!=this->end())
        {
            float tmp=pcl::squaredEuclideanDistance((*it)->getCenter(), point);
            if (tmp < min_dist)
            {
                min_dist = tmp;
                ret = (*it);
            }
            it++;
        }
    }
    return ret;
}

NDTCell* CellVector::addPoint(const pcl::PointXYZ &point)
{
    return NULL;
    // Do nothing...
}

void CellVector::addCellPoints(pcl::PointCloud<pcl::PointXYZ> pc, const std::vector<size_t> &indices)
{
    auto tmp = protoType->clone();
    activeCells.insert(tmp);
    for (size_t i = 0; i < indices.size(); i++) {
        tmp->addPoint(pc[indices[i]]); // Add the point to the cell.
    }
    treeUpdated = false;
}

void CellVector::addCell(NDTCell* cell)
{
    activeCells.insert(cell);
}

void CellVector::addNDTCell(NDTCell* cell)
{
    this->addCell(cell);
}

typename SpatialIndex::CellVectorItr CellVector::begin()
{
    //cout<<"active cells "<<activeCells.size()<<endl;
    return activeCells.begin();
}

typename SpatialIndex::CellVectorConstItr CellVector::begin() const
{
    //cout<<"active cells "<<activeCells.size()<<endl;
    return activeCells.begin();
}

typename SpatialIndex::CellVectorItr CellVector::end()
{
    return activeCells.end();
}

typename SpatialIndex::CellVectorConstItr CellVector::end() const
{
    return activeCells.end();
}

int CellVector::size()
{
    return activeCells.size();
}

SpatialIndex* CellVector::clone() const
{
    return new CellVector();
}

SpatialIndex* CellVector::copy() const
{
    //std::cout<<"COPY CELL VECTOR\n";
    //assert(false); // This needs to be updated!
    CellVector *ret = new CellVector();
    for(unsigned int i =0; i< activeCells.size(); i++)
    {
        NDTCell* r = (*std::next(activeCells.begin(), i))->copy();
        if(r == NULL) continue;
        for(size_t i=0; i<r->points_.size(); i++)
        {
            ret->activeCells.insert(r->copy());
        }
    }
    return ret;
}

void CellVector::getNeighbors(const pcl::PointXYZ &point, const double &radius, std::vector<NDTCell*> &cells)
{

    if (treeUpdated)
    {
        std::vector<int> id;
        std::vector<float> dist;
        int NCELLS = 4;
        id.reserve(NCELLS);
        dist.reserve(NCELLS);
        const pcl::PointXYZ pt(point);

        if(!meansTree.nearestKSearch(pt,NCELLS,id,dist)) return;

        for(int i=0; i<NCELLS; i++)
        {
            NDTCell* tmp = *std::next(activeCells.begin(), id[i]);
            if (tmp != NULL)
                cells.push_back(tmp);
        }
    }
    else
    {
        float radius_sqr = radius*radius;
        typename SpatialIndex::CellVectorItr it = this->begin();
        while(it!=this->end())
        {
            float tmp=pcl::squaredEuclideanDistance((*it)->getCenter(), point);
            if (tmp < radius_sqr)
            {
                cells.push_back(*it);
            }
        }
    }
}

void CellVector::initKDTree()
{

    NDTCell* ndcell = NULL;
    pcl::PointXYZ curr;
    Eigen::Vector3d m;
    pcl::PointCloud<pcl::PointXYZ> mc;

    for(SpatialIndex::CellVectorItr it=activeCells.begin(); it!=activeCells.end(); ++it)
    {
        ndcell = (*it);
        if(ndcell == NULL) continue;
        if(!ndcell->hasGaussian_) continue;
        m = ndcell->getMean();
        curr.x = m(0);
        curr.y = m(1);
        curr.z = m(2);
        mc.push_back(curr);
    }

    if(mc.points.size() > 0)
    {
        *mp = mc;
        meansTree.setInputCloud(mp);
    }

    //++++++++++++++++++treeUpdated = true;
}

void CellVector::setCellType(NDTCell *type)
{
    if(type!=NULL)
    {
        protoType = type->clone();
    }
}

NDTCell* CellVector::getClosestNDTCell(const double* const point)
{
    return getCellForPoint(point);
}

std::vector<NDTCell*> CellVector::getClosestNDTCells(const pcl::PointXYZ &point, double &radius)
{
    std::vector<NDTCell*> ret;
    getNeighbors(point, radius, ret);
    return ret;
}

NDTCell* CellVector::getCellIdx(unsigned int idx) const
{
    if (idx >= activeCells.size())
        return NULL;
    perception_oru::SpatialIndex::CellVectorItr it = activeCells.begin();
    std::advance(it,idx);
    return *it;
}

void CellVector::cleanCellsAboveSize(double size)
{
    //clean cells with variance more then x meters in any direction
    Eigen::Vector3d evals;
    perception_oru::SpatialIndex::CellVectorItr it = this->begin();
    perception_oru::SpatialIndex::CellVectorItr it_tmp;
    while(it!=this->end())
    {
        perception_oru::NDTCell *ndcell = (*it);
        if(ndcell != NULL)
        {
            if(ndcell->hasGaussian_)
            {
                evals = ndcell->getEvals();
                if(sqrt(evals(2)) < size)
                {
                    it++;
                    continue;
                }
                //std::cout<<"rem cell at "<<ndcell->getMean().transpose()<<" evals are "<<evals.transpose()<<std::endl;
                ndcell->hasGaussian_ = false;
            }
            delete ndcell;
            ndcell = NULL;
        }
        it_tmp = it;
        it_tmp--;
        this->activeCells.erase(it);
        it = it_tmp;
        it++;
    }

}
int CellVector::loadFromJFF(FILE * jffin)
{
    NDTCell prototype_;
    if(fread(&prototype_, sizeof(NDTCell), 1, jffin) <= 0)
    {
        JFFERR("reading prototype_ failed");
    }
    protoType = prototype_.clone();
    // load all cells
    while (1)
    {
        if(prototype_.loadFromJFF(jffin) < 0)
        {
            if(feof(jffin))
            {
                break;
            }
            else
            {
                JFFERR("loading cell failed");
            }
        }

        if(!feof(jffin))
        {
            // std::cout << prototype_.getOccupancy() << std::endl; /* for debugging */
        }
        else
        {
            break;
        }
        //initialize cell
        this->addCell(prototype_.copy());
    }

    this->initKDTree();

    return 0;
}

}
