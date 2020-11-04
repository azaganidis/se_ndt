#include <cstring>
#include <sys/stat.h>
#include <cstdio>
#include <ndt_map/lazy_grid.h>

#ifdef BIN_ARCHIVE
#include <boost/archive/binary_iarchive.hpp>
#else
#include <boost/archive/text_iarchive.hpp>
#endif
#include <boost/iostreams/filter/zlib.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>
#define JFFERR(x) std::cerr << x << std::endl; return -1;

namespace perception_oru
{

LazyGrid::LazyGrid(double cellSize)  : protoType(NULL)
{
    initialized = false;
    centerIsSet = false;
    sizeIsSet = false;
    cellSizeX = cellSizeY = cellSizeZ = cellSize;

}

  LazyGrid::LazyGrid(double cellSizeX_, double cellSizeY_, double cellSizeZ_)
{
    initialized = false;
    centerIsSet = false;
    sizeIsSet = false;
    cellSizeX = cellSizeX_;
    cellSizeY = cellSizeY_;
   cellSizeZ = cellSizeZ_;

}


LazyGrid::LazyGrid(LazyGrid *prot)
{

    sizeXmeters = prot->sizeXmeters;
    sizeYmeters = prot->sizeYmeters;
    sizeZmeters = prot->sizeZmeters;

    cellSizeX = prot->cellSizeX;
    cellSizeY = prot->cellSizeY;
    cellSizeZ = prot->cellSizeZ;

    sizeX = abs(ceil(sizeXmeters/cellSizeX));
    sizeY = abs(ceil(sizeYmeters/cellSizeY));
    sizeZ = abs(ceil(sizeZmeters/cellSizeZ));

    centerX = prot->centerX;
    centerY = prot->centerY;
    centerZ = prot->centerZ;

    protoType = prot->protoType->clone();
    initialize();

}

LazyGrid::LazyGrid(double _sizeXmeters, double _sizeYmeters, double _sizeZmeters,
                           double _cellSizeX, double _cellSizeY, double _cellSizeZ,
                           double _centerX, double _centerY, double _centerZ,
                           NDTCell *cellPrototype ) 
{

    sizeXmeters = _sizeXmeters;
    sizeYmeters = _sizeYmeters;
    sizeZmeters = _sizeZmeters;

    cellSizeX = _cellSizeX;
    cellSizeY = _cellSizeY;
    cellSizeZ = _cellSizeZ;

    sizeX = abs(ceil(sizeXmeters/cellSizeX));
    sizeY = abs(ceil(sizeYmeters/cellSizeY));
    sizeZ = abs(ceil(sizeZmeters/cellSizeZ));

    centerX = _centerX;
    centerY = _centerY;
    centerZ = _centerZ;

    protoType = cellPrototype->clone();
    initialize();


}

void LazyGrid::setCenter(const double &cx, const double &cy, const double &cz)
{
    centerX = cx;
    centerY = cy;
    centerZ = cz;

    centerIsSet = true;
    if(sizeIsSet)
    {
        initialize();
    }
}

void LazyGrid::setSize(const double &sx, const double &sy, const double &sz)
{

    sizeXmeters = sx;
    sizeYmeters = sy;
    sizeZmeters = sz;

    sizeX = abs(ceil(sizeXmeters/cellSizeX));
    sizeY = abs(ceil(sizeYmeters/cellSizeY));
    sizeZ = abs(ceil(sizeZmeters/cellSizeZ));

    sizeIsSet = true;
    if(centerIsSet)
    {
        initialize();
    }
}

void makeFolder(std::string path){
    if (mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
    {
        if( errno == EEXIST ) {
           // alredy exists
        } else {
           // something else
            std::cout << "cannot create map folder error:" << strerror(errno) << std::endl;
        }
    }
}


void LazyGrid::initialize()
{
    sizeX=ceil(sizeX/2)*2;
    sizeY=ceil(sizeY/2)*2;
    sizeZ=ceil(sizeZ/2)*2;

    clearCells();
    dataArray.initialize(sizeX,sizeY,sizeZ);
    initialized = true;
    makeFolder("/tmp/maps");
    finalpath="/tmp/maps/s_"+std::to_string(semantic_index);
    makeFolder(finalpath);
    finalpath=finalpath+"/r_"+std::to_string(cellSizeX);
    makeFolder(finalpath);
}
void LazyGrid::clearCells(){
        for(std::set<NDTCell*>::iterator it=activeCells.begin();it!=activeCells.end();++it)
            if(*it)
                delete *it;
        activeCells.clear();
        if(dataArray.array!=NULL)
            free( dataArray.array);
}

LazyGrid::~LazyGrid()
{
    if(initialized)
    {
        clearCells();
        if(protoType!=NULL)
        {
            delete protoType;
        }
        //fprintf(stderr,"Deleted %d cells and array of (%d x %d)!!!\n",cnt, sizeX, sizeY);
    }
}
void LazyGrid::addNDTCell(NDTCell* cell)
{
    pcl::PointXYZ cellCenter = cell->getCenter();
    if(inRange(cellCenter))
    {
        //std::cout<<".";
        int indX,indY,indZ;
        this->getIndexForPoint(cellCenter,indX,indY,indZ);
        if(dataArray(indX, indY, indZ)!=NULL)
        {
            deallocateCell(indX,indY,indZ);
            //std::cerr<<"DEALLOCED\n";
        }
        dataArray.set(indX,indY,indZ, cell);
        activeCells.insert(cell);
        return; //cell must not be deleted here.
    }
    delete cell;
}
void LazyGrid::loadCells(int index_start, int index_end)
{
    boost::filesystem::path p(finalpath);
    boost::filesystem::directory_iterator it{p};
//    std::cerr<<index_start<<" "<<index_end<<std::endl;
    while(it!=boost::filesystem::directory_iterator{})
    {
        std::string fname = it->path().filename().string();
//        std::cerr<<"FNAME: "<<fname<<std::endl;
        int t = stoi(fname.substr(0,fname.find('_')));
        if(t<=index_end&&t>=index_start)
        {
#ifdef BIN_ARCHIVE
            std::ifstream ifs(it->path().string(),std::ios::binary);
            //std::cerr<<it->path().string()<<std::endl;
            assert(ifs.good());

            //compression
            /*
            boost::iostreams::filtering_istreambuf fis;
            fis.push(boost::iostreams::zlib_decompressor());
            fis.push(ifs);
            */

            boost::archive::binary_iarchive ia(ifs);
            std::streampos aOffset= ifs.tellg();
            std::streampos strmEnd= ifs.seekg(0,std::ios_base::end).tellg();
            ifs.seekg(aOffset);
            while(ifs.tellg()<strmEnd)
#else
            std::ifstream ifs(it->path().string());
            assert(ifs.good());
            boost::archive::text_iarchive ia(ifs);
            while(ifs.good())
#endif
            {
                NDTCell* cell = protoType->clone();
                ia>>*cell;
                if(cell->cloud_index>=index_start && cell->cloud_index<=index_end )
                    addNDTCell(cell);
                else
                    delete cell;
            }
        }
        it++;
    }
//    std::cerr<<"D"<<std::endl;
}
void LazyGrid::deallocateCell(int i,int j, int k)
{
    SpatialIndex::CellVectorItr del_iter=activeCells.find(dataArray(i,j,k));
    if(del_iter!=activeCells.end())
        activeCells.erase(del_iter);
    delete dataArray(i,j,k);
    dataArray.set(i,j,k, NULL);
 
}
#ifdef BIN_ARCHIVE
void LazyGrid::dealocateCells(int dim_ind, int newP, boost::archive::binary_oarchive& oa, unsigned int &min_index)// newP only -1 or +1
#else
void LazyGrid::dealocateCells(int dim_ind, int newP, boost::archive::text_oarchive& oa, unsigned int &min_index)// newP only -1 or +1
#endif
{
    if(newP>0)
        newP=1;
    else if(newP<0)
        newP=-1;
    else
        return;
    int start_val[3]={0,0,0};
    int end_val[3]={0,0,0};
    int bH[3]={0,0,0};
    bH[dim_ind]=newP;
    for(int i=0;i<3;i++)
    {
        //end_val[i]=dataArray.size[i]/2 + sensor_pose[i];
        //start_val[i]=end_val[i] - dataArray.size[i];
        start_val[i]=0;
        end_val[i]=dataArray.size[i];
    }
    if(newP==1)
    {
        start_val[dim_ind]=dataArray.size[dim_ind]/2+sensor_pose[dim_ind];
        end_val[dim_ind]=start_val[dim_ind]+1;
    }
    if(newP==-1)
    {
        start_val[dim_ind]=dataArray.size[dim_ind]/2+sensor_pose[dim_ind]-1;
        end_val[dim_ind]=start_val[dim_ind]+1;
    }
    for(int i=start_val[0];i<end_val[0];i++)
        for(int j=start_val[1];j<end_val[1];j++)
            for(int k=start_val[2];k<end_val[2];k++)
            {
                if(dataArray(i,j,k)==NULL)
                    continue;
                if(true || dataArray(i,j,k)->hasGaussian_)
                {
                    if ( dataArray(i,j,k)->cloud_index < min_index)
                        min_index = dataArray(i,j,k)->cloud_index ;
                    oa << *dataArray(i,j,k);
                }
                deallocateCell(i,j,k);
            }
    sensor_pose[dim_ind]+=newP;
    return ;
}
bool fexists(const std::string& filename)
{
    std::ifstream ifile(filename.c_str());
    return (bool)ifile;
}

void LazyGrid::setSensorPose(const double *pose)
{
    cloud_index++;
    int indPose[3];
    translation<<pose[0],pose[1],pose[2];
    indPose[0] = floor(pose[0]/cellSizeX);
    indPose[1] = floor(pose[1]/cellSizeY);
    indPose[2] = floor(pose[2]/cellSizeZ);
    unsigned int min_index=UINT_MAX;
    std::string tmp_name  = finalpath+"/tmp_" + std::to_string(rand()%1000);
#ifdef BIN_ARCHIVE
    std::ofstream ofs(tmp_name,std::ios::binary);

    /*
    //Compression
    boost::iostreams::filtering_ostreambuf fosb;
    fosb.push(boost::iostreams::zlib_compressor(boost::iostreams::zlib::best_speed));
    fosb.push(ofs);
    */

    boost::archive::binary_oarchive oa(ofs);
#else
    std::ofstream ofs(tmp_name);
    boost::archive::text_oarchive oa(ofs);
#endif
    for(int i=0;i<3;i++)
        while(sensor_pose[i]!=indPose[i])
            dealocateCells(i, indPose[i]-sensor_pose[i], oa, min_index);
    ofs.close();
    if(min_index<UINT_MAX)
    {
        std::string fname=(finalpath+"/"+std::to_string(min_index))+"_";
        int i=0;
        while(fexists(fname+std::to_string(i)))
            i++;
        fname = fname+std::to_string(i);
        std::rename(tmp_name.c_str(), fname.c_str());
    }
    else
        std::remove(tmp_name.c_str());
    //std::cerr<<"POSE SET"<<std::endl;
    //loadCells(0, cloud_index);

}
bool LazyGrid::isValid(const pcl::PointXYZ &p, NDTCell* cell)
{
    auto c = cell->getCenter();
    if(abs(c.x-p.x)>cellSizeX ||abs(c.y-p.y)>cellSizeY ||abs(c.z-p.z)>cellSizeZ)
        return false;
    return true;
}
NDTCell* LazyGrid::getCellForPoint(const pcl::PointXYZ &point)
{

    int indX,indY,indZ;
    this->getIndexForPoint(point,indX,indY,indZ);

    if(!initialized) return NULL;
    if(dataArray.array==NULL) return NULL;
    //    cout<<"LZ: "<<indX<<" "<<indY<<" "<<indZ<<endl;
    auto rval = dataArray(indX, indY, indZ);
    if(!isValid(point, rval))
        return NULL;
    std::cout<< rval->cloud_index<<std::endl;
    return rval;
}


bool LazyGrid::getCellAtAllocate(const pcl::PointXYZ &pt, NDTCell* &cell)//RETURNS TRUE IF CELL ALLOCATED
{
    bool allocated=false;
    if(!inRange(pt))
        return allocated;
    cell = NULL;
    pcl::PointXYZ point = pt;
    if(std::isnan(point.x) ||std::isnan(point.y) ||std::isnan(point.z))
    {
        return allocated;
    }
    int indX,indY,indZ;
    this->getIndexForPoint(point,indX,indY,indZ);
    pcl::PointXYZ centerCell;
    centerCell.x = floor(pt.x/cellSizeX)*cellSizeX + cellSizeX/2;
    centerCell.y = floor(pt.y/cellSizeY)*cellSizeY + cellSizeY/2;
    centerCell.z = floor(pt.z/cellSizeZ)*cellSizeZ + cellSizeZ/2;

    //if(indX >= sizeX || indY >= sizeY || indZ >= sizeZ || indX<0 || indY<0 || indZ<0)

    if(dataArray.array==NULL) return allocated;
    if(!initialized) return allocated;

    if(dataArray(indX,indY,indZ)==NULL)
    {
        allocated=true;
        //initialize cell
        dataArray.set(indX,indY,indZ, protoType->clone());
        activeCells.insert(dataArray(indX,indY,indZ));
        dataArray(indX,indY,indZ)->setDimensions(cellSizeX,cellSizeY,cellSizeZ);

        dataArray(indX,indY,indZ)->setCenter(centerCell);
        dataArray(indX,indY,indZ)->cloud_index = cloud_index;
        /*
           cout<<"center: "<<centerX<<" "<<centerY<<" "<<centerZ<<endl;
           cout<<"size  : "<<sizeX<<" "<<sizeY<<" "<<sizeZ<<endl;
           cout<<"p  : "<<point.x<<" "<<point.y<<" "<<point.z<<endl;
           cout<<"c  : "<<centerCell.x<<" "<<centerCell.y<<" "<<centerCell.z<<endl;
           cout<<"id : "<<indX<<" "<<indY<<" "<<indZ<<endl;
           cout<<"cs : "<<cellSizeX<<" "<<cellSizeY<<" "<<cellSizeZ<<endl;
         */
    }
    else if (!isValid(point, dataArray(indX,indY,indZ)))
    {
        std::cout<<pt<<std::endl;
        std::cout<<dataArray(indX,indY,indZ)->getCenter()<<std::endl;
        assert(false);
        return allocated;
    }
    cell = dataArray(indX,indY,indZ);
    return allocated;
}


NDTCell* LazyGrid::addPoint(const pcl::PointXYZ &point_c)
{
  NDTCell* cell=NULL;
  this->getCellAtAllocate(point_c, cell);
  if (cell != NULL)
    cell->addPoint(point_c);
  return cell;
}

typename SpatialIndex::CellVectorItr LazyGrid::begin()
{
    return activeCells.begin();
}

typename SpatialIndex::CellVectorConstItr LazyGrid::begin() const
{
    return activeCells.begin();
}

typename SpatialIndex::CellVectorItr LazyGrid::end()
{
    return activeCells.end();
}

typename SpatialIndex::CellVectorConstItr LazyGrid::end() const
{
    return activeCells.end();
}

int LazyGrid::size()
{
    return activeCells.size();
}

SpatialIndex* LazyGrid::clone() const
{
    LazyGrid * rt = new LazyGrid(cellSizeX);
    rt->semantic_index = semantic_index;
    return dynamic_cast<SpatialIndex *>(rt);
}

SpatialIndex* LazyGrid::copy() const
{
    LazyGrid *ret = new LazyGrid(cellSizeX);
    typename std::vector<NDTCell*>::const_iterator it;
    for(std::set<NDTCell*>::iterator it=activeCells.begin();it!=activeCells.end();++it)
    {
        NDTCell* r = (*it);
        if(r == NULL) continue;
        for(unsigned int i=0; i<r->points_.size(); i++)
        {
            ret->addPoint(r->points_[i]);
        }
    }
    return ret;
}

void LazyGrid::getNeighbors(const pcl::PointXYZ &point, const double &radius, std::vector<NDTCell*> &cells)
{
    ///NOT TESTED, not modified
    std::cerr<<"OH NO! SHOULDNT BE HERE! LazyGrid:N4578"<<std::endl;
    int indX,indY,indZ;
    this->getIndexForPoint(point, indX,indY,indZ);
    if(indX >= sizeX || indY >= sizeY || indZ >= sizeZ)
    {
        cells.clear();
        return;
    }

    for(int x = indX - radius/cellSizeX; x<=indX+radius/cellSizeX; x++)
    {
        if(x < 0 || x >= sizeX) continue;
        for(int y = indY - radius/cellSizeY; y<=indY+radius/cellSizeY; y++)
        {
            if(y < 0 || y >= sizeY) continue;
            for(int z = indZ - radius/cellSizeZ; z<=indZ+radius/cellSizeZ; z++)
            {
                if(z < 0 || z >= sizeZ) continue;
                if(dataArray(x,y,z)==NULL) continue;
                cells.push_back(dataArray(x,y,z));
            }
        }
    }

}

void LazyGrid::getNeighborsShared(const pcl::PointXYZ &point, const double &radius, std::vector<boost::shared_ptr< NDTCell > > &cells)
{
    ///NOT TESTED, not modified
    std::cerr<<"OH NO! SHOULDNT BE HERE! LazyGrid:N4579"<<std::endl;
    int indX,indY,indZ;
    this->getIndexForPoint(point, indX,indY,indZ);
    if(indX >= sizeX || indY >= sizeY || indZ >= sizeZ)
    {
        cells.clear();
        return;
    }

    for(int x = indX - radius/cellSizeX; x<=indX+radius/cellSizeX; x++)
    {
        if(x < 0 || x >= sizeX) continue;
        for(int y = indY - radius/cellSizeY; y<=indY+radius/cellSizeY; y++)
        {
            if(y < 0 || y >= sizeY) continue;
            for(int z = indZ - radius/cellSizeZ; z<=indZ+radius/cellSizeZ; z++)
            {
                if(z < 0 || z >= sizeZ) continue;
                if(dataArray(x,y,z)==NULL) continue;
				NDTCell* nd = dataArray(x,y,z)->copy();
				boost::shared_ptr< NDTCell > smart_pointer(nd);
				cells.push_back(smart_pointer);
            }
        }
    }

}
bool LazyGrid::inRange(const pcl::PointXYZ& p)
{
    if( p.x < (sensor_pose[0]+sizeX/2)*cellSizeX &&
        p.x > (sensor_pose[0]-sizeX/2)*cellSizeX &&
        p.y < (sensor_pose[1]+sizeY/2)*cellSizeY &&
        p.y > (sensor_pose[1]-sizeY/2)*cellSizeY &&
        p.z < (sensor_pose[2]+sizeZ/2)*cellSizeZ &&
        p.z > (sensor_pose[2]-sizeZ/2)*cellSizeZ )
        return true;
    else
        return false;
}

void LazyGrid::getIndexForPoint(const pcl::PointXYZ& point, int &indX, int &indY, int &indZ)
{
    indX = floor(point.x/cellSizeX);
    indY = floor(point.y/cellSizeY);
    indZ = floor(point.z/cellSizeZ);
}

std::vector<NDTCell*> LazyGrid::getClosestCells(const pcl::PointXYZ &pt)
{
    int indXn,indYn,indZn;
    int indX,indY,indZ;
    this->getIndexForPoint(pt,indX,indY,indZ);
    std::vector<NDTCell*> cells;

    int i =2; //how many cells on each side
    //the strange thing for the indeces is for convenience of writing
    //basicly, we go through 2* the number of cells and use parity to
    //decide if we subtract or add. should work nicely
    for(int x=1; x<2*i+2; x++)
    {
        indXn = (x%2 == 0) ? indX+x/2 : indX-x/2;
        for(int y=1; y<2*i+2; y++)
        {
            indYn = (y%2 == 0) ? indY+y/2 : indY-y/2;
            for(int z=1; z<2*i+2; z++)
            {
                indZn = (z%2 == 0) ? indZ+z/2 : indZ-z/2;
                if(checkCellforNDT(indXn,indYn,indZn,true))
                {
                    cells.push_back(dataArray(indXn,indYn,indZn));
                }
            }
        }
    }
    return cells;
}

std::vector< NDTCell* > LazyGrid::getClosestNDTCells(const pcl::PointXYZ &point, int &n_neigh, bool checkForGaussian)
{
    std::vector<NDTCell*> cells;
    if(!inRange(point))
        return cells;

    int indXn,indYn,indZn;
    int indX,indY,indZ;
    this->getIndexForPoint(point,indX,indY,indZ);

    int i = n_neigh; //how many cells on each side

    //for(int i=1; i<= maxNumberOfCells; i++) {
    //the strange thing for the indeces is for convenience of writing
    //basicly, we go through 2* the number of cells and use parity to
    //decide if we subtract or add. should work nicely
    for(int x=1; x<2*i+2; x++)
    {
        indXn = (x%2 == 0) ? indX+x/2 : indX-x/2;
        for(int y=1; y<2*i+2; y++)
        {
            indYn = (y%2 == 0) ? indY+y/2 : indY-y/2;
            for(int z=1; z<2*i+2; z++)
            {
                indZn = (z%2 == 0) ? indZ+z/2 : indZ-z/2;
                if(checkCellforNDT(indXn,indYn,indZn,checkForGaussian))
                {
                    if(pcl::geometry::distance(point, dataArray(indXn,indYn,indZn)->getCenter())<(1+n_neigh)*cellSizeX)
                        cells.push_back(dataArray(indXn,indYn,indZn));
                }
            }
        }
    }

    return cells;
}

std::vector<boost::shared_ptr< NDTCell > > LazyGrid::getClosestNDTCellsShared(const pcl::PointXYZ &point, int &n_neigh, bool checkForGaussian)
{

    int indXn,indYn,indZn;
    int indX,indY,indZ;
    this->getIndexForPoint(point,indX,indY,indZ);
    std::vector<boost::shared_ptr< NDTCell > > cells;

    int i = n_neigh; //how many cells on each side

    //for(int i=1; i<= maxNumberOfCells; i++) {
    //the strange thing for the indeces is for convenience of writing
    //basicly, we go through 2* the number of cells and use parity to
    //decide if we subtract or add. should work nicely
    for(int x=1; x<2*i+2; x++)
    {
        indXn = (x%2 == 0) ? indX+x/2 : indX-x/2;
        for(int y=1; y<2*i+2; y++)
        {
            indYn = (y%2 == 0) ? indY+y/2 : indY-y/2;
            for(int z=1; z<2*i+2; z++)
            {
                indZn = (z%2 == 0) ? indZ+z/2 : indZ-z/2;
                if(checkCellforNDT(indXn,indYn,indZn,checkForGaussian))
                    if(pcl::geometry::distance(point, dataArray(indXn,indYn,indZn)->getCenter())<(1+n_neigh)*cellSizeX)
                    {
                        NDTCell* nd = dataArray(indXn,indYn,indZn)->copy();
                        boost::shared_ptr< NDTCell > smart_pointer(nd);
                        cells.push_back(smart_pointer);
                    }
            }
        }
    }

    return cells;
}

NDTCell* LazyGrid::getClosestNDTCell(const pcl::PointXYZ &point, bool checkForGaussian)
{
    int indXn,indYn,indZn;
    int indX,indY,indZ;
    this->getIndexForPoint(point,indX,indY,indZ);
    NDTCell *ret = NULL;
    std::vector<NDTCell*> cells;

    if(!checkForGaussian)
    {
        //just give me whatever is in this cell
        if(checkCellforNDT(indX,indY,indZ,checkForGaussian))
            if(pcl::geometry::distance(point, dataArray(indX,indY,indZ)->getCenter())<(2)*cellSizeX)
            {
                ret = (dataArray(indX,indY,indZ));
            }
        return ret;
    }

    int i =1; //how many cells on each side

    //for(int i=1; i<= maxNumberOfCells; i++) {
    //the strange thing for the indeces is for convenience of writing
    //basicly, we go through 2* the number of cells and use parity to
    //decide if we subtract or add. should work nicely
    for(int x=1; x<2*i+2; x++)
    {
        indXn = (x%2 == 0) ? indX+x/2 : indX-x/2;
        for(int y=1; y<2*i+2; y++)
        {
            indYn = (y%2 == 0) ? indY+y/2 : indY-y/2;
            for(int z=1; z<2*i+2; z++)
            {
                indZn = (z%2 == 0) ? indZ+z/2 : indZ-z/2;
                if(checkCellforNDT(indXn,indYn,indZn))
                    if(pcl::geometry::distance(point, dataArray(indXn,indYn,indZn)->getCenter())<(2*i+2)*cellSizeX)
                    {
                        ret = (dataArray(indXn,indYn,indZn));
                        cells.push_back(ret);
                    }
            }
        }
    }

    double minDist = INT_MAX;
    Eigen::Vector3d tmean;
    pcl::PointXYZ pt = point;
    for(unsigned int i=0; i<cells.size(); i++)
    {
        tmean = cells[i]->getMean();
        tmean(0) -= pt.x;
        tmean(1) -= pt.y;
        tmean(2) -= pt.z;
        double d = tmean.norm();
        if(d<minDist)
        {
            minDist = d;
            ret = cells[i];
        }
    }
    cells.clear();
    return ret;
}

bool LazyGrid::checkCellforNDT(int indX, int indY, int indZ, bool checkForGaussian)
{

    if(dataArray(indX,indY,indZ)!=NULL)
    {
	    if( dataArray(indX,indY,indZ)->hasGaussian_ || (!checkForGaussian))
            return true;
    }
    return false;
}

void LazyGrid::setCellType(NDTCell *type)
{
    if(type!=NULL)
    {
      if (protoType != NULL) {
        delete protoType;
      }
        protoType = type->clone();
    }
}

void LazyGrid::getCellSize(double &cx, double &cy, double &cz)
{
    cx = cellSizeX;
    cy = cellSizeY;
    cz = cellSizeZ;
}

void LazyGrid::getCenter(double &cx, double &cy, double &cz)
{
    cx = centerX;
    cy = centerY;
    cz = centerZ;

}

void LazyGrid::getGridSize(int &cx, int &cy, int &cz)
{
    cx = sizeX;
    cy = sizeY;
    cz = sizeZ;
}

void LazyGrid::getGridSizeInMeters(double &cx, double &cy, double &cz)
{
    cx = sizeXmeters;
    cy = sizeYmeters;
    cz = sizeZmeters;
}

int LazyGrid::loadFromJFF(FILE * jffin)
{
    if(initialized)
    {
        for(std::set<NDTCell*>::iterator it=activeCells.begin();it!=activeCells.end();++it)
            if(*it)
                delete *it;
		activeCells.clear();
		free( dataArray.array);
		if(protoType!=NULL)
		{
			delete protoType;
		}
		//fprintf(stderr,"Deleted %d cells and array of (%d x %d)!!!\n",cnt, sizeX, sizeY);
    }
    double lazyGridData[9]; // = { sizeXmeters, sizeYmeters, sizeZmeters,
    //     cellSizeX,   cellSizeY,   cellSizeZ,
    //     centerX,     centerY,     centerZ };
    NDTCell prototype_;
    if(fread(&lazyGridData, sizeof(double), 9, jffin) <= 0)
    {
        JFFERR("reading lazyGridData failed");
    }
    if(fread(&prototype_, sizeof(NDTCell), 1, jffin) <= 0)
    {
        JFFERR("reading prototype_ failed");
    }

    // just in case someone was messing around with the new NDTMap
    centerIsSet = false;
    sizeIsSet = false;

    protoType = prototype_.clone();

    //std::cerr<<"size meters: "<<lazyGridData[0]<<" "<<lazyGridData[1]<<" "<<lazyGridData[2]<<std::endl;
    //std::cerr<<"cell size: "<<lazyGridData[3]<<" "<<lazyGridData[4]<<" "<<lazyGridData[5]<<std::endl;
    //std::cerr<<"center meters: "<<lazyGridData[6]<<" "<<lazyGridData[7]<<" "<<lazyGridData[8]<<std::endl;

    this->setSize(lazyGridData[0], lazyGridData[1], lazyGridData[2]);

    cellSizeX = lazyGridData[3];
    cellSizeY = lazyGridData[4];
    cellSizeZ = lazyGridData[5];

    this->setCenter(lazyGridData[6], lazyGridData[7], lazyGridData[8]);

    int indX, indY, indZ;
    pcl::PointXYZ centerCell;

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

        if(feof(jffin)) break;
		centerCell = prototype_.getCenter();
        this->getIndexForPoint(centerCell, indX, indY, indZ);

        if(indX < 0 || indX >= sizeX) continue; 
        if(indY < 0 || indY >= sizeY) continue; 
        if(indZ < 0 || indZ >= sizeZ) continue; 
		if(!initialized) return -1;
        if(dataArray(indX,indY,indZ) != NULL)
			delete dataArray(indX,indY,indZ);

	    dataArray.set(indX,indY,indZ, prototype_.copy());
	    activeCells.insert(dataArray(indX,indY,indZ));
    }

    return 0;
}

bool LazyGrid::traceLine(const Eigen::Vector3d &origin, const pcl::PointXYZ &endpoint,const Eigen::Vector3d &diff_ ,
	const double& maxz, std::vector<NDTCell*> &cells)
{
	if(endpoint.z>maxz)
	{
		return false;
	}
	
	double min1 = std::min(cellSizeX,cellSizeY);
	double min2 = std::min(cellSizeZ,cellSizeY);
	double resolution = std::min(min1,min2); ///Select the smallest resolution
	
	if(resolution<0.01)
	{
		fprintf(stderr,"Resolution very very small (%lf) :( \n",resolution);
		return false;
	}
	double l = diff_.norm();
	int N = l / (resolution);
	//if(N <= 0)
	//{
	//	//fprintf(stderr,"N=%d (r=%lf l=%lf) :( ",N,resolution,l);
	//	return false;
	//}
	
	NDTCell* ptCell = NULL;    
	pcl::PointXYZ pt;
	pcl::PointXYZ po;
	po.x = origin(0); po.y = origin(1); po.z = origin(2);
	Eigen::Vector3d diff = diff_/(float)N;
	
	int idxo=0, idyo=0,idzo=0;
	for(int i=0; i<N-2; i++)
	{
		pt.x = origin(0) + ((float)(i+1)) *diff(0);
		pt.y = origin(1) + ((float)(i+1)) *diff(1);
		pt.z = origin(2) + ((float)(i+1)) *diff(2);
		int idx,idy,idz;
        this->getIndexForPoint(pt, idx, idy, idz);
		///We only want to check every cell once, so
		///increase the index if we are still in the same cell
		if(idx == idxo && idy==idyo && idz ==idzo)
		{
			continue;
		}
		else
		{
			idxo = idx;
			idyo = idy;
			idzo = idz;
		}
        bool cell_allocated = this->getCellAtAllocate(pt, ptCell);
        if (ptCell != NULL)
        {
            if(cell_allocated)
                ptCell->addPoint(pt);
            cells.push_back(ptCell);
        }
	}
	return true;
	
}


bool LazyGrid::insertCell(NDTCell cell){
    pcl::PointXYZ centerCell;
    int indX,indY,indZ;
    float r,g,b;
    double xs,ys,zs;
    float occ;
    unsigned int N;
    centerCell = cell.getCenter();
    this->getIndexForPoint(centerCell, indX, indY, indZ);
    if(!initialized) return false;
    if(dataArray.array == NULL) return false;

    if(dataArray(indX,indY,indZ) != NULL)
    {
      NDTCell* ret = dataArray(indX,indY,indZ);
      cell.getRGB(r,g,b);
      cell.getDimensions(xs,ys,zs);

      ret->setDimensions(xs,ys,zs);
      ret->setCenter(centerCell);
      ret->setMean(cell.getMean());
      ret->setCov(cell.getCov());
      ret->setRGB(r,g,b);
      ret->setOccupancy(cell.getOccupancy());
      ret->setEmptyval(cell.getEmptyval());
      ret->setEventData(cell.getEventData());
      ret->setN(cell.getN());
      ret->isEmpty = cell.isEmpty;
      ret->hasGaussian_ = cell.hasGaussian_;
      ret->consistency_score = cell.consistency_score;

    } else {
      //initialize cell
     // std::cerr<<"NEW CELL\n";
      dataArray.set(indX,indY,indZ,cell.copy());
      activeCells.insert(dataArray(indX,indY,indZ));
    }
    return true;
    }
void LazyGrid::InitializeDefaultValues(){

  protoType  =  NULL;
  activeCells.clear();
  centerIsSet = sizeIsSet   = initialized = false;
  sizeXmeters = sizeYmeters = sizeZmeters=0;
  cellSizeX   = cellSizeY   = cellSizeZ=0;
  centerX     = centerY     = centerZ =0;
  sizeX       = sizeY       = sizeZ=0;
}
std::string LazyGrid::GetDataString()
{
  std::stringstream ss;
  //  for(int i=0;i<activeCells.size();i++){
  //    if(activeCells[i]->hasGaussian_)
  //      ss<<activeCells[i]->ToString()<<std::endl;
  //  }
  for(int i=0; i<sizeX; i++)
  {
    for(int j=0; j<sizeY; j++)
    {
      for(int k=0; k<sizeZ; k++)
      {
        if(dataArray(i,j,k)!=NULL){
          if(dataArray(i,j,k)->hasGaussian_){
            ss<<"mean=\n"<< dataArray(i,j,k)->getMean()<<"\ncov=\n"<<dataArray(i,j,k)->getCov()<<std::endl;
          }
        }


      }
    }
  }
  return ss.str();
}
std::string LazyGrid::ToString(){
  std::stringstream ss;
  ss <<"<<LazyGrid:initialized="<<initialized<<", active cells="<<activeCells.size()<<"\nCenterIsSet="<<centerIsSet<< "\nSizeIsSet="<<sizeIsSet<<"\nsizeXmeters="<<sizeXmeters<<", sizeYmeters="<<sizeYmeters<<", sizeZmeters="<<sizeZmeters<<"\ncellSizeX="<<cellSizeX<<", cellSizeY="<<cellSizeY<<", cellSizeZ="<<cellSizeZ<<std::endl;
  ss<<"centerX="<<centerX<<", centerY="<<centerY<<", centerZ="<<centerZ<<"\nsizeX="<<sizeX<<", sizeY="<<sizeY<<", sizeZ="<<sizeZ<<std::endl;
  ss<<"prototype="<<protoType->ToString()<<"lazygrid>>"<<std::endl;
  ss<<"datastring=\n"<<GetDataString()<<std::endl;
  return ss.str();
}
} //end namespace
