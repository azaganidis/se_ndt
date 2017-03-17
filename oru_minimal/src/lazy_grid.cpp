#include <cstring>
#include <cstdio>
#include <ndt_map/lazy_grid.h>

#define JFFERR(x) std::cerr << x << std::endl; return -1;

namespace lslgeneric
{

LazyGrid::LazyGrid(double cellSize) 
{
    initialized = false;
    centerIsSet = false;
    sizeIsSet = false;
    cellSizeX = cellSizeY = cellSizeZ = cellSize;

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

void LazyGrid::initializeAll()
{
    if(!initialized)
    {
        this->initialize();
    }
    pcl::PointXYZ centerCell;
    for(int i=0; i<sizeX; i++)
    {
        for(int j=0; j<sizeY; j++)
        {
            for(int k=0; k<sizeZ; k++)
            {
                dataArray[i][j][k] = new NDTCell();
                dataArray[i][j][k]->setDimensions(cellSizeX,cellSizeY,cellSizeZ);

                int idcX, idcY, idcZ;
                pcl::PointXYZ center;
                center.x = centerX;
                center.y = centerY;
                center.z = centerZ;
                this->getIndexForPoint(center, idcX,idcY,idcZ);
                centerCell.x = centerX + (i-idcX)*cellSizeX;
                centerCell.y = centerY + (j-idcY)*cellSizeY;
                centerCell.z = centerZ + (k-idcZ)*cellSizeZ;
                dataArray[i][j][k]->setCenter(centerCell);
                activeCells.push_back(dataArray[i][j][k]);
            }
        }
    }
}

void LazyGrid::initialize()
{

    dataArray = new NDTCell***[sizeX];
    for(int i=0; i<sizeX; i++)
    {
        dataArray[i] = new NDTCell**[sizeY];
        for(int j=0; j<sizeY; j++)
        {
            dataArray[i][j] = new NDTCell*[sizeZ];
            //set all cells to NULL
            memset(dataArray[i][j],0,sizeZ*sizeof(NDTCell*));
        }
    }
    initialized = true;
}

LazyGrid::~LazyGrid()
{
    if(initialized)
    {
	//fprintf(stderr,"LAZY GRID DESTRUCTION -- ");
	int cnt = 0;
	//go through all cells and delete the non-NULL ones
	for(unsigned int i=0; i<activeCells.size(); ++i)
	{
	    if(activeCells[i])
	    {
		delete activeCells[i];
		cnt++;
	    }
	}
	for(int i=0; i<sizeX; i++)
	{
	    for(int j=0; j<sizeY; j++)
	    {
		delete[] dataArray[i][j];
	    }
	    delete[] dataArray[i];
	}
	delete[] dataArray;
	if(protoType!=NULL)
	{
	    delete protoType;
	}
	//fprintf(stderr,"Deleted %d cells and array of (%d x %d)!!!\n",cnt, sizeX, sizeY);
    }
}

NDTCell* LazyGrid::getCellForPoint(const pcl::PointXYZ &point)
{

    int indX,indY,indZ;
    this->getIndexForPoint(point,indX,indY,indZ);

    if(indX >= sizeX || indY >= sizeY || indZ >= sizeZ || indX<0 || indY<0 || indZ<0) return NULL;
    if(!initialized) return NULL;
    if(dataArray == NULL) return NULL;
    if(dataArray[indX] == NULL) return NULL;
    if(dataArray[indX][indY] == NULL) return NULL;

    //    cout<<"LZ: "<<indX<<" "<<indY<<" "<<indZ<<endl;
    return dataArray[indX][indY][indZ];
}

NDTCell* LazyGrid::addPoint(const pcl::PointXYZ &point_c)
{

    pcl::PointXYZ point = point_c;
    if(std::isnan(point.x) ||std::isnan(point.y) ||std::isnan(point.z))
    {
        return NULL;
    }
    int indX,indY,indZ;
    this->getIndexForPoint(point,indX,indY,indZ);
    pcl::PointXYZ centerCell;

    if(indX >= sizeX || indY >= sizeY || indZ >= sizeZ || indX<0 || indY<0 || indZ<0)
    {
        return NULL;
    }

    if(!initialized) return NULL;
    if(dataArray == NULL) return NULL;
    if(dataArray[indX] == NULL) return NULL;
    if(dataArray[indX][indY] == NULL) return NULL;

    if(dataArray[indX][indY][indZ]==NULL)
    {
        //initialize cell
        dataArray[indX][indY][indZ] = protoType->clone();
        dataArray[indX][indY][indZ]->setDimensions(cellSizeX,cellSizeY,cellSizeZ);

        int idcX, idcY, idcZ;
        pcl::PointXYZ center;
        center.x = centerX;
        center.y = centerY;
        center.z = centerZ;
        this->getIndexForPoint(center, idcX,idcY,idcZ);
        centerCell.x = centerX + (indX-idcX)*cellSizeX;
        centerCell.y = centerY + (indY-idcY)*cellSizeY;
        centerCell.z = centerZ + (indZ-idcZ)*cellSizeZ;
        dataArray[indX][indY][indZ]->setCenter(centerCell);
        /*
           cout<<"center: "<<centerX<<" "<<centerY<<" "<<centerZ<<endl;
           cout<<"size  : "<<sizeX<<" "<<sizeY<<" "<<sizeZ<<endl;
           cout<<"p  : "<<point.x<<" "<<point.y<<" "<<point.z<<endl;
           cout<<"c  : "<<centerCell.x<<" "<<centerCell.y<<" "<<centerCell.z<<endl;
           cout<<"id : "<<indX<<" "<<indY<<" "<<indZ<<endl;
           cout<<"cs : "<<cellSizeX<<" "<<cellSizeY<<" "<<cellSizeZ<<endl;
         */
        activeCells.push_back(dataArray[indX][indY][indZ]);
    }
    dataArray[indX][indY][indZ]->addPoint(point);
    return dataArray[indX][indY][indZ];
}

typename SpatialIndex::CellVectorItr LazyGrid::begin()
{
    return activeCells.begin();
}

typename SpatialIndex::CellVectorItr LazyGrid::end()
{
    return activeCells.end();
}

int LazyGrid::size()
{
    return activeCells.size();
}

SpatialIndex* LazyGrid::clone() const
{
    return new LazyGrid(cellSizeX);
}

SpatialIndex* LazyGrid::copy() const
{
    LazyGrid *ret = new LazyGrid(cellSizeX);
    typename std::vector<NDTCell*>::const_iterator it;
    it = activeCells.begin();
    while(it!=activeCells.end())
    {
        NDTCell* r = (*it);
        if(r == NULL) continue;
        for(unsigned int i=0; i<r->points_.size(); i++)
        {
            ret->addPoint(r->points_[i]);
        }
        it++;
    }
    return ret;
}

void LazyGrid::getNeighbors(const pcl::PointXYZ &point, const double &radius, std::vector<NDTCell*> &cells)
{
    int indX,indY,indZ;
    this->getIndexForPoint(point, indX,indY,indZ);
    if(indX >= sizeX || indY >= sizeY || indZ >= sizeZ)
    {
        cells.clear();
        return;
    }

    for(int x = indX - radius/cellSizeX; x<indX+radius/cellSizeX; x++)
    {
        if(x < 0 || x >= sizeX) continue;
        for(int y = indY - radius/cellSizeY; y<indY+radius/cellSizeY; y++)
        {
            if(y < 0 || y >= sizeY) continue;
            for(int z = indZ - radius/cellSizeZ; z<indZ+radius/cellSizeZ; z++)
            {
                if(z < 0 || z >= sizeZ) continue;
                if(dataArray[x][y][z]==NULL) continue;
                cells.push_back(dataArray[x][y][z]);
            }
        }
    }

}

void LazyGrid::getIndexForPoint(const pcl::PointXYZ& point, int &indX, int &indY, int &indZ)
{
    indX = floor((point.x - centerX)/cellSizeX+0.5) + sizeX/2.0;
    indY = floor((point.y - centerY)/cellSizeY+0.5) + sizeY/2.0;
    indZ = floor((point.z - centerZ)/cellSizeZ+0.5) + sizeZ/2.0;
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
                    cells.push_back(dataArray[indXn][indYn][indZn]);
                }
            }
        }
    }
    return cells;
}

std::vector<NDTCell*> LazyGrid::getClosestNDTCells(const pcl::PointXYZ &point, int &n_neigh, bool checkForGaussian)
{

    int indXn,indYn,indZn;
    int indX,indY,indZ;
    this->getIndexForPoint(point,indX,indY,indZ);
    std::vector<NDTCell*> cells;

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
                    cells.push_back(dataArray[indXn][indYn][indZn]);
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
        {
            ret = (dataArray[indX][indY][indZ]);
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
                {
                    ret = (dataArray[indXn][indYn][indZn]);
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

    if(indX < sizeX && indY < sizeY && indZ < sizeZ &&
            indX >=0 && indY >=0 && indZ >=0)
    {
        if(dataArray[indX][indY][indZ]!=NULL)
        {
	    if(dataArray[indX][indY][indZ]->hasGaussian_ || (!checkForGaussian))
	    {
		return true;
	    }
        }
    }
    return false;
}

void LazyGrid::setCellType(NDTCell *type)
{
    if(type!=NULL)
    {
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

    std::cerr<<"size meters: "<<lazyGridData[0]<<" "<<lazyGridData[1]<<" "<<lazyGridData[2]<<std::endl;
    std::cerr<<"cell size: "<<lazyGridData[3]<<" "<<lazyGridData[4]<<" "<<lazyGridData[5]<<std::endl;
    std::cerr<<"center meters: "<<lazyGridData[6]<<" "<<lazyGridData[7]<<" "<<lazyGridData[8]<<std::endl;

    this->setSize(lazyGridData[0], lazyGridData[1], lazyGridData[2]);

    cellSizeX = lazyGridData[3];
    cellSizeY = lazyGridData[4];
    cellSizeZ = lazyGridData[5];

    this->setCenter(lazyGridData[6], lazyGridData[7], lazyGridData[8]);

    this->initializeAll();
    int indX, indY, indZ;
    float r,g,b;
    double xs,ys,zs;
    pcl::PointXYZ centerCell;
    float occ;
    unsigned int N;

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
	centerCell = prototype_.getCenter();
        this->getIndexForPoint(centerCell, indX, indY, indZ);
        if(indX < 0 || indX >= sizeX) continue; 
        if(indY < 0 || indY >= sizeY) continue; 
        if(indZ < 0 || indZ >= sizeZ) continue; 
	if(!initialized) return -1;
        if(dataArray == NULL) return -1;
        if(dataArray[indX] == NULL) return -1;
        if(dataArray[indX][indY] == NULL) return -1;

        if(dataArray[indX][indY][indZ] != NULL)
        {
            NDTCell* ret = dataArray[indX][indY][indZ];
	    prototype_.getRGB(r,g,b);
	    prototype_.getDimensions(xs,ys,zs);

	    ret->setDimensions(xs,ys,zs);
	    ret->setCenter(centerCell);
	    ret->setMean(prototype_.getMean());
	    ret->setCov(prototype_.getCov());
	    ret->setRGB(r,g,b);
	    ret->setOccupancy(prototype_.getOccupancy());
	    ret->setEmptyval(prototype_.getEmptyval());
	    ret->setEventData(prototype_.getEventData());
	    ret->setN(prototype_.getN());
	    ret->isEmpty = prototype_.isEmpty;
	    ret->hasGaussian_ = prototype_.hasGaussian_;
	    ret->consistency_score = prototype_.consistency_score;

	} else {
	    //initialize cell
	    std::cerr<<"NEW CELL\n";
	    dataArray[indX][indY][indZ] = prototype_.copy();
	    activeCells.push_back(dataArray[indX][indY][indZ]);
	}
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
		idx = floor((pt.x - centerX)/cellSizeX+0.5) + sizeX/2.0;
		idy = floor((pt.y - centerY)/cellSizeY+0.5) + sizeY/2.0;
		idz = floor((pt.z - centerZ)/cellSizeZ+0.5) + sizeZ/2.0;
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
		
		if(idx < sizeX && idy < sizeY && idz < sizeZ && idx >=0 && idy >=0 && idz >=0){
			ptCell = dataArray[idx][idy][idz];
			if(ptCell !=NULL) {
				cells.push_back(ptCell);
			} else {
				this->addPoint(pt); ///Add fake point to initialize!
			}
		}
	}
	return true;
	
}

bool LazyGrid::traceLineWithEndpoint(const Eigen::Vector3d &origin, const pcl::PointXYZ &endpoint,const Eigen::Vector3d &diff_ ,const double& maxz, std::vector<NDTCell*> &cells, Eigen::Vector3d &final_point)
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
    NDTCell* ptCell = NULL;    
    pcl::PointXYZ pt;
    pcl::PointXYZ po;
    po.x = origin(0); po.y = origin(1); po.z = origin(2);
    if(N == 0)
    {
	//fprintf(stderr,"N=%d (r=%lf l=%lf) :( ",N,resolution,l);
	//return false;
	this->getNDTCellAt(po,ptCell);
	if(ptCell!=NULL) {
	    cells.push_back(ptCell);
	}
	return true;
    }

    Eigen::Vector3d diff = diff_/(float)N;

    int idxo=0, idyo=0,idzo=0;
    bool complete = true;
    for(int i=0; i<N-2; i++)
    {
	pt.x = origin(0) + ((float)(i+1)) *diff(0);
	pt.y = origin(1) + ((float)(i+1)) *diff(1);
	pt.z = origin(2) + ((float)(i+1)) *diff(2);
	int idx,idy,idz;
	idx = floor((pt.x - centerX)/cellSizeX+0.5) + sizeX/2.0;
	idy = floor((pt.y - centerY)/cellSizeY+0.5) + sizeY/2.0;
	idz = floor((pt.z - centerZ)/cellSizeZ+0.5) + sizeZ/2.0;
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

	if(idx < sizeX && idy < sizeY && idz < sizeZ && idx >=0 && idy >=0 && idz >=0){
	    ptCell = dataArray[idx][idy][idz];
	    if(ptCell !=NULL) {
		cells.push_back(ptCell);
	    } else {
		this->addPoint(pt); ///Add fake point to initialize!
	    }
	} else {
	    //out of the map, we won't be coming back any time soon
	    complete = false;
	    final_point = origin+(float(i))*diff;
	    break;
	}
    }
    if(complete) final_point = origin+diff_;
    return true;

}

} //end namespace
