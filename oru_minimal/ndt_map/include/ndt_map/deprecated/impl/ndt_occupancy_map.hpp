#include <Eigen/Eigen>
#include <pcl/point_cloud.h>
#include <pcl/features/feature.h>

#include <string>
#include <climits>

#include<ndt_map/oc_tree.h>
#include<ndt_map/lazy_grid.h>
#include<ndt_map/cell_vector.h>

namespace lslgeneric
{

/** Main implemented method for now. Loads a single point cloud into an NDTOccupancyMap.
 */
template<typename PointT>
void NDTOccupancyMap<PointT>::loadPointCloud(const Eigen::Vector3d &origin, const pcl::PointCloud<PointT> &pc)
{
    if(index_ != NULL)
    {
        //std::cout<<"CLONE INDEX\n";
        SpatialIndex<PointT> *si = index_->clone();
        //cout<<"allocating index\n";
        if(!isFirstLoad_)
        {
            //std::cout<<"deleting old index\n";
            delete index_;
        }
        isFirstLoad_ = false;
        index_ = si;
    }
    else
    {
        //NULL index in constructor, abort!
        //ERR("constructor must specify a non-NULL spatial index\n");
        return;
    }

    if(index_ == NULL)
    {
        //ERR("Problem creating index, unimplemented method\n");
        return;
    }

    double maxDist = 0;//, distCeil = 200;

    typename pcl::PointCloud<PointT>::const_iterator it = pc.points.begin();
    Eigen::Vector3d centroid(0,0,0);
    int npts = 0;
    while(it!=pc.points.end())
    {
        Eigen::Vector3d d;
        if(std::isnan(it->x) ||std::isnan(it->y) ||std::isnan(it->z))
        {
            it++;
            continue;
        }
        d << it->x, it->y, it->z;
        centroid += d;
        it++;
        npts++;
    }

    centroid /= (double)npts;

    //Eigen::Vector4f centroid(0,0,0,0);
    //pcl::compute3DCentroid(pc,centroid);

    //compute distance to furthest point
    it = pc.points.begin();
    double maxz=0, minz=10000;

    while(it!=pc.points.end())
    {
        Eigen::Vector3d d;
        if(std::isnan(it->x) ||std::isnan(it->y) ||std::isnan(it->z))
        {
            it++;
            continue;
        }
        d << centroid(0)-it->x, centroid(1)-it->y, centroid(2)-it->z;
        double dist = d.norm();
        maxDist = (dist > maxDist) ? dist : maxDist;
        maxz = ((centroid(2)-it->z) > maxz) ? (centroid(2)-it->z) : maxz;
        minz = ((centroid(2)-it->z) < minz) ? (centroid(2)-it->z) : minz;
        it++;
    }
    fprintf(stderr,"minz=%lf maxz=%lf diff=%lf",minz,maxz,maxz-minz);

    // cout<<"Points = " <<pc.points.size()<<" maxDist = "<<maxDist<<endl;

    NDTCell<PointT> *ptCell = new NDTCell<PointT>();
    index_->setCellType(ptCell);
    delete ptCell;

    index_->setCenter(centroid(0),centroid(1),centroid(2));
    index_->setSize(3*maxDist,3*maxDist,2*(maxz-minz));

    /*
    		LazyGrid<PointT> *gr = dynamic_cast<LazyGrid<PointT>*>(index_);
    		if(gr!=NULL) {
    			fprintf(stderr,"Lazy Grid initialization!!\n");
    			gr->initializeAll();
    			fprintf(stderr,"Done!");
    		}else{
    			fprintf(stderr,"ERROR: This implementation is only for the Lazy Grid - Sorry!");
    			exit(1);
    		}
    */

//    ROS_INFO("centroid is %f,%f,%f", centroid(0),centroid(1),centroid(2));
//    ROS_INFO("maxDist is %lf", maxDist);

    it = pc.points.begin();
    while(it!=pc.points.end())
    {
        if(std::isnan(it->x) ||std::isnan(it->y) ||std::isnan(it->z))
        {
            it++;
            continue;
        }

        Eigen::Vector3d diff;
        diff << it->x-origin(0), it->y-origin(1), it->z-origin(2);
        double l = diff.norm();
        unsigned int N = l / resolution;

        diff = diff/(float)N;
        PointT pt;
        for(unsigned int i=0; i<N; i++)
        {
            pt.x = origin(0) + ((float)(i+1)) *diff(0);
            pt.y = origin(1) + ((float)(i+1)) *diff(1);
            pt.z = origin(2) + ((float)(i+1)) *diff(2);
            //fprintf(stderr,"Requesting point for (%.3f %.3f %.3f) target = (%.3f %.3f %.3f).... ",pt.x,pt.y,pt.z,it->x, it->y, it->z );
            NDTCell<PointT> *ptCell = dynamic_cast<NDTCell<PointT> *> (index_->getCellForPoint(pt));
            //fprintf(stderr,"Done\n");
            if(ptCell != NULL)
            {
                //fprintf(stderr,"Updating point for (%.3f %.3f %.3f) target = (%.3f %.3f %.3f).... ",pt.x,pt.y,pt.z,it->x, it->y, it->z );
                ptCell->updateEmpty();
            }
            else
            {
                //fprintf(stderr,"The cell is NULL\n");
            }
        }
        //pt.x = it->x; pt.y = it->y; pt.z = it->z;
        //NDTCell<PointT> *ptCell = dynamic_cast<NDTCell<PointT> *> (index_->getCellForPoint(pt));
        //ptCell->updateOccupancy(0.5);


        index_->addPoint(*it);
        it++;
    }

    isFirstLoad_ = false;
}


/**
* Adds a new cloud - NOTE::highly Experimental
*
*/
template<typename PointT>
void NDTOccupancyMap<PointT>::addPointCloud(const Eigen::Vector3d &origin, const pcl::PointCloud<PointT> &pc)
{
    if(isFirstLoad_)
    {
        loadPointCloud(origin, pc);
    }

    if(index_ == NULL)
    {
        //ERR("Problem creating index, unimplemented method\n");
        return;
    }
    typename pcl::PointCloud<PointT>::const_iterator it = pc.points.begin();
    it = pc.points.begin();
    while(it!=pc.points.end())
    {
        if(std::isnan(it->x) ||std::isnan(it->y) ||std::isnan(it->z))
        {
            it++;
            continue;
        }
        Eigen::Vector3d diff;
        diff << it->x-origin(0), it->y-origin(1), it->z-origin(2);
        double l = diff.norm();
        unsigned int N = l / (resolution);

        diff = diff/(float)N;
        PointT pt;
        for(unsigned int i=0; i<N; i++)
        {
            pt.x = origin(0) + ((float)(i+1)) *diff(0);
            pt.y = origin(1) + ((float)(i+1)) *diff(1);
            pt.z = origin(2) + ((float)(i+1)) *diff(2);
            NDTCell<PointT> *ptCell = dynamic_cast<NDTCell<PointT> *> (index_->getCellForPoint(pt));

            if(ptCell != NULL)
            {
                //fprintf(stderr,"Updating point for (%.3f %.3f %.3f) target = (%.3f %.3f %.3f).... ",pt.x,pt.y,pt.z,it->x, it->y, it->z );
                ptCell->updateEmpty();
            }
            else
            {
                //fprintf(stderr,"The cell is NULL\n");
            }
        }


        index_->addPoint(*it);
        it++;
    }
    isFirstLoad_ = false;
}





template<typename PointT>
void NDTOccupancyMap<PointT>::loadPointCloud(const pcl::PointCloud<PointT> &pc, const std::vector<std::vector<size_t> > &indices)
{

    loadPointCloud(pc);
    // Specific function related to CellVector
    CellVector<PointT> *cl = dynamic_cast<CellVector<PointT>*>(index_);
    if (cl != NULL)
    {
        for (size_t i = 0; i < indices.size(); i++)
        {
            cl->addCellPoints(pc, indices[i]);
        }

    }
    else
    {
        //ERR("loading point clouds using indices are currently supported in CellVector index_.");
    }
}

template<typename PointT>
void NDTOccupancyMap<PointT>::loadDepthImage(const cv::Mat& depthImage, DepthCamera<PointT> &cameraParams)
{
    pcl::PointCloud<PointT> pc;
    cameraParams.convertDepthImageToPointCloud(depthImage, pc);
    this->loadPointCloud(pc);
}

template<typename PointT>
pcl::PointCloud<PointT> NDTOccupancyMap<PointT>::loadDepthImageFeatures(const cv::Mat& depthImage, std::vector<cv::KeyPoint> &keypoints,
        size_t &supportSize, double maxVar, DepthCamera<PointT> &cameraParams, bool estimateParamsDI, bool nonMean)
{
    std::vector<cv::KeyPoint> good_keypoints;
    Eigen::Vector3d mean;
    PointT mn;
    pcl::PointCloud<PointT> cloudOut;
    CellVector<PointT> *cl = dynamic_cast<CellVector<PointT>*>(index_);
    if(cl==NULL)
    {
        std::cerr<<"wrong index type!\n";
        return cloudOut;
    }
    for (size_t i=0; i<keypoints.size(); i++)
    {
        if(!estimateParamsDI)
        {
            pcl::PointCloud<PointT> points;
            PointT center;
            cameraParams.computePointsAtIndex(depthImage,keypoints[i],supportSize,points,center);
            NDTCell<PointT> *ndcell = new NDTCell<PointT>();
            typename pcl::PointCloud<PointT>::iterator it = points.points.begin();
            while (it!= points.points.end() )
            {
                if(std::isnan(it->x) ||std::isnan(it->y) ||std::isnan(it->z))
                {
                    it++;
                    continue;
                }
                ndcell->addPoint(*it);
                it++;
            }
            ndcell->computeGaussian();
            if(ndcell->hasGaussian_)
            {
                Eigen::Vector3d evals = ndcell->getEvals();
                if(sqrt(evals(2)) < maxVar)
                {
                    if (nonMean)
                    {
                        if(std::isnan(center.x) ||std::isnan(center.y) ||std::isnan(center.z))
                        {
                            continue;
                        }
                        mn = center;
                    }
                    else
                    {
                        mean = ndcell->getMean();
                        mn.x = mean(0);
                        mn.y = mean(1);
                        mn.z = mean(2);
                    }
                    cloudOut.points.push_back(mn);
                    ndcell->setCenter(mn);
                    cl->addCell(ndcell);
                    good_keypoints.push_back(keypoints[i]);
                }
            }
        }
        else
        {
            assert(nonMean = false); // Not implemented / used.
            Eigen::Vector3d mean;
            Eigen::Matrix3d cov;
            cameraParams.computeParamsAtIndex(depthImage,keypoints[i],supportSize,mean,cov);
            NDTCell<PointT> *ndcell = new NDTCell<PointT>();
            ndcell->setMean(mean);
            ndcell->setCov(cov);

            if(ndcell->hasGaussian_)
            {
                Eigen::Vector3d evals = ndcell->getEvals();
                //std::cout<<evals.transpose()<<std::endl;
                if(sqrt(evals(2)) < maxVar)
                {
                    mean = ndcell->getMean();
                    mn.x = mean(0);
                    mn.y = mean(1);
                    mn.z = mean(2);
                    cloudOut.points.push_back(mn);
                    ndcell->setCenter(mn);
                    cl->addCell(ndcell);
                    good_keypoints.push_back(keypoints[i]);
                }
            }

        }
    }

    //TODO
    keypoints = good_keypoints;
    return cloudOut;
}

/** Helper function, computes the  NDTCells
  */
template<typename PointT>
void NDTOccupancyMap<PointT>::computeNDTCells(int cellupdatemode)
{
    CellVector<PointT> *cv = dynamic_cast<CellVector<PointT>*>(index_);

    typename SpatialIndex<PointT>::CellVectorItr it = index_->begin();

    while (it != index_->end())
    {
        NDTCell<PointT> *cell = dynamic_cast<NDTCell<PointT>*> (*it);
        //fprintf(stderr,"*");
        if(cell!=NULL)
        {
            //fprintf(stderr,"Something is not null!\n");
            cell->computeGaussian(cellupdatemode);
            if (cv!=NULL)
            {
                // Set the mean to the cell's centre.
                Eigen::Vector3d mean = cell->getMean();
                //cout << "mean : " << mean << std::endl;
                PointT pt;
                pt.x = mean[0];
                pt.y = mean[1];
                pt.z = mean[2];

                cell->setCenter(pt);
            }
        }
        else
        {
            //ERR("problem casting cell to NDT!\n");
        }
        it++;
    }

    LazyGrid<PointT> *lz = dynamic_cast<LazyGrid<PointT>*>(index_);
    if(lz!=NULL)
    {
        lz->initKDTree();
    }

    CellVector<PointT> *cl = dynamic_cast<CellVector<PointT>*>(index_);
    if(cl!=NULL)
    {
        cl->initKDTree();
    }
}

/** output method, saves as vrml the oct tree and all the ellipsoids
 */
template<typename PointT>
void NDTOccupancyMap<PointT>::writeToVRML(const char* filename)
{

    if(filename == NULL)
    {
        //ERR("problem outputing to vrml\n");
        return;
    }
    /*
       std::string fn(filename);
       fn = "oct_"+fn;

       FILE *fo = fopen(fn.c_str(),"w");
       if(fo == NULL) {
       ERR("problem outputing to vrml\n");
       return;
       }
       fprintf(fo,"#VRML V2.0 utf8\n");
       writeToVRML(fo,true);
       fclose(fo);
     */

    FILE *fout = fopen(filename,"w");
    if(fout == NULL)
    {
        //ERR("problem outputing to vrml\n");
        return;
    }
    fprintf(fout,"#VRML V2.0 utf8\n");
    writeToVRML(fout);
    fclose(fout);
}

/** helper output method
 */
template<typename PointT>
void NDTOccupancyMap<PointT>::writeToVRML(FILE* fout)
{
    if(fout == NULL)
    {
        //ERR("problem outputing to vrml\n");
        return;
    }

    //move the ellipsoid stuff to NDTCell
    typename SpatialIndex<PointT>::CellVectorItr it = index_->begin();
    while (it != index_->end())
    {
        NDTCell<PointT> *cell = dynamic_cast<NDTCell<PointT>*> (*it);
        //	double xs,ys,zs;
        if(cell!=NULL)
        {
            cell->writeToVRML(fout);
        }
        else
        {
            //	    ERR("problem casting cell to NDT!\n");
        }
        it++;
    }

}
template<typename PointT>
void NDTOccupancyMap<PointT>::writeToVRML(FILE* fout, Eigen::Vector3d col)
{
    if(fout == NULL)
    {
        //ERR("problem outputing to vrml\n");
        return;
    }

    //move the ellipsoid stuff to NDTCell
    typename SpatialIndex<PointT>::CellVectorItr it = index_->begin();
    while (it != index_->end())
    {
        NDTCell<PointT> *cell = dynamic_cast<NDTCell<PointT>*> (*it);
        if(cell!=NULL)
        {
            cell->writeToVRML(fout,col);
        }
        else
        {
        }
        it++;
    }
}

/// returns the current spatial index as a string (debugging function)
template<typename PointT>
std::string NDTOccupancyMap<PointT>::getMyIndexStr() const
{
    CellVector<PointT>* cl = dynamic_cast<CellVector<PointT> * >(index_);
    if(cl!=NULL)
    {
        return std::string("CellVector");
    }
    OctTree<PointT>* tr = dynamic_cast<OctTree<PointT>*>(index_);
    if(tr!=NULL)
    {
        return std::string("OctTree");
    }
    LazyGrid<PointT> *gr = dynamic_cast<LazyGrid<PointT>*>(index_);
    if(gr!=NULL)
    {
        return std::string("LazyGrid<PointT>");
    }

    return std::string("Unknown index type");
}

//computes the *negative log likelihood* of a single observation
template<typename PointT>
double NDTOccupancyMap<PointT>::getLikelihoodForPoint(PointT pt)
{
    //assert(false);
    double uniform=0.00100;
    NDTCell<PointT>* ndCell = NULL;
    OctTree<PointT>* tr = dynamic_cast<OctTree<PointT>*>(index_);

    if(tr==NULL)
    {
        LazyGrid<PointT> *gr = dynamic_cast<LazyGrid<PointT>*>(index_);
        if(gr==NULL)
        {
            //cout<<"bad index - getLikelihoodForPoint\n";
            return uniform;
        }
        ndCell = gr->getClosestNDTCell(pt);
    }
    else
    {
        ndCell = tr->getClosestNDTCell(pt);
    }
    if(ndCell == NULL) return uniform;

    double prob = ndCell->getLikelihood(pt);
    prob = (prob<0) ? 0 : prob; //uniform!! TSV
    return prob;
}

/*
//use trilinear interpolation from available immediate neighbors
template<typename PointT>
double NDTOccupancyMap<PointT>::getLikelihoodForPointWithInterpolation(PointT pt) {

    //ATM only for grid map
    //     tll------tlr
    //     /|       /|
    //    / |      / |
    //  tul------tur |    z
    //   | bll----|-blr   ^  y
    //   | /      | /     | /
    //   |/       |/      |/
    //  bul------bur      ---> x
    double uniform=0;//0.00100;
    Cell* cell = NULL;
    NDTCell<PointT>* ndCell = NULL;
    double cumProb = 0;
    double weight = 0;
    int evals = 1;

    LazyGrid<PointT> *gr = dynamic_cast<LazyGrid<PointT>*>(index_);
    if(gr==NULL) {
	//cout<<"bad index - getLikelihoodForPointWithInterpolation\n";
	return uniform;
    }
    cell = gr->getCellForPoint(pt);
    if(cell == NULL) return uniform;


    //get coordinates of cell
    int indXn, indYn, indZn;
    PointT centerGrid, sizeCell, centerCell;
    int sizeGridX, sizeGridY,sizeGridZ;
    centerCell = cell->getCenter();
    gr->getCenter(centerGrid.x,centerGrid.y,centerGrid.z);
    gr->getGridSize(sizeGridX,sizeGridY,sizeGridZ);
    gr->getCellSize(sizeCell.x,sizeCell.y,sizeCell.z);
    gr->getIndexForPoint(pt,indXn,indYn,indZn);

    double x,y,z;
    x = (pt.x - centerCell.x)/sizeCell.x;
    y = (pt.y - centerCell.y)/sizeCell.y;
    z = (pt.z - centerCell.z)/sizeCell.z;
    if(x <0 ) x = 0;
    if(y <0 ) y = 0;
    if(z <0 ) z = 0;
    if(x >1 ) x = 1;
    if(y >1 ) y = 1;
    if(z >1 ) z = 1;

    //bul
    double prob = 0;
    ndCell = dynamic_cast<NDTCell<PointT>*> (cell);
    if(ndCell != NULL) {
	prob = ndCell->getLikelihood(pt);
	prob = (prob<0) ? uniform : prob;
	weight = (1 - x + 1 - y + 1 - z)/(3.0);
	if(weight < 0) cerr<<weight<<endl;
	cumProb += prob*weight;
	//cout<<"\t"<<weight<<" "<<prob<<" --> "<<cumProb<<endl;
	evals++;
    }

    //tul
    Cell* c = gr->getCellAt(indXn,indYn,indZn+1);
    if(c != NULL) {
	ndCell = dynamic_cast<NDTCell<PointT>*> (c);
	if(ndCell != NULL) {
	    prob = ndCell->getLikelihood(pt);
	    prob = (prob<0) ? uniform : prob;
	    weight = (1 - x + 1 - y + z)/(3.0);
	    if(weight < 0) cerr<<weight<<endl;
	    cumProb += prob*weight;
	    //cout<<"\t"<<weight<<" "<<prob<<" --> "<<cumProb<<endl;
	    evals++;
	}
    }
    //tur
    c = gr->getCellAt(indXn+1,indYn,indZn+1);
    if(c != NULL) {
	ndCell = dynamic_cast<NDTCell<PointT>*> (c);
	if(ndCell != NULL) {
	    prob = ndCell->getLikelihood(pt);
	    prob = (prob<0) ? uniform : prob;
	    weight = (x + 1-y + z )/(3.0);
	    if(weight < 0) cerr<<weight<<endl;
	    cumProb += prob*weight;
	    //cout<<"\t"<<weight<<" "<<prob<<" --> "<<cumProb<<endl;
	    evals++;
	}
    }
    //tll
    c = gr->getCellAt(indXn,indYn+1,indZn+1);
    if(c != NULL) {
	ndCell = dynamic_cast<NDTCell<PointT>*> (c);
	if(ndCell != NULL) {
	    prob = ndCell->getLikelihood(pt);
	    prob = (prob<0) ? uniform : prob;
	    weight = (1-x + y + z )/(3.0);
	    if(weight < 0) cerr<<weight<<endl;
	    cumProb += prob*weight;
	    //cout<<"\t"<<weight<<" "<<prob<<" --> "<<cumProb<<endl;
	    evals++;
	}
    }
    //tlr
    c = gr->getCellAt(indXn+1,indYn+1,indZn+1);
    if(c != NULL) {
	ndCell = dynamic_cast<NDTCell<PointT>*> (c);
	if(ndCell != NULL) {
	    prob = ndCell->getLikelihood(pt);
	    prob = (prob<0) ? uniform : prob;
	    weight = (x + y + z )/(3.0);
	    if(weight < 0) cerr<<weight<<endl;
	    cumProb += prob*weight;
	    //cout<<"\t"<<weight<<" "<<prob<<" --> "<<cumProb<<endl;
	    evals++;
	}
    }
    //bur
    c = gr->getCellAt(indXn+1,indYn,indZn);
    if(c != NULL) {
	ndCell = dynamic_cast<NDTCell<PointT>*> (c);
	if(ndCell != NULL) {
	    prob = ndCell->getLikelihood(pt);
	    prob = (prob<0) ? uniform : prob;
	    weight = (x + 1-y + 1-z )/(3.0);
	    if(weight < 0) cerr<<weight<<endl;
	    cumProb += prob*weight;
	    //cout<<"\t"<<weight<<" "<<prob<<" --> "<<cumProb<<endl;
	    evals++;
	}
    }
    //bll
    c = gr->getCellAt(indXn,indYn+1,indZn);
    if(c != NULL) {
	ndCell = dynamic_cast<NDTCell<PointT>*> (c);
	if(ndCell != NULL) {
	    prob = ndCell->getLikelihood(pt);
	    prob = (prob<0) ? uniform : prob;
	    weight = (1-x + y + 1-z )/(3.0);
	    if(weight < 0) cerr<<weight<<endl;
	    cumProb += prob*weight;
	    //cout<<"\t"<<weight<<" "<<prob<<" --> "<<cumProb<<endl;
	    evals++;
	}
    }
    //blr
    c = gr->getCellAt(indXn+1,indYn+1,indZn);
    if(c != NULL) {
	ndCell = dynamic_cast<NDTCell<PointT>*> (c);
	if(ndCell != NULL) {
	    prob = ndCell->getLikelihood(pt);
	    prob = (prob<0) ? uniform : prob;
	    weight = (x + y + 1-z )/(3.0);
	    if(weight < 0) cerr<<weight<<endl;
	    cumProb += prob*weight;
	    //cout<<"\t"<<weight<<" "<<prob<<" --> "<<cumProb<<endl;
	    evals++;
	}
    }

    //cout<<"== "<<cumProb<<endl;
    return cumProb;
}
*/

template<typename PointT>
std::vector<NDTCell<PointT>*> NDTOccupancyMap<PointT>::getCellsForPoint(const PointT pt, double radius)
{
    //assert(false);
    std::vector<NDTCell<PointT>*> cells;
    OctTree<PointT>* tr = dynamic_cast<OctTree<PointT>*>(index_);
    if(tr==NULL)
    {
        LazyGrid<PointT> *gr = dynamic_cast<LazyGrid<PointT>*>(index_);
        if(gr==NULL)
        {
            //cout<<"bad index - getCellsForPoint\n";
            return cells;
        }
        cells = gr->getClosestNDTCells(pt,radius);
        return cells;
    }
    //TODO:implement for ocTree
    return cells;
}

template<typename PointT>
bool NDTOccupancyMap<PointT>::getCellForPoint(const PointT &pt, NDTCell<PointT> *&out_cell)
{

#if 0
    OctTree<PointT>* tr = dynamic_cast<OctTree<PointT>*>(index);
    if(tr==NULL)
    {
        LazyGrid<PointT> *gr = dynamic_cast<LazyGrid<PointT>*>(index);
        if(gr==NULL)
        {
            CellVector *cl = dynamic_cast<CellVector*>(index);
            if(cl==NULL)
            {
                cout<<"bad index - getCellForPoint\n";
                return false;
            }
            out_cell = cl->getClosestNDTCell(pt);
            return true;
        }
        out_cell = gr->getClosestNDTCell(pt);
        return true;
    }
    out_cell = tr->getClosestNDTCell(pt);
    return true;
#endif
    CellVector<PointT> *cl = dynamic_cast<CellVector<PointT>*>(index_);
    if(cl!=NULL)
    {
        out_cell = cl->getClosestNDTCell(pt);
        return true;
    }
    OctTree<PointT>* tr = dynamic_cast<OctTree<PointT>*>(index_);
    if(tr!=NULL)
    {
        out_cell = tr->getClosestNDTCell(pt);
        return true;
    }

    LazyGrid<PointT> *gr = dynamic_cast<LazyGrid<PointT>*>(index_);
    if(gr!=NULL)
    {
        //fprintf(stderr,"Lazy!!\n");
        out_cell = gr->getClosestNDTCell(pt);
        return true;
    }
    //cout<<"bad index - getCellForPoint\n";
    return false;

    /*
    if(c==NULL) {
    //we have to find the closest leaf cell
    c = tr->getClosestLeafCell(pt);
    }
    if(c==NULL) {
    cout<<"null cell\n";
    return false;
    }
    NDTCell* ndCell = dynamic_cast<NDTCell*>(c);
    if(ndCell == NULL) {
    cout<<"not ndt cell\n";
    return false;
    }
    if(!ndCell->hasGaussian) {
    Eigen::Vector3d dim;
    ndCell->getDimensions(dim(0),dim(1),dim(2));
    typename SpatialIndex<PointT>::CellVectorItr neigh;
    index->getNeighbors(pt,20*dim.norm(),neigh);
    double minDist = INT_MAX;
    cout<<"neighs "<<neigh.size()<<endl;
    for(int i=0; i<neigh.size(); i++) {
        NDTCell *n = dynamic_cast<NDTCell*>(neigh[i]);
        if(n==NULL) continue;
        if(n->hasGaussian) {
    	double d = lslgeneric::geomDist(n->getCenter(),pt);
    	if(d<minDist) {
    	    ndCell = n;
    	}
        }
    }
    if(!ndCell->hasGaussian) {
        cout<<"no gaussian\n";
        return false;
    }
    }
    out_cell = ndCell;
    return true;
    */
}

template<typename PointT>
void NDTOccupancyMap<PointT>::debugToVRML(const char* fname, pcl::PointCloud<PointT> &pc)
{

    FILE* fout = fopen(fname, "w");

    fprintf(fout,"#VRML V2.0 utf8\n");
    this->writeToVRML(fout);
    lslgeneric::writeToVRML(fout,pc,Eigen::Vector3d(1,0,0));

    fprintf(fout,"Shape {\n\tgeometry IndexedLineSet {\n\tcoord Coordinate {\n\t point [\n\t");

    int n_lines = 0;
    PointT centerCell;
    for(size_t i=0; i<pc.points.size(); i++)
    {
        NDTCell<PointT>* link;
        if(this->getCellForPoint(pc.points[i], link))
        {
            if(link == NULL) continue;
            centerCell = link->getCenter();
            if(link->hasGaussian_)
            {
                centerCell.x = link->getMean()(0);
                centerCell.y = link->getMean()(1);
                centerCell.z = link->getMean()(2);
            }
            fprintf(fout,"%lf %lf %lf\n\t%lf %lf %lf\n\t",
                    pc.points[i].x, pc.points[i].y,pc.points[i].z,
                    centerCell.x, centerCell.y, centerCell.z);
            n_lines++;
        }
    }

    fprintf(fout, "]\n\t}\n\tcoordIndex [\n\t");
    for(int i = 0; i<n_lines; i++)
    {
        fprintf(fout,"%d, %d, -1\n\t",2*i, 2*i+1);
    }
    fprintf(fout, "]\n}\n}\n");


    fclose(fout);
}

template<typename PointT>
std::vector<NDTCell<PointT>*> NDTOccupancyMap<PointT>::pseudoTransformNDT(Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> T)
{

    std::vector<NDTCell<PointT>*> ret;
    typename SpatialIndex<PointT>::CellVectorItr it = index_->begin();
    while (it != index_->end())
    {
        NDTCell<PointT> *cell = dynamic_cast<NDTCell<PointT>*> (*it);
        if(cell!=NULL)
        {
            if(cell->hasGaussian_)
            {
                Eigen::Vector3d mean = cell->getMean();
                Eigen::Matrix3d cov = cell->getCov();
                mean = T*mean;
                cov = T.rotation().transpose()*cov*T.rotation();
                NDTCell<PointT>* nd = (NDTCell<PointT>*)cell->clone();
                nd->setMean(mean);
                nd->setCov(cov);
                ret.push_back(nd);
            }
        }
        else
        {
            //ERR("problem casting cell to NDT!\n");
        }
        it++;
    }
    return ret;
}

template<typename PointT>
NDTCell<PointT>*
NDTOccupancyMap<PointT>::getCellIdx(unsigned int idx)
{
    CellVector<PointT> *cl = dynamic_cast<CellVector<PointT>*>(index_);
    if (cl != NULL)
    {
        return cl->getCellIdx(idx);
    }
    return NULL;
}


/**
* Returns dynamic cells from the map
 */
template<typename PointT>
std::vector<lslgeneric::NDTCell<PointT>*>  NDTOccupancyMap<PointT>::getDynamicCells(unsigned int Timescale, float threshold)
{
    std::vector<NDTCell<PointT>*> ret;
    typename SpatialIndex<PointT>::CellVectorItr it = index_->begin();
    while (it != index_->end())
    {
        NDTCell<PointT> *cell = dynamic_cast<NDTCell<PointT>*> (*it);
        if(cell!=NULL)
        {

            if(cell->getDynamicLikelihood(Timescale)>threshold)
            {
                NDTCell<PointT>* nd = (NDTCell<PointT>*)cell->copy();
                ret.push_back(nd);
            }
        }
        else
        {
        }
        it++;
    }
    return ret;
}


template<typename PointT>
//std::vector<lslgeneric::NDTCell<pcl::PointXYZ>*> NDTOccupancyMap<PointT>::getAllCells(){
std::vector<lslgeneric::NDTCell<PointT>*> NDTOccupancyMap<PointT>::getAllCells()
{

    std::vector<NDTCell<PointT>*> ret;
    typename SpatialIndex<PointT>::CellVectorItr it = index_->begin();
    while (it != index_->end())
    {
        NDTCell<PointT> *cell = dynamic_cast<NDTCell<PointT>*> (*it);
        if(cell!=NULL)
        {
            if(cell->hasGaussian_)
            {
                NDTCell<PointT>* nd = (NDTCell<PointT>*)cell->copy();
                ret.push_back(nd);
            }
        }
        else
        {
        }
        it++;
    }
    return ret;
}

template<typename PointT>
int NDTOccupancyMap<PointT>::numberOfActiveCells()
{
    int ret = 0;
    if(index_ == NULL) return ret;
    typename SpatialIndex<PointT>::CellVectorItr it = index_->begin();
    while (it != index_->end())
    {
        NDTCell<PointT> *cell = dynamic_cast<NDTCell<PointT>*> (*it);
        if(cell!=NULL)
        {
            if(cell->hasGaussian_)
            {
                ret++;
            }
        }
        it++;
    }
    return ret;
}

}
