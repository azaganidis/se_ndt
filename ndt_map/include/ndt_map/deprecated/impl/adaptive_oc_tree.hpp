#include <boost/math/distributions/chi_squared.hpp>
#include <Eigen/Eigen>
//#include <ros/ros.h>

//using boost::math::chi_squared;
//using boost::math::cdf;

namespace lslgeneric
{

template <typename PointT>
void AdaptiveOctTree<PointT>::setParameters(bool _useDornikHansen,
        bool _useFlatness,
        double _RSS_THRESHOLD,
        double _DH_SIGNIFICANCE_LVL,
        double _MIN_CELL_SIZE,
        double _FLAT_FACTOR,
        double _BIG_CELL_SIZE,
        double _SMALL_CELL_SIZE
                                           )
{

    this->BIG_CELL_SIZE = _BIG_CELL_SIZE;
    this->SMALL_CELL_SIZE = _SMALL_CELL_SIZE;

    useDornikHansen      = _useDornikHansen;
    useFlatness          = _useFlatness;
    RSS_THRESHOLD        = _RSS_THRESHOLD;
    DH_SIGNIFICANCE_LVL  = _DH_SIGNIFICANCE_LVL;
    MIN_CELL_SIZE        = _MIN_CELL_SIZE;
    FLAT_FACTOR          = _FLAT_FACTOR;

    parametersSet = true;
}
/** empty! default constructor
  */
template <typename PointT>
AdaptiveOctTree<PointT>::AdaptiveOctTree() : OctTree<PointT>()
{
    if(!parametersSet)
    {
        setParameters();
    }
}

/** constructor, calls parent OctTree constructor
  */
template <typename PointT>
AdaptiveOctTree<PointT>::AdaptiveOctTree(pcl::PointXYZ center, double xsize, double ysize,
        double zsize, NDTCell<PointT>* type, OctTree<PointT> *_parent, unsigned int _depth) :
    OctTree<PointT>(center,xsize,ysize,zsize,type,_parent,_depth)
{
    if(!parametersSet)
    {
        setParameters();
    }
}

/** empty destructor, all data are deallocated by parent class
  */
template <typename PointT>
AdaptiveOctTree<PointT>::~AdaptiveOctTree()
{
}

/**
  finds all leafs of this tree and fills the vector of pointers to the leafs
  */
template <typename PointT>
void AdaptiveOctTree<PointT>::computeTreeLeafs()
{
    if(this->isLeaf())
    {
        myTreeLeafs.push_back(this);
        return;
    }

    myTreeLeafs.clear();
    std::vector<OctTree<PointT>*> next;
    next.push_back(this);

    while(next.size()>0)
    {
        OctTree<PointT> *cur = next.front();
        if(cur!=NULL)
        {
            if(cur->isLeaf())
            {
                myTreeLeafs.push_back(cur);
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

/**
go through all leafs and split the ones with high residuals
  */
template <typename PointT>
void AdaptiveOctTree<PointT>::postProcessPoints()
{

    //compute leafs as OctTree*
    computeTreeLeafs();
    //cout<<"leafs :"<<myTreeLeafs.size()<<endl;

    for(unsigned int i=0; i<myTreeLeafs.size(); i++)
    {
        NDTCell<PointT> * nd = dynamic_cast<NDTCell<PointT>*>((myTreeLeafs[i])->myCell_);
        if(nd == NULL) continue;
        if(useDornikHansen)
        {
            double significance = computeDornikHansen(nd);
            if(significance < 0)
            {
                //there were not enough points, let's clear the cell!
                //(myTreeLeafs[i])->myCell->points.clear();
                continue;
            }
            if(significance < DH_SIGNIFICANCE_LVL)
            {
                //split the leafs
                std::vector<OctTree<PointT>*> newLeafs = splitTree(myTreeLeafs[i]);
                myTreeLeafs.insert(myTreeLeafs.end(),newLeafs.begin(),newLeafs.end());
            }
        }
        else if(useFlatness)
        {
            nd->computeGaussian();
            if(!nd->hasGaussian_)
            {
                continue;
            }

            Eigen::Vector3d evals = nd->getEvals();
            int idMin, idMax;
            double minEval = evals.minCoeff(&idMin);
            double maxEval = evals.maxCoeff(&idMax);
            int idMiddle = -1;
            for(int j=0; j<3; j++)
            {
                if(j!=idMin && j!=idMax)
                {
                    idMiddle =  j;
                }
            }
            if(idMiddle < 0) continue;

            if(minEval*FLAT_FACTOR > evals(idMiddle))
            {
                std::vector<OctTree<PointT>*> newLeafs = splitTree(myTreeLeafs[i]);
                myTreeLeafs.insert(myTreeLeafs.end(),newLeafs.begin(),newLeafs.end());
            }


        }
        else
        {
            double rss = computeResidualSquare(nd);
            if(rss > RSS_THRESHOLD)
            {
                //split the leafs
                std::vector<OctTree<PointT>*> newLeafs = splitTree(myTreeLeafs[i]);
                myTreeLeafs.insert(myTreeLeafs.end(),newLeafs.begin(),newLeafs.end());
            }
        }
    }

    this->leafsCached_ = false;
}

/**
  performs the Dornik-Hansen Omnibus normality test
*/
template <typename PointT>
double AdaptiveOctTree<PointT>::computeDornikHansen(NDTCell<PointT> *cell)
{
    double pval = 1;
    double Ep = 0;
    //test statistics breaks down for n<=7
    if(cell->points_.size() <= 7) return -1;
    cell->computeGaussian();

    //degree is 2*number_of_dimensions
    boost::math::chi_squared dist(6);

    Eigen::Vector3d mean = cell->getMean();
    Eigen::Matrix3d C = cell->getCov();

    Eigen::MatrixXd Xhat(cell->points_.size(),3); //, tmpMat;

    for(unsigned int i=0; i< cell->points_.size(); i++)
    {
        Xhat(i,0) = cell->points_[i].x - mean(0);
        Xhat(i,1) = cell->points_[i].y - mean(1);
        Xhat(i,2) = cell->points_[i].z - mean(2);
    }

    //Compute transform for observed points
    Eigen::Matrix3d V;
    Eigen::Matrix3d H;
    Eigen::Vector3d lambda;
    Eigen::Matrix3d L;
    Eigen::Matrix3d R1;

    Eigen::MatrixXd R;

    V(0,0) = 1/sqrt(C(0,0));
    V(1,1) = 1/sqrt(C(1,1));
    V(2,2) = 1/sqrt(C(2,2));

    C = V*C*V;
    Eigen::EigenSolver<Eigen::Matrix3d> eig(C);
    H = eig.eigenvectors().real();
    lambda = eig.eigenvalues().real();

    //covariance is not positive semidefinate
    if(lambda.minCoeff() <= 0) return -1;

    L(0,0) = 1/sqrt(lambda(0));
    L(1,1) = 1/sqrt(lambda(1));
    L(2,2) = 1/sqrt(lambda(2));

    //transform observations
    R1 = H*L*H.transpose()*V;
    R = R1*Xhat.transpose();

    //compute skewness and kurtois of new observations (in each dimension)
    //samples are zero mean, compute for each dimension standard deviation
    Ep = 0.0;
    double n = R.cols();
    for(unsigned int dim = 0; dim < 3; dim ++)
    {
        double m2 = 0, m3 = 0, m4 =0;
        double b1 = 0, b2;
        double beta, omega2, gamma, y;
        double gamma2, a,c, k, alpha, chi;
        double z1,z2;

        for(unsigned int i=0; i<R.cols(); i++)
        {
            m2 += pow(R(dim,i),2);
            m3 += pow(R(dim,i),3);
            m4 += pow(R(dim,i),4);
        }
        m2 /= R.cols();
        m3 /= R.cols();
        m4 /= R.cols();
//	cout <<"dim "<<dim<<" m2 "<<m2<<" m3 "<<m3<<" m4 "<<m4<<endl;
        b1 = m3/(pow(m2,1.5));
        b2 = m4/(pow(m2,2));

//	cout<<"b1 "<<b1<<" b2 "<<b2<<endl;

        //compute Z1 and Z2
        beta = 3*(n*n + 27*n -70)*(n+1)*(n+3)/((n-2)*(n+5)*(n+7)*(n+9));
        omega2 = -1+sqrt(2*(beta-1));
        gamma = 1/sqrt(log(sqrt(omega2)));
        y = b1*sqrt((omega2-1)*(n+1)*(n+3)/(12*(n-2)));
        z1 = gamma*log(y+sqrt(y*y+1));

        gamma2 = (n-3)*(n+1)*(n*n+15*n-4);
        a = (n-2)*(n+5)*(n+7)*(n*n+27*n-70)/(6*gamma2);
        c = (n-7)*(n+5)*(n+7)*(n*n+2*n-5)/(6*gamma2);
        k = (n+5)*(n+7)*(pow(n,3)+37*n*n+11*n-313)/(12*gamma2);
        alpha = a + b1*b1*c;
        chi = (b2-1-b1*b1)*2*k;
        z2 = (pow((chi/(2*alpha)),(1./3.))-1+(1./(9*alpha)))*sqrt(9*alpha);

//	cout<<"z1: "<<z1<<" z2: "<<z2<<endl;

        //compute Ep
        Ep += z1*z1 + z2*z2;
    }

    //compute probability from chi square cdf
    pval = 1-boost::math::cdf(dist,Ep);

//    cout<<"Ep "<<Ep<<endl;
//    cout<<"P: "<<pval<<endl;
    if(pval>DH_SIGNIFICANCE_LVL)
    {
        double dx,dy,dz;
        cell->getDimensions(dx,dy,dz);
        std::cout<<"final split was at ("<<dx<<","<<dy<<","<<dz<<"); pval is "<<pval<<std::endl;
    }
    return 0;

}
/**
fits a 3d gaussian in the cell and computes the residual squares sum
*/
template <typename PointT>
double AdaptiveOctTree<PointT>::computeResidualSquare(NDTCell<PointT> *cell)
{
    double rss = 0; //residual sum squared
    double meanResidual = 0;
    double residualVar = 0;
    if(cell->points_.size() <= 3) return 0;

    cell->computeGaussian();

    Eigen::Vector3d cur, curProd;
    Eigen::Matrix3d cov = cell->getCov();
    Eigen::Vector3d mean = cell->getMean();
    //Eigen::LLT<Eigen::Matrix3d> lltOfCov = cov.llt();

    for(unsigned int i=0; i< cell->points_.size(); i++)
    {

        cur(0) = cell->points_[i].x - mean(0);
        cur(1) = cell->points_[i].y - mean(1);
        cur(2) = cell->points_[i].z - mean(2);

        //lltOfCov.solve(cur,&curProd);
        //rss += cur.dot(curProd);
        rss += cur.dot(cur);
        meanResidual += cur.norm()/cell->points_.size();
    }

    for(unsigned int i=0; i< cell->points_.size(); i++)
    {
        cur(0) = cell->points_[i].x - mean(0);
        cur(1) = cell->points_[i].y - mean(1);
        cur(2) = cell->points_[i].z - mean(2);

        residualVar += pow(cur.norm()-meanResidual,2)/(cell->points_.size()-1);
    }
    double bic = rss/residualVar + log(cell->points_.size());
//    ROS_INFO("rss %lf mean %lf var %lf bic %lf",rss,meanResidual,residualVar,bic);
    return bic;
}

/**
  splits a cell and returns a vector of the newly created children
  iterates points downwards
  */
template <typename PointT>
std::vector<OctTree<PointT>*> AdaptiveOctTree<PointT>::splitTree(OctTree<PointT> *octLeaf)
{
    std::vector<OctTree<PointT>*> newLeafs;

    if(octLeaf->isLeaf())
    {
        double xs,ys,zs;
        octLeaf->myCell_->getDimensions(xs,ys,zs);

        double cellSize = (xs+ys+zs)/3.; //average for now

        if(octLeaf->depth_ > this->MAX_DEPTH || cellSize <= this->MIN_CELL_SIZE )
        {
            //just store point, we can't split any more
            return newLeafs;
        }

        pcl::PointXYZ myCenter = octLeaf->myCell_->getCenter();

        //branch leaf
        for(unsigned int it=0; it<8; it++)
        {

            pcl::PointXYZ newCenter;

            //computes the center of the it'th child
            newCenter.x = (myCenter.x + pow(-1.,it/4)*xs/4.);
            newCenter.y = (myCenter.y + pow(-1.,it/2)*ys/4.);
            newCenter.z = (myCenter.z + pow(-1.,(it+1)/2)*zs/4.);

            octLeaf->children_[it] = new OctTree<PointT>(newCenter,xs/2,ys/2,
                    zs/2, octLeaf->myCell_, this, this->depth_+1);
            newLeafs.push_back(octLeaf->children_[it]);
        }
        //add current points
        for(unsigned int jt=0; jt<octLeaf->myCell_->points_.size(); jt++)
        {
            size_t ind = octLeaf->getIndexForPoint(octLeaf->myCell_->points_[jt]);
            octLeaf->children_[ind]->addPoint(octLeaf->myCell_->points_[jt]);
        }
        octLeaf->leaf_=false;
        octLeaf->myCell_->points_.clear();
    }

    return newLeafs;
}

/**
  creates an oct tree with the same parameters.
  \note the points are not copied in teh returned instance
  */
template <typename PointT>
SpatialIndex<PointT>* AdaptiveOctTree<PointT>::clone()
{
    if(this->myCell_ == NULL)
    {
        return new AdaptiveOctTree();
    }
    double sx,sy,sz;
    this->myCell_->getDimensions(sx,sy,sz);
    AdaptiveOctTree<PointT> *tr = new AdaptiveOctTree<PointT>(this->myCell_->getCenter(),sx,sy,sz,this->myCell_);
    return tr;
}

}
