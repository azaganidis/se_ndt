#include "ndt_registration/ndt_matcher_d2d_2d.h"
#include "ndt_map/ndt_cell.h"
#include "ndt_map/lazy_grid.h"
#include "ndt_map/pointcloud_utils.h"
#include "Eigen/Eigen"
#include <fstream>
#include <omp.h>
#include <sys/time.h>
namespace lslgeneric
{

//#define DO_DEBUG_PROC

void NDTMatcherD2D_2D::init(bool _isIrregularGrid,
        bool useDefaultGridResolutions, std::vector<double> _resolutions)
{
    Jest.setZero();
    Jest.block<2,2>(0,0).setIdentity();
    Hest.setZero();
    Zest.setZero();
    ZHest.setZero();

    isIrregularGrid = _isIrregularGrid;
    if(useDefaultGridResolutions)
    {
        resolutions.push_back(0.2);
        resolutions.push_back(0.5);
        resolutions.push_back(1);
        resolutions.push_back(2);
    }
    else
    {
        resolutions = _resolutions;
    }

    current_resolution = 0.1; // Argggg!!! This is very important to have initiated! (one day debugging later) :-) TODO do we need to set it to anything better?
    lfd1 = 1; //lfd1/(double)sourceNDT.getMyIndex()->size(); //current_resolution*2.5;
    lfd2 = 0.05; //0.1/current_resolution;
    ITR_MAX = 50;
    DELTA_SCORE = 10e-3*current_resolution;
    step_control = true;

}

bool NDTMatcherD2D_2D::match( pcl::PointCloud<pcl::PointXYZ>& target,
        pcl::PointCloud<pcl::PointXYZ>& source,
        Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T ,
        bool useInitialGuess)
{

    struct timeval tv_start, tv_end;
    struct timeval tv_start0, tv_end0;
    double time_load =0, time_match=0, time_combined=0;

    gettimeofday(&tv_start0,NULL);

    //initial guess
    pcl::PointCloud<pcl::PointXYZ> sourceCloud = source;
    
    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> Temp, Tinit;
    Tinit.setIdentity();
    if(useInitialGuess)
    {
        lslgeneric::transformPointCloudInPlace(T,sourceCloud);
	Tinit = T;
    }

    T.setIdentity();
    bool ret = false;

#if 0
    if(isIrregularGrid)
    {

        OctTree<PointTarget> pr1;
        NDTMap<PointTarget> targetNDT( &pr1 );
        targetNDT.loadPointCloud( target );
        targetNDT.computeNDTCells();

        OctTree<PointSource> pr2;
        NDTMap<PointSource> sourceNDT( &pr2 );
        sourceNDT.loadPointCloud( source );
        sourceNDT.computeNDTCells();

        ret = this->match( targetNDT, sourceNDT, T );

    }
    else
#endif
    {

        //iterative regular grid
        for(int r_ctr = resolutions.size()-1; r_ctr >=0;  r_ctr--)
        {

            current_resolution = resolutions[r_ctr];

            LazyGrid prototypeSource(current_resolution);
            LazyGrid prototypeTarget(current_resolution);

            gettimeofday(&tv_start,NULL);
            NDTMap targetNDT( &prototypeTarget );
            targetNDT.loadPointCloud( target );
            targetNDT.computeNDTCells();

            NDTMap sourceNDT( &prototypeSource );
            sourceNDT.loadPointCloud( sourceCloud );
            sourceNDT.computeNDTCells();
            gettimeofday(&tv_end,NULL);

            time_load += (tv_end.tv_sec-tv_start.tv_sec)*1000.+(tv_end.tv_usec-tv_start.tv_usec)/1000.;
            Temp.setIdentity();

            gettimeofday(&tv_start,NULL);
            ret = this->match( targetNDT, sourceNDT, Temp );
            lslgeneric::transformPointCloudInPlace(Temp,sourceCloud);
            gettimeofday(&tv_end,NULL);

            time_match += (tv_end.tv_sec-tv_start.tv_sec)*1000.+(tv_end.tv_usec-tv_start.tv_usec)/1000.;

            //transform moving
            T = Temp*T;

#ifdef DO_DEBUG_PROC
            std::cout<<"RESOLUTION: "<<current_resolution<<std::endl;
            std::cout<<"rotation   : "<<Temp.rotation().eulerAngles(0,1,2).transpose()<<std::endl;
            std::cout<<"translation: "<<Temp.translation().transpose()<<std::endl;
            std::cout<<"--------------------------------------------------------\nOverall Transform:\n";
            std::cout<<"rotation   : "<<T.rotation().eulerAngles(0,1,2).transpose()<<std::endl;
            std::cout<<"translation: "<<T.translation().transpose()<<std::endl;

#endif
        }
    }
    gettimeofday(&tv_end0,NULL);
    time_combined = (tv_end0.tv_sec-tv_start0.tv_sec)*1000.+(tv_end0.tv_usec-tv_start0.tv_usec)/1000.;
    std::cout<<"load: "<<time_load<<" match "<<time_match<<" combined "<<time_combined<<std::endl;
    if(useInitialGuess)
    {
	T = T*Tinit;
    }
    return ret;
}

bool NDTMatcherD2D_2D::match( NDTMap& targetNDT,
        NDTMap& sourceNDT,
        Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T ,
        bool useInitialGuess)
{

    //locals
    bool convergence = false;
    int itr_ctr = 0;
    double step_size = 1;
    Eigen::Matrix<double,3,1>  pose_increment_v, scg;
    //Eigen::MatrixXd Hessian(6,6), score_gradient(6,1), H(3,3); //column vectors, pose_increment_v(6,1)
    Eigen::MatrixXd H(3,3), score_gradient_2d(3,1);
    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> TR, Tbest;
    Eigen::Vector3d transformed_vec, mean;
    bool ret = true;
    double score_best = INT_MAX;
    if(!useInitialGuess)
    {
        T.setIdentity();
    }
    Tbest = T;

    std::vector<NDTCell*> nextNDT = sourceNDT.pseudoTransformNDT(T);
    while(!convergence)
    {
        TR.setIdentity();
        H.setZero();
        score_gradient_2d.setZero();

//        double score_here = derivativesNDT(nextNDT,targetNDT,score_gradient,Hessian,true);
        double score_here = derivativesNDT_2d(nextNDT,targetNDT,score_gradient_2d,H,true);
	scg = score_gradient_2d;
//	std::cout<<"itr "<<itr_ctr<<" score "<<score_here<<std::endl;
	if(score_here < score_best) 
	{
	    Tbest = T;
	    score_best = score_here;
//	    std::cout<<"best score "<<score_best<<" at "<<itr_ctr<<std::endl;
	}

        if (score_gradient_2d.norm()<= 10e-2*DELTA_SCORE)
        {
            std::cout<<"\%gradient vanished\n";
            for(unsigned int i=0; i<nextNDT.size(); i++)
            {
                if(nextNDT[i]!=NULL)
                    delete nextNDT[i];
            }
	    if(score_here > score_best) 
	    {
		std::cout<<"crap iterations, best was "<<score_best<<" last was "<<score_here<<std::endl;
		T = Tbest;
	    }
            return true;
        }
        //std::cout<<"Hh(:,:,"<<itr_ctr+1<<")  =  ["<< H<<"];\n"<<std::endl;				  //
        //std::cout<<"gradh (:,"<<itr_ctr+1<<")= ["<<scg.transpose()<<"];"<<std::endl;         //
        
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,3,3> > Sol (H);
        Eigen::Matrix<double,3,1> evals = Sol.eigenvalues().real();
        double minCoeff = evals.minCoeff();
        double maxCoeff = evals.maxCoeff();
        if(minCoeff < 10e-3 )   //|| evals.minCoeff() < 0.001*evals.maxCoeff()) {
        {
            Eigen::Matrix<double,3,3> evecs = Sol.eigenvectors().real();
            double regularizer = score_gradient_2d.norm();
	    regularizer = regularizer + minCoeff > 0 ? regularizer : 0.001*maxCoeff - minCoeff;
            //double regularizer = 0.001*maxCoeff - minCoeff;
            Eigen::Matrix<double,3,1> reg;
            //ugly
            reg<<regularizer,regularizer,regularizer;
            evals += reg;
            Eigen::Matrix<double,3,3> Lam;
            Lam = evals.asDiagonal();
            H = evecs*Lam*(evecs.transpose());
            //std::cerr<<"regularizing\n";
        }
//        std::cout<<"Hh(:,:,"<<itr_ctr+1<<")  =  ["<< H<<"];\n"<<std::endl;				  //
//        std::cout<<"gradh (:,"<<itr_ctr+1<<")= ["<<scg.transpose()<<"];"<<std::endl;         //

        pose_increment_v = -H.ldlt().solve(scg);
        double dginit = pose_increment_v.dot(scg);
        if(dginit > 0)
        {
            std::cout<<"\%can't decrease in this direction any more, done \n";
            for(unsigned int i=0; i<nextNDT.size(); i++)
            {
                if(nextNDT[i]!=NULL)
                    delete nextNDT[i];
            }
	    if(score_here > score_best) 
	    {
		//std::cout<<"crap iterations, best was "<<score_best<<" last was "<<score_here<<std::endl;
		T = Tbest;
	    }
            return true;
        }
	if(step_control) {
	    step_size = lineSearch2D(pose_increment_v,nextNDT,targetNDT);
	} else {
	    step_size = 1;
	}
        pose_increment_v = step_size*pose_increment_v;
        TR.setIdentity();
        TR =  Eigen::Translation<double,3>(pose_increment_v(0),pose_increment_v(1),0)*
              Eigen::AngleAxis<double>(pose_increment_v(2),Eigen::Vector3d::UnitZ()) ;
        T = TR*T;
//        std::cout<<"step size = "<<step_size<<std::endl;
//        std::cout<<"incr(:,"<<itr_ctr+1<<") = ["<<pose_increment_v.transpose()<<"]';\n";
//        std::cout<<"pose(:,"<<itr_ctr+2<<") = ["<<T.translation().transpose()<<" "<<T.rotation().eulerAngles(0,1,2).transpose()<<"]';\n";
        /*for(unsigned int i=0; i<nextNDT.size(); i++)
        {
            if(nextNDT[i]!=NULL)
                delete nextNDT[i];
        }
	*/
        for(unsigned int i=0; i<nextNDT.size(); i++)
        {
	    //TRANSFORM
	    Eigen::Vector3d meanC = nextNDT[i]->getMean();
	    Eigen::Matrix3d covC = nextNDT[i]->getCov();
	    meanC = TR*meanC;
	    covC = TR.rotation()*covC*TR.rotation().transpose();
	    nextNDT[i]->setMean(meanC);
	    nextNDT[i]->setCov(covC);
        }

        if(itr_ctr>0)
        {
            convergence = ((pose_increment_v.norm()) < DELTA_SCORE);
        }
        if(itr_ctr>ITR_MAX)
        {
            convergence = true;
            ret = false;
        }
        itr_ctr++;
    }
    
    //std::vector<NDTCell<PointSource>*> nextNDT = sourceNDT.pseudoTransformNDT(T);
    score_gradient_2d.setZero();
    //double score_here = derivativesNDT(nextNDT,targetNDT,score_gradient,Hessian,false);
    double score_here = derivativesNDT_2d(nextNDT,targetNDT,score_gradient_2d,H,false);
    if(score_here > score_best) 
    {
	//std::cout<<"crap iterations, best was "<<score_best<<" last was "<<score_here<<std::endl;
	T = Tbest;
    }
    for(unsigned int i=0; i<nextNDT.size(); i++)
    {
	if(nextNDT[i]!=NULL)
	    delete nextNDT[i];
    }
    
    return ret;
}

//iteratively update the score gradient and hessian (2d version)
bool NDTMatcherD2D_2D::update_gradient_hessian_local_2d(
        Eigen::MatrixXd &score_gradient,
        Eigen::MatrixXd &Hessian,
	const Eigen::Vector3d & x,
	const Eigen::Matrix3d & B,
        const double &likelihood,
        const Eigen::Matrix<double,3,3> &_Jest,
        const Eigen::Matrix<double,9,3> &_Hest,
        const Eigen::Matrix<double,3,9> &_Zest,
        const Eigen::Matrix<double,9,9> &_ZHest,
        bool computeHessian) 
{

    //vars for gradient
    Eigen::Matrix<double,3,1> _xtBJ, _xtBZBx, _Q;
    //vars for hessian
    Eigen::Matrix<double,3,3> _xtBZBJ, _xtBH, _xtBZBZBx, _xtBZhBx;
    Eigen::Matrix<double,1,3> _TMP1, _xtB;

    _xtBJ.setZero();
    _xtBZBx.setZero();
    _Q.setZero();
    _xtBZBJ.setZero();
    _xtBH.setZero();
    _xtBZBZBx.setZero();
    _xtBZhBx.setZero();
    _TMP1.setZero();
    _xtB.setZero();

    _xtB = x.transpose()*B;
    _xtBJ = _xtB*_Jest;

//    for(unsigned int i=0; i<3; i++)
//    {
        _TMP1 = _xtB*_Zest.block<3,3>(0,6)*B;
        _xtBZBx(2) = _TMP1*x;
        if(computeHessian)
        {
            _xtBZBJ.col(2) = (_TMP1*_Jest).transpose(); //-
            for(unsigned int j=0; j<3; j++)
            {
                _xtBH(2,j) = _xtB*_Hest.block<3,1>(6,j);
                _xtBZBZBx(2,j) = _TMP1*_Zest.block<3,3>(0,3*j)*B*x;
                _xtBZhBx(2,j) = _xtB*_ZHest.block<3,3>(6,3*j)*B*x;
            }
        }
//    }
    _Q = 2*_xtBJ-_xtBZBx;
    double factor = -(lfd2/2)*likelihood;
    score_gradient += _Q*factor;

    if(computeHessian)
    {
        Hessian += factor*(2*_Jest.transpose()*B*_Jest+2*_xtBH -_xtBZhBx -2*_xtBZBJ.transpose()
                           -2*_xtBZBJ +_xtBZBZBx +_xtBZBZBx.transpose() -lfd2*_Q*_Q.transpose()/2 ); // + Eigen::Matrix<double,6,6>::Identity();

    }
    return true;
}


bool NDTMatcherD2D_2D::update_gradient_hessian_2d(
    Eigen::MatrixXd &score_gradient,
    Eigen::MatrixXd &Hessian,
    const Eigen::Vector3d & x,
    const Eigen::Matrix3d & B,
    const double &likelihood,
    bool computeHessian)
{

    xtBJ.setZero();
    xtBZBx.setZero();
    Q.setZero();
    JtBJ.setZero();
    xtBZBJ.setZero();
    xtBH.setZero();
    xtBZBZBx.setZero();
    xtBZhBx.setZero();
    TMP1.setZero();
    xtB.setZero();
    
    xtB = x.transpose()*B;
    xtBJ = xtB*Jest;

    TMP1 = xtB*Zest.block<3,3>(0,6)*B;
    xtBZBx(2) = TMP1*x;
    if(computeHessian)
    {
	xtBZBJ.col(2) = (TMP1*Jest).transpose(); //-
	for(unsigned int j=0; j<3; j++)
	{
	    xtBH(2,j) = xtB*Hest.block<3,1>(6,j);
	    xtBZBZBx(2,j) = TMP1*Zest.block<3,3>(0,3*j)*B*x;
	    xtBZhBx(2,j) = xtB*ZHest.block<3,3>(6,3*j)*B*x;
	}
    }
    Q = 2*xtBJ-xtBZBx;
    double factor = -(lfd2/2)*likelihood;
    score_gradient += Q*factor;

    if(computeHessian)
    {
	Hessian += factor*(2*Jest.transpose()*B*Jest+2*xtBH -xtBZhBx -2*xtBZBJ.transpose()
		-2*xtBZBJ +xtBZBZBx +xtBZBZBx.transpose() -lfd2*Q*Q.transpose()/2 );

    }
    return true;
}

//pre-computes the derivative matrices Jest, Hest, Zest, ZHest
void NDTMatcherD2D_2D::computeDerivativesLocal_2d(Eigen::Vector3d &x, Eigen::Matrix3d C1,
                                 Eigen::Matrix<double,3,3> &_Jest,
                                 Eigen::Matrix<double,9,3> &_Hest,
                                 Eigen::Matrix<double,3,9> &_Zest,
                                 Eigen::Matrix<double,9,9> &_ZHest,
                                 bool computeHessian) 
{
    _Jest(0,2) = -x(1);
    _Jest(1,2) = x(0);

    Eigen::Matrix3d myBlock;
    //_Zest
    _Zest.block<3,3>(0,6)<<
	-2*C1(0,1), -C1(1,1) + C1(0,0),  -C1(1,2),
	-C1(1,1) + C1(0,0),    2*C1(0,1), C1(0,2),
	-C1(1,2),      C1(0,2),    0;

    if(computeHessian)
    {
	_Hest.block<3,1>(6,2)<<-x(0),-x(1),0;
	_ZHest.block<3,3>(6,6)<<
	    2*C1(1,1) - 2*C1(0,0),        -4*C1(0,1), -C1(0,2),
	    -4*C1(0,1), 2*C1(0,0) - 2*C1(1,1), -C1(1,2),
	    -C1(0,2),          -C1(1,2),    0;
    }
}

void NDTMatcherD2D_2D::computeDerivatives_2d(Eigen::Vector3d &x, Eigen::Matrix3d C1, bool computeHessian)
{

    Jest(0,2) = -x(1);
    Jest(1,2) = x(0);

    Eigen::Matrix3d myBlock;
    //_Zest
    Zest.block<3,3>(0,6)<<
	-2*C1(0,1), -C1(1,1) + C1(0,0),  -C1(1,2),
	-C1(1,1) + C1(0,0),    2*C1(0,1), C1(0,2),
	-C1(1,2),      C1(0,2),    0;

    if(computeHessian)
    {
	Hest.block<3,1>(6,2)<<-x(0),-x(1),0;
	ZHest.block<3,3>(6,6)<<
	    2*C1(1,1) - 2*C1(0,0),        -4*C1(0,1), -C1(0,2),
	    -4*C1(0,1), 2*C1(0,0) - 2*C1(1,1), -C1(1,2),
	    -C1(0,2),          -C1(1,2),    0;
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#define USE_OMP_NDT_MATCHER_D2D_2D
//compute the score gradient of a point cloud + transformation to an NDT
double NDTMatcherD2D_2D::derivativesNDT_2d(
        const std::vector<NDTCell*> &sourceNDT,
        const NDTMap &targetNDT,
        Eigen::MatrixXd &score_gradient,
        Eigen::MatrixXd &Hessian,
        bool computeHessian
    )
{

//    struct timeval tv_start, tv_end;
    double score_here = 0;
    int n_dimensions = score_gradient.rows();

//    gettimeofday(&tv_start,NULL);
    NUMBER_OF_ACTIVE_CELLS = 0;
    score_gradient.setZero();
    Hessian.setZero();

#ifdef USE_OMP_NDT_MATCHER_D2D_2D
    Eigen::MatrixXd score_gradient_omp;
    Eigen::MatrixXd score_here_omp;
    Eigen::MatrixXd Hessian_omp;

#define N_THREADS_2D 2

    //n_threads = omp_get_num_threads();
    score_gradient_omp.resize(n_dimensions,N_THREADS_2D);
    score_here_omp.resize(1,N_THREADS_2D);
    Hessian_omp.resize(n_dimensions,n_dimensions*N_THREADS_2D);

    score_gradient_omp.setZero();
    score_here_omp.setZero();
    Hessian_omp.setZero();
    //std::cout<<n_threads<<" "<<omp_get_thread_num()<<std::endl;

    #pragma omp parallel num_threads(N_THREADS_2D)
    {
        #pragma omp for
        for(unsigned int i=0; i<sourceNDT.size(); i++)
        {
            pcl::PointXYZ point;
            Eigen::Vector3d transformed;
            Eigen::Vector3d meanMoving, meanFixed;
            Eigen::Matrix3d CMoving, CFixed, CSum, Cinv, R;
            Eigen::MatrixXd score_gradient_omp_loc(n_dimensions,1);
            Eigen::MatrixXd Hessian_omp_loc(n_dimensions,n_dimensions);
            Eigen::Matrix<double,3,3> _Jest;
            Eigen::Matrix<double,9,3> _Hest;
            Eigen::Matrix<double,3,9> _Zest;
            Eigen::Matrix<double,9,9> _ZHest;
            double score_here_loc=0;
            int thread_id = omp_get_thread_num();
            NDTCell *cell;
            bool exists = false;
            double det = 0;


            score_gradient_omp_loc.setZero();
            Hessian_omp_loc.setZero();
            _Jest.setZero();
            _Jest.block<2,2>(0,0).setIdentity();
            _Hest.setZero();
            _Zest.setZero();
            _ZHest.setZero();

            meanMoving = sourceNDT[i]->getMean();
            CMoving= sourceNDT[i]->getCov();
            computeDerivativesLocal_2d(meanMoving, CMoving, _Jest, _Hest, _Zest, _ZHest, computeHessian);

            point.x = meanMoving(0);
            point.y = meanMoving(1);
            point.z = meanMoving(2);
            std::vector<NDTCell*> cells = targetNDT.getCellsForPoint(point,2); //targetNDT.getAllCells(); //
            for(unsigned int j=0; j<cells.size(); j++)
            {
                cell = cells[j];
                if(cell == NULL)
                {
                    continue;
                }
                if(cell->hasGaussian_)
                {
                    transformed = meanMoving - cell->getMean();
                    CFixed = cell->getCov();
                    CSum = (CFixed+CMoving);
                    CSum.computeInverseAndDetWithCheck(Cinv,det,exists);
                    if(!exists)
                    {
                        continue;
                    }
                    double l = (transformed).dot(Cinv*(transformed));
                    if(l*0 != 0)
                    {
                        continue;
                    }
                    //if(l > 120) continue;
                    double sh = -lfd1*(exp(-lfd2*l/2));
                    if(!update_gradient_hessian_local_2d(score_gradient_omp_loc,Hessian_omp_loc,transformed, Cinv, sh,
                                                      _Jest, _Hest, _Zest, _ZHest, computeHessian))
                    {
                        continue;
                    }
                    score_here_loc += sh;
                    cell = NULL;
                }
            }
            //score_gradient_omp.block(0,thread_id,n_dimensions,1) += score_gradient_omp_loc;
            score_gradient_omp.col(thread_id) += score_gradient_omp_loc;
            Hessian_omp.block(0,n_dimensions*thread_id,n_dimensions,n_dimensions) += Hessian_omp_loc;
            score_here_omp(0,thread_id) += score_here_loc;

        }
    } //end pragma block
    //std::cout<<"sgomp: "<<score_gradient_omp<<std::endl;
    //std::cout<<"somp: "<<score_here_omp<<std::endl;

    score_gradient = score_gradient_omp.rowwise().sum();
    score_here = score_here_omp.sum();
    if(computeHessian)
    {
        //std::cout<<"Homp: "<<Hessian_omp<<std::endl;
        for(int i=0; i<N_THREADS_2D; ++i)
        {
            Hessian += Hessian_omp.block(0,n_dimensions*i,n_dimensions,n_dimensions);
        }
    }

#else
    pcl::PointXYZ point;
    Eigen::Vector3d transformed;
    Eigen::Vector3d meanMoving, meanFixed;
    Eigen::Matrix3d CMoving, CFixed, CSum, Cinv, R;
    NDTCell *cell;
    bool exists = false;
    double det = 0;
    for(unsigned int i=0; i<sourceNDT.size(); i++)
    {
        meanMoving = sourceNDT[i]->getMean();
        CMoving= sourceNDT[i]->getCov();
        computeDerivatives_2d(meanMoving, CMoving, computeHessian);

        point.x = meanMoving(0);
        point.y = meanMoving(1);
        point.z = meanMoving(2);
        std::vector<NDTCell*> cells = targetNDT.getCellsForPoint(point,2); //targetNDT.getAllCells(); //
        for(int j=0; j<cells.size(); j++)
        {
            cell = cells[j];
            if(cell == NULL)
            {
                continue;
            }
            if(cell->hasGaussian_)
            {
                transformed = meanMoving - cell->getMean();
                CFixed = cell->getCov();
                CSum = (CFixed+CMoving);
                CSum.computeInverseAndDetWithCheck(Cinv,det,exists);
                if(!exists)
                {
                    //delete cell;
                    continue;
                }
                double l = (transformed).dot(Cinv*(transformed));
                if(l*0 != 0)
                {
                    //delete cell;
                    continue;
                }
                //if(l > 120) continue;
                double sh = -lfd1*(exp(-lfd2*l/2));
                //compute Jest, Hest, Zest, ZHest
                //update score gradient
                //std::cout<<"m1 = ["<<meanMoving.transpose()<<"]';\n m2 = ["<<cell->getMean().transpose()<<"]';\n";
                //std::cout<<"C1 = ["<<CMoving<<"];\n C2 = ["<<CFixed<<"];\n";
//		    if(!update_gradient_hessian(score_gradient, Hessian, transformed, Cinv, sh, computeHessian))
                if(!update_gradient_hessian_2d(score_gradient,Hessian,transformed, Cinv, sh, computeHessian))
                {
                    //delete cell;
                    continue;
                }
                //NUMBER_OF_ACTIVE_CELLS++;
                score_here += sh;
                //delete cell;
                cell = NULL;
            }
        }
    }
#endif

//    gettimeofday(&tv_end,NULL);
//    double time_load = (tv_end.tv_sec-tv_start.tv_sec)*1000.+(tv_end.tv_usec-tv_start.tv_usec)/1000.;
//    std::cout<<"time derivatives took is: "<<time_load<<std::endl;
    return score_here;
}

double NDTMatcherD2D_2D::lineSearch2D(
    Eigen::Matrix<double,3,1> &increment,
    std::vector<NDTCell*> &sourceNDT,
    NDTMap &targetNDT
)
{

    // default params
    double stp = 1.0; //default step
    double recoverystep = 0.01;
    double dginit = 0.0;
    double ftol = 0.11111; //epsilon 1
    double gtol = 0.99999; //epsilon 2
    double stpmax = 4.0;
    double stpmin = 0.0001;
    int maxfev = 40; //max function evaluations
    double xtol = 0.01; //window of uncertainty around the optimal step

    //my temporary variables
    std::vector<NDTCell*> sourceNDTHere;
    double score_init = 0.0;

    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> ps;
    ps.setIdentity();

    Eigen::Matrix<double,3,1> scg_here;
    Eigen::MatrixXd pincr(3,1), score_gradient_here(3,1);
    Eigen::MatrixXd pseudoH(3,3);
    Eigen::Vector3d eulerAngles;
    /////

    int info = 0;			// return code
    int infoc = 1;		// return code for subroutine cstep

    // Compute the initial gradient in the search direction and check
    // that s is a descent direction.

    //we want to maximize s, so we should minimize -s
    //score_init = scoreNDT(sourceNDT,targetNDT);

    //gradient directions are opposite for the negated function
    //score_gradient_init = -score_gradient_init;

//  cout<<"score_init "<<score_init<<endl;
//  cout<<"score_gradient_init "<<score_gradient_init.transpose()<<endl;
//  cout<<"increment "<<increment.transpose()<<endl;

    score_gradient_here.setZero();
    score_init = derivativesNDT_2d(sourceNDT,targetNDT,score_gradient_here,pseudoH,false);
    scg_here = score_gradient_here;

    //scg_here(0,0) = score_gradient_here(0,0);
    //scg_here(1,0) = score_gradient_here(1,0);
    //scg_here(2,0) = score_gradient_here(5,0);

    dginit = increment.dot(scg_here);
//  cout<<"dginit "<<dginit<<endl;

    if (dginit >= 0.0)
    {
        std::cout << "MoreThuente::cvsrch - wrong direction (dginit = " << dginit << ")" << std::endl;
        //return recoverystep; //TODO TSV -1; //
        //return -1;

        increment = -increment;
        dginit = -dginit;

        if (dginit >= 0.0)
        {
            for(unsigned int i=0; i<sourceNDTHere.size(); i++)
            {
                if(sourceNDTHere[i]!=NULL)
                    delete sourceNDTHere[i];
            }
            return recoverystep;
        }
    }
    else
    {
//     cout<<"correct direction (dginit = " << dginit << ")" << endl;
    }

    // Initialize local variables.

    bool brackt = false;		// has the soln been bracketed?
    bool stage1 = true;		// are we in stage 1?
    int nfev = 0;			// number of function evaluations
    double dgtest = ftol * dginit; // f for curvature condition
    double width = stpmax - stpmin; // interval width
    double width1 = 2 * width;	// ???

    //cout<<"dgtest "<<dgtest<<endl;
    // initial function value
    double finit = 0.0;
    finit = score_init;

    // The variables stx, fx, dgx contain the values of the step,
    // function, and directional derivative at the best step.  The
    // variables sty, fy, dgy contain the value of the step, function,
    // and derivative at the other endpoint of the interval of
    // uncertainty.  The variables stp, f, dg contain the values of the
    // step, function, and derivative at the current step.

    double stx = 0.0;
    double fx = finit;
    double dgx = dginit;
    double sty = 0.0;
    double fy = finit;
    double dgy = dginit;

    // Get the linear solve tolerance for adjustable forcing term
    //double eta_original = -1.0;
    //double eta = 0.0;
    //eta = eta_original;

    // Start of iteration.

    double stmin, stmax;
    double fm, fxm, fym, dgm, dgxm, dgym;

    while (1)
    {
        // Set the minimum and maximum steps to correspond to the present
        // interval of uncertainty.
        if (brackt)
        {
            stmin = MoreThuente::min(stx, sty);
            stmax = MoreThuente::max(stx, sty);
        }
        else
        {
            stmin = stx;
            stmax = stp + 4 * (stp - stx);
        }

        // Force the step to be within the bounds stpmax and stpmin.
        stp = MoreThuente::max(stp, stpmin);
        stp = MoreThuente::min(stp, stpmax);

        // If an unusual termination is to occur then let stp be the
        // lowest point obtained so far.

        if ((brackt && ((stp <= stmin) || (stp >= stmax))) ||
                (nfev >= maxfev - 1) || (infoc == 0) ||
                (brackt && (stmax - stmin <= xtol * stmax)))
        {
            stp = stx;
        }

        // Evaluate the function and gradient at stp
        // and compute the directional derivative.
        ///////////////////////////////////////////////////////////////////////////

        pincr = stp*increment;

        ps = Eigen::Translation<double,3>(pincr(0),pincr(1),0)*
             Eigen::AngleAxisd(pincr(2),Eigen::Vector3d::UnitZ());

        for(unsigned int i=0; i<sourceNDTHere.size(); i++)
        {
            if(sourceNDTHere[i]!=NULL)
                delete sourceNDTHere[i];
        }
        sourceNDTHere.clear();
        for(unsigned int i=0; i<sourceNDT.size(); i++)
        {
            NDTCell *cell = sourceNDT[i];
            if(cell!=NULL)
            {
                Eigen::Vector3d mean = cell->getMean();
                Eigen::Matrix3d cov = cell->getCov();
                mean = ps*mean;
                cov = ps.rotation()*cov*ps.rotation().transpose();
                NDTCell* nd = (NDTCell*)cell->copy();
                nd->setMean(mean);
                nd->setCov(cov);
                sourceNDTHere.push_back(nd);
            }
        }

        double f = 0.0;
        score_gradient_here.setZero();

        /*f = scoreNDT(sourceNDT,targetNDT,ps);
        derivativesNDT(sourceNDT,targetNDT,ps,score_gradient_here,pseudoH,false);
        std::cout<<"scg1  " <<score_gradient_here.transpose()<<std::endl;
        */

        //option 2:
        //f = scoreNDT(sourceNDTHere,targetNDT);
        f = derivativesNDT_2d(sourceNDTHere,targetNDT,score_gradient_here,pseudoH,false);
        //std::cout<<"scg2  " <<score_gradient_here.transpose()<<std::endl;


        //cout<<"incr " <<pincr.transpose()<<endl;
        //cout<<"score (f) "<<f<<endl;

        double dg = 0.0;
        scg_here = score_gradient_here;
        //scg_here(0,0) = score_gradient_here(0,0);
        //scg_here(1,0) = score_gradient_here(1,0);
        //scg_here(2,0) = score_gradient_here(5,0);
        dg = increment.dot(scg_here);


        //VALGRIND_CHECK_VALUE_IS_DEFINED(dg);
        //cout<<"dg = "<<dg<<endl;
        nfev ++;

///////////////////////////////////////////////////////////////////////////

        //cout<<"consider step "<<stp<<endl;
        // Armijo-Goldstein sufficient decrease
        double ftest1 = finit + stp * dgtest;
        //cout<<"ftest1 is "<<ftest1<<endl;

        // Test for convergence.

        if ((brackt && ((stp <= stmin) || (stp >= stmax))) || (infoc == 0))
            info = 6;			// Rounding errors

        if ((stp == stpmax) && (f <= ftest1) && (dg <= dgtest))
            info = 5;			// stp=stpmax

        if ((stp == stpmin) && ((f > ftest1) || (dg >= dgtest)))
            info = 4;			// stp=stpmin

        if (nfev >= maxfev)
            info = 3;			// max'd out on fevals

        if (brackt && (stmax-stmin <= xtol*stmax))
            info = 2;			// bracketed soln

        // RPP sufficient decrease test can be different
        bool sufficientDecreaseTest = false;
        sufficientDecreaseTest = (f <= ftest1);  // Armijo-Golstein

        //cout<<"ftest2 "<<gtol*(-dginit)<<endl;
        //cout<<"sufficientDecrease? "<<sufficientDecreaseTest<<endl;
        //cout<<"curvature ok? "<<(fabs(dg) <= gtol*(-dginit))<<endl;
        if ((sufficientDecreaseTest) && (fabs(dg) <= gtol*(-dginit)))
            info = 1;			// Success!!!!

        if (info != 0) 		// Line search is done
        {
            if (info != 1) 		// Line search failed
            {
                // RPP add
                // counter.incrementNumFailedLineSearches();

                //if (recoveryStepType == Constant)
                stp = recoverystep;

                //newgrp.computeX(oldgrp, dir, stp);

                //message = "(USING RECOVERY STEP!)";

            }
            else 			// Line search succeeded
            {
                //message = "(STEP ACCEPTED!)";
            }

            //print.printStep(nfev, stp, finit, f, message);

            // Returning the line search flag
            //cout<<"LineSearch::"<<message<<" info "<<info<<endl;
            for(unsigned int i=0; i<sourceNDTHere.size(); i++)
            {
                if(sourceNDTHere[i]!=NULL)
                    delete sourceNDTHere[i];
            }
            //std::cout<<"nfev = "<<nfev<<std::endl;
            return stp;

        } // info != 0

        // RPP add
        //counter.incrementNumIterations();

        // In the first stage we seek a step for which the modified
        // function has a nonpositive value and nonnegative derivative.

        if (stage1 && (f <= ftest1) && (dg >= MoreThuente::min(ftol, gtol) * dginit))
        {
            stage1 = false;
        }

        // A modified function is used to predict the step only if we have
        // not obtained a step for which the modified function has a
        // nonpositive function value and nonnegative derivative, and if a
        // lower function value has been obtained but the decrease is not
        // sufficient.

        if (stage1 && (f <= fx) && (f > ftest1))
        {

            // Define the modified function and derivative values.

            fm = f - stp * dgtest;
            fxm = fx - stx * dgtest;
            fym = fy - sty * dgtest;
            dgm = dg - dgtest;
            dgxm = dgx - dgtest;
            dgym = dgy - dgtest;

            // Call cstep to update the interval of uncertainty
            // and to compute the new step.

            //VALGRIND_CHECK_VALUE_IS_DEFINED(dgm);
            infoc = MoreThuente::cstep(stx,fxm,dgxm,sty,fym,dgym,stp,fm,dgm,
                                       brackt,stmin,stmax);

            // Reset the function and gradient values for f.

            fx = fxm + stx*dgtest;
            fy = fym + sty*dgtest;
            dgx = dgxm + dgtest;
            dgy = dgym + dgtest;

        }

        else
        {

            // Call cstep to update the interval of uncertainty
            // and to compute the new step.

            //VALGRIND_CHECK_VALUE_IS_DEFINED(dg);
            infoc = MoreThuente::cstep(stx,fx,dgx,sty,fy,dgy,stp,f,dg,
                                       brackt,stmin,stmax);

        }

        // Force a sufficient decrease in the size of the
        // interval of uncertainty.

        if (brackt)
        {
            if (fabs(sty - stx) >= 0.66 * width1)
                stp = stx + 0.5 * (sty - stx);
            width1 = width;
            width = fabs(sty-stx);
        }

    } // while-loop

}

int NDTMatcherD2D_2D::MoreThuente::cstep(double& stx, double& fx, double& dx,
        double& sty, double& fy, double& dy,
        double& stp, double& fp, double& dp,
        bool& brackt, double stmin, double stmax)
{
    int info = 0;

    // Check the input parameters for errors.

    if ((brackt && ((stp <= MoreThuente::min(stx, sty)) || (stp >= MoreThuente::max(stx, sty)))) ||
            (dx * (stp - stx) >= 0.0) || (stmax < stmin))
        return info;

    // Determine if the derivatives have opposite sign.

    double sgnd = dp * (dx / fabs(dx));

    // First case. A higher function value.  The minimum is
    // bracketed. If the cubic step is closer to stx than the quadratic
    // step, the cubic step is taken, else the average of the cubic and
    // quadratic steps is taken.

    bool bound;
    double theta;
    double s;
    double gamma;
    double p,q,r;
    double stpc, stpq, stpf;

    if (fp > fx)
    {
        info = 1;
        bound = 1;
        theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
        //VALGRIND_CHECK_VALUE_IS_DEFINED(theta);
        //VALGRIND_CHECK_VALUE_IS_DEFINED(dx);
        //VALGRIND_CHECK_VALUE_IS_DEFINED(dp);
        s = MoreThuente::absmax(theta, dx, dp);
        gamma = s * sqrt(((theta / s) * (theta / s)) - (dx / s) * (dp / s));
        if (stp < stx)
            gamma = -gamma;

        p = (gamma - dx) + theta;
        q = ((gamma - dx) + gamma) + dp;
        r = p / q;
        stpc = stx + r * (stp - stx);
        stpq = stx + ((dx / ((fx - fp) / (stp - stx) + dx)) / 2) * (stp - stx);
        if (fabs(stpc - stx) < fabs(stpq - stx))
            stpf = stpc;
        else
            stpf = stpc + (stpq - stpc) / 2;

        brackt = true;
    }

    // Second case. A lower function value and derivatives of opposite
    // sign. The minimum is bracketed. If the cubic step is closer to
    // stx than the quadratic (secant) step, the cubic step is taken,
    // else the quadratic step is taken.

    else if (sgnd < 0.0)
    {
        info = 2;
        bound = false;
        theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
        s = MoreThuente::absmax(theta,dx,dp);
        gamma = s * sqrt(((theta/s) * (theta/s)) - (dx / s) * (dp / s));
        if (stp > stx)
            gamma = -gamma;
        p = (gamma - dp) + theta;
        q = ((gamma - dp) + gamma) + dx;
        r = p / q;
        stpc = stp + r * (stx - stp);
        stpq = stp + (dp / (dp - dx)) * (stx - stp);
        if (fabs(stpc - stp) > fabs(stpq - stp))
            stpf = stpc;
        else
            stpf = stpq;
        brackt = true;
    }

    // Third case. A lower function value, derivatives of the same sign,
    // and the magnitude of the derivative decreases.  The cubic step is
    // only used if the cubic tends to infinity in the direction of the
    // step or if the minimum of the cubic is beyond stp. Otherwise the
    // cubic step is defined to be either stmin or stmax. The
    // quadratic (secant) step is also computed and if the minimum is
    // bracketed then the the step closest to stx is taken, else the
    // step farthest away is taken.

    else if (fabs(dp) < fabs(dx))
    {
        info = 3;
        bound = true;
        theta = 3 * (fx - fp) / (stp - stx) + dx + dp;
        s = MoreThuente::absmax(theta, dx, dp);

        // The case gamma = 0 only arises if the cubic does not tend
        // to infinity in the direction of the step.

        gamma = s * sqrt(max(0,(theta / s) * (theta / s) - (dx / s) * (dp / s)));
        if (stp > stx)
            gamma = -gamma;

        p = (gamma - dp) + theta;
        q = (gamma + (dx - dp)) + gamma;
        r = p / q;
        if ((r < 0.0) && (gamma != 0.0))
            stpc = stp + r * (stx - stp);
        else if (stp > stx)
            stpc = stmax;
        else
            stpc = stmin;

        stpq = stp + (dp/ (dp - dx)) * (stx - stp);
        if (brackt)
        {
            if (fabs(stp - stpc) < fabs(stp - stpq))
                stpf = stpc;
            else
                stpf = stpq;
        }
        else
        {
            if (fabs(stp - stpc) > fabs(stp - stpq))
                stpf = stpc;
            else
                stpf = stpq;
        }
    }

    // Fourth case. A lower function value, derivatives of the same
    // sign, and the magnitude of the derivative does not decrease. If
    // the minimum is not bracketed, the step is either stmin or
    // stmax, else the cubic step is taken.

    else
    {
        info = 4;
        bound = false;
        if (brackt)
        {
            theta = 3 * (fp - fy) / (sty - stp) + dy + dp;
            s = MoreThuente::absmax(theta, dy, dp);
            gamma = s * sqrt(((theta/s)*(theta/s)) - (dy / s) * (dp / s));
            if (stp > sty)
                gamma = -gamma;
            p = (gamma - dp) + theta;
            q = ((gamma - dp) + gamma) + dy;
            r = p / q;
            stpc = stp + r * (sty - stp);
            stpf = stpc;
        }
        else if (stp > stx)
            stpf = stmax;
        else
            stpf = stmin;
    }

    // Update the interval of uncertainty. This update does not depend
    // on the new step or the case analysis above.

    if (fp > fx)
    {
        sty = stp;
        fy = fp;
        dy = dp;
    }
    else
    {
        if (sgnd < 0.0)
        {
            sty = stx;
            fy = fx;
            dy = dx;
        }
        stx = stp;
        fx = fp;
        dx = dp;
    }

    // Compute the new step and safeguard it.

    stpf = MoreThuente::min(stmax, stpf);
    stpf = MoreThuente::max(stmin, stpf);
    stp = stpf;
    if (brackt && bound)
    {
        if (sty > stx)
            stp = min(stx + 0.66 * (sty - stx), stp);
        else
            stp = max(stx + 0.66 * (sty - stx), stp);
    }

    return info;

}

double NDTMatcherD2D_2D::MoreThuente::min(double a, double b)
{
    return (a < b ? a : b);
}

double NDTMatcherD2D_2D::MoreThuente::max(double a, double b)
{
    return (a > b ? a : b);
}

double NDTMatcherD2D_2D::MoreThuente::absmax(double a, double b, double c)
{
    a = fabs(a);
    b = fabs(b);
    c = fabs(c);

    if (a > b)
        return (a > c) ? a : c;
    else
        return (b > c) ? b : c;
}

double NDTMatcherD2D_2D::normalizeAngle(double a)
{
    //set the angle between -M_PI and M_PI
    return atan2(sin(a), cos(a));

}

}
