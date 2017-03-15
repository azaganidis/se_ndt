#include "Eigen/Eigen"
#include "ndt_map/ndt_cell.h"
#include <fstream>
#include <vector>
#include <ndt_map/lazy_grid.h>
#include <ndt_registration/ndt_matcher_p2dl.h>

namespace lslgeneric
{

//#define DO_DEBUG_PROC

void NDTMatcherP2DL::init()
{

    ///////////
    double lfc1,lfc2,lfd3;
    double integral, outlier_ratio, support_size;
    integral = 0.1;
    outlier_ratio = 0.01;
    //outlier_ratio = 0.5;
    support_size = 4; //???
    lfc1 = (1-outlier_ratio)/integral;
    lfc2 = outlier_ratio/pow(support_size,3);
    lfd3 = -log(lfc2);
    lfd1 = -log( lfc1 + lfc2 ) - lfd3;
    lfd2 = -log((-log( lfc1 * exp( -0.5 ) + lfc2 ) - lfd3 ) / lfd1);
    //d1 = 0.1;
    //d2 = 0.001;
    //cout<<lfd1<<" "<<lfd2<<endl;
    ///////////
    useSimpleDerivatives = false;
    NUMBER_OF_ACTIVE_CELLS = 0;
    ITR_MAX = 100;
    subsample_size = 0.4;

    Eigen::Vector3d tt;
    tt.setZero();
    precomputeAngleDerivatives(tt);
	current_resolution=3;
}
#define NUM_MAX 50

bool NDTMatcherP2DL::match( NDTMap **targetNDT,
        pcl::PointCloud<pcl::PointXYZ> *source,
        Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T )
{
    ///////////
    double lfc1,lfc2,lfd3;
    double integral, outlier_ratio, support_size;
    integral = 0.1;
    outlier_ratio = 0.35;
    //outlier_ratio = 0.35;
    support_size =4; //???
    lfc1 = (1-outlier_ratio)/integral;
    lfc2 = outlier_ratio/pow(support_size,3);
    lfd3 = -log(lfc2);
    lfd1 = -log( lfc1 + lfc2 ) - lfd3;
    lfd2 = -log((-log( lfc1 * exp( -0.5 ) + lfc2 ) - lfd3 ) / lfd1);
    //d1 = 0.1;
    //d2 = 0.001;
    //cout<<lfd1<<" "<<lfd2<<endl;
    ///////////
    useSimpleDerivatives = true;

    //locals
    //int ITR_MAX = 10;
    bool convergence = false;
    double score=0;
    double DELTA_SCORE = 0.0001;
//    double NORM_MAX = support_size/8, ROT_MAX = M_PI/18; //
    double NORM_MAX = 4*support_size, ROT_MAX = M_PI/4; //
    int itr_ctr = 0;
    double step_size = 1;
    Eigen::Matrix<double,6,1> pose_increment_v, pose_increment_reg_v, score_gradient, scg; //column vectors
    Eigen::Matrix<double,6,6> Hessian;
    Eigen::Matrix3d cov;
    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> TR, Tbest;
    Eigen::Vector3d transformed_vec, mean;
    bool ret = true;

    pcl::PointCloud<pcl::PointXYZ> prevCloud[NumInputs];
	for(int i=0;i<NumInputs;i++)
		prevCloud[i] = source[i];
    T.setIdentity();
    TR.setIdentity();
    Eigen::Vector3d eulerAngles = T.rotation().eulerAngles(0,1,2);

    double scoreP = 0;
    double score_best = INT_MAX;
    
    while(!convergence)
    {

        score_gradient.setZero();
        Hessian.setZero();
        //derivativesPointCloud(source,targetNDT,T,score_gradient,Hessian,true);

        TR.setIdentity();
        derivativesPointCloud(prevCloud,targetNDT,TR,score_gradient,Hessian,true);
	scg = score_gradient;

	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,6,6> > Sol (Hessian);
        Eigen::Matrix<double,6,1> evals = Sol.eigenvalues().real();
        double minCoeff = evals.minCoeff();
        double maxCoeff = evals.maxCoeff();
	if(minCoeff < 0)  //|| evals.minCoeff()) // < 10e-5*evals.maxCoeff()) 
	{
	//    std::cerr<<"Hessian near singular "<<evals.transpose()<<std::endl;
	    Eigen::Matrix<double,6,6> evecs = Sol.eigenvectors().real();
	    double regularizer = score_gradient.norm();
	    regularizer = regularizer + minCoeff > 0 ? regularizer : 0.001*maxCoeff - minCoeff;
	    //double regularizer = 0.001*maxCoeff - minCoeff;
	    Eigen::Matrix<double,6,1> reg;
	    //ugly
	    reg<<regularizer,regularizer,regularizer,regularizer,regularizer,regularizer;
	    evals += reg;
	    Eigen::Matrix<double,6,6> Lam;
	    Lam = evals.asDiagonal();
	    Hessian = evecs*Lam*(evecs.transpose());
	}
        
	pose_increment_v = Hessian.ldlt().solve(-score_gradient);

        score = scorePointCloud(prevCloud,targetNDT);
	if(score < score_best) 
	{
	    Tbest = T;
	    score_best = score;
	}
	//std::cerr<<"iteration "<<itr_ctr<<" pose norm "<<(pose_increment_v.norm())<<" score "<<score<<std::endl;
        

//step control...
#if 0
       	double pnorm = sqrt(pose_increment_v(0)*pose_increment_v(0) + pose_increment_v(1)*pose_increment_v(1)
                            +pose_increment_v(2)*pose_increment_v(2));
        if(pnorm > NORM_MAX)
        {
            pose_increment_v(0) = NORM_MAX*pose_increment_v(0)/pnorm;
            pose_increment_v(1) = NORM_MAX*pose_increment_v(1)/pnorm;
            pose_increment_v(2) = NORM_MAX*pose_increment_v(2)/pnorm;
        }
        pose_increment_v(3) = normalizeAngle(pose_increment_v(3));
        pose_increment_v(3) = (pose_increment_v(3) > ROT_MAX) ? ROT_MAX : pose_increment_v(3);
        pose_increment_v(3) = (pose_increment_v(3) < -ROT_MAX) ? -ROT_MAX : pose_increment_v(3);
        pose_increment_v(4) = normalizeAngle(pose_increment_v(4));
        pose_increment_v(4) = (pose_increment_v(4) > ROT_MAX) ? ROT_MAX : pose_increment_v(4);
        pose_increment_v(4) = (pose_increment_v(4) < -ROT_MAX) ? -ROT_MAX : pose_increment_v(4);
        pose_increment_v(5) = normalizeAngle(pose_increment_v(5));
        pose_increment_v(5) = (pose_increment_v(5) > ROT_MAX) ? ROT_MAX : pose_increment_v(5);
        pose_increment_v(5) = (pose_increment_v(5) < -ROT_MAX) ? -ROT_MAX : pose_increment_v(5);
//	cout<<"H  =  ["<<Hessian<<"]"<<endl;
//	cout<<"grad= ["<<score_gradient.transpose()<<"]"<<endl;
//	cout<<"dg    "<<pose_increment_v.dot(score_gradient)<<endl;
#endif
        TR.setIdentity();
        TR =  Eigen::Translation<double,3>(pose_increment_v(0),pose_increment_v(1),pose_increment_v(2))*
              Eigen::AngleAxis<double>(pose_increment_v(3),Eigen::Vector3d::UnitX()) *
              Eigen::AngleAxis<double>(pose_increment_v(4),Eigen::Vector3d::UnitY()) *
              Eigen::AngleAxis<double>(pose_increment_v(5),Eigen::Vector3d::UnitZ()) ;

    /*    double dginit = pose_increment_v.dot(scg);
        if (score_gradient.norm()<= DELTA_SCORE || dginit > 0)
        {
	    //std::cerr<<"Termination\n";
	    if(score > score_best) 
	    {
		T = Tbest;
	    }
	    return true;
	}
    */
        step_size = lineSearchMT(score_gradient,pose_increment_v,prevCloud,TR,targetNDT);
        if(step_size < 0)
        {
            //    cout<<"can't decrease in this direction any more, done \n";
            return true;
        }
        pose_increment_v = step_size*pose_increment_v;


//	cout<<"incr= ["<<pose_increment_v.transpose()<<"]"<<endl;
        TR.setIdentity();
        TR =  Eigen::Translation<double,3>(pose_increment_v(0),pose_increment_v(1),pose_increment_v(2))*
              Eigen::AngleAxis<double>(pose_increment_v(3),Eigen::Vector3d::UnitX()) *
              Eigen::AngleAxis<double>(pose_increment_v(4),Eigen::Vector3d::UnitY()) *
              Eigen::AngleAxis<double>(pose_increment_v(5),Eigen::Vector3d::UnitZ()) ;
        T = TR*T;

        //eulerAngles<<pose_increment_v(3),pose_increment_v(4),pose_increment_v(5);
        //eulerAngles = T.rotation().eulerAngles(0,1,2);

		for(int nIn=0;nIn<NumInputs;nIn++)
			prevCloud[nIn] = lslgeneric::transformPointCloud<pcl::PointXYZ>(T,source[nIn]);
        scoreP = score;
        score = scorePointCloud(prevCloud,targetNDT);
	if(score < score_best) 
	{
	    Tbest = T;
	    score_best = score;
	}

		std::cout<<"iteration "<<itr_ctr<<" pose norm "<<(pose_increment_v.norm())<<" score_prev "<<scoreP<<" scoreN "<<score<<std::endl;
		std::cout<<"step size "<<step_size<<std::endl;

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
	std::cout<<"res: "<<current_resolution<<" itr "<<itr_ctr<<std::endl;
	std::cout<<"T: \n t = "<<T.translation().transpose()<<std::endl;
	std::cout<<"r= \n"<<T.rotation()<<std::endl;

    T = Tbest;
    this->finalscore = score/NUMBER_OF_ACTIVE_CELLS;
    return ret;
}

bool NDTMatcherP2DL::update_score_gradient(Eigen::Matrix<double,6,1> &score_gradient,
        Eigen::Vector3d &transformed,
        Eigen::Matrix3d & Cinv,
        Eigen::Matrix<double,3,6> &_Jest)
{

    Eigen::Vector3d CxX, CxdX;
    double factor = (-lfd2*transformed.dot(Cinv*transformed)/2);

    //these conditions were copied from martin's code
    if(factor < -120)
    {
        return false;
    }
    factor = lfd2*exp(factor);
    if(factor > 1 || factor < 0 || factor*0 !=0)
    {
        return false;
    }
    factor *=lfd1;

    for(int i=0; i<6; i++)
    {
        CxdX = Cinv*_Jest.col(i);
        score_gradient(i) += transformed.dot(CxdX)*factor;
    }
    return true;

}

void NDTMatcherP2DL::update_hessian(Eigen::Matrix<double,6,6> &Hessian,
        Eigen::Vector3d &transformed,
        Eigen::Matrix3d & Cinv,
        Eigen::Matrix<double,18,6> &_Hest,
        Eigen::Matrix<double,3,6> &_Jest)
{

    Eigen::Vector3d CxX, CxdXdI, CxdXdJ, CxSecondOrder;
    CxX = Cinv*transformed;
    double factor = lfd1*lfd2*exp(-lfd2*transformed.dot(CxX)/2);
    for(int i=0; i<Hessian.rows(); i++)
    {
        for(int j=0; j<Hessian.cols(); j++)
        {

            CxdXdI = Cinv*_Jest.col(i);
            CxdXdJ = Cinv*_Jest.col(j);
            CxSecondOrder = Cinv*_Hest.block<3,1>(3*i,j);

            Hessian(i,j) += factor*(-lfd2*transformed.dot(CxdXdI)*transformed.dot(CxdXdJ) +
                                    transformed.dot(CxSecondOrder) +
                                    _Jest.col(j).dot(CxdXdI) );
        }
    }

}

void NDTMatcherP2DL::precomputeAngleDerivatives(Eigen::Vector3d &eulerAngles)
{
    if(fabsf(eulerAngles(0)) < 10e-5) eulerAngles(0) = 0;
    if(fabsf(eulerAngles(1)) < 10e-5) eulerAngles(1) = 0;
    if(fabsf(eulerAngles(2)) < 10e-5) eulerAngles(2) = 0;
    double cx,cy,cz, sx,sy,sz;
    cx = cos(eulerAngles(0));
    cy = cos(eulerAngles(1));
    cz = cos(eulerAngles(2));
    sx = sin(eulerAngles(0));
    sy = sin(eulerAngles(1));
    sz = sin(eulerAngles(2));

    jest13 << (-sx*sz+cx*sy*cz) , (-sx*cz - cx*sy*sz) , (-cx*cy) ;
    jest23 << (cx*sz+sx*sy*cz) , (cx*cz-sx*sy*sz) , (-sx*cy);
    jest04 << (-sy*cz) , sy*sz , cy;
    jest14 << sx*cy*cz , (-sx*cy*sz) , sx*sy;
    jest24 << (-cx*cy*cz) , cx*cy*sz , (-cx*sy);
    jest05 << (-cy*sz) , (-cy*cz), 0;
    jest15 << (cx*cz-sx*sy*sz) , (-cx*sz - sx*sy*cz), 0;
    jest25 << (sx*cz + cx*sy*sz) ,(cx*sy*cz - sx*sz), 0;
/*
    std::cerr<<"jest13 "<<jest13.transpose() <<std::endl;
    std::cerr<<"jest23 "<<jest23.transpose() <<std::endl;
    std::cerr<<"jest04 "<<jest04.transpose() <<std::endl;
    std::cerr<<"jest14 "<<jest14.transpose() <<std::endl;
    std::cerr<<"jest24 "<<jest24.transpose() <<std::endl;
    std::cerr<<"jest05 "<<jest05.transpose() <<std::endl;
    std::cerr<<"jest15 "<<jest15.transpose() <<std::endl;
    std::cerr<<"jest25 "<<jest25.transpose() <<std::endl;
*/
    a2 << (-cx*sz-sx*sy*cz),(-cx*cz+sx*sy*sz),sx*cy;
    a3 << (-sx*sz+cx*sy*cz),(-cx*sy*sz-sx*cz),(-cx*cy);
    b2 << (cx*cy*cz),(-cx*cy*sz),(cx*sy);
    b3 << (sx*cy*cz),(-sx*cy*sz),(sx*sy);
    c2 << (-sx*cz-cx*sy*sz),(sx*sz-cx*sy*cz),0;
    c3 << (cx*cz-sx*sy*sz),(-sx*sy*cz-cx*sz),0;
    d1 << (-cy*cz),(cy*sz),(sy);
    d2 << (-sx*sy*cz),(sx*sy*sz),(sx*cy);
    d3 << (cx*sy*cz),(-cx*sy*sz),(-cx*cy);
    e1 << (sy*sz),(sy*cz),0;
    e2 << (-sx*cy*sz),(-sx*cy*cz),0;
    e3 << (cx*cy*sz),(cx*cy*cz),0;
    f1 << (-cy*cz),(cy*sz),0;
    f2 << (-cx*sz -sx*sy*cz),(-cx*cz+sx*sy*sz),0;
    f3 << (-sx*sz+cx*sy*cz),(-cx*sy*sz-sx*cz),0;
/*
    std::cerr<<"a2 "<<a2.transpose() <<std::endl;
    std::cerr<<"a3 "<<a3.transpose() <<std::endl;
    std::cerr<<"b2 "<<b2.transpose() <<std::endl;
    std::cerr<<"b3 "<<b3.transpose() <<std::endl;
    std::cerr<<"c2 "<<c2.transpose() <<std::endl;
    std::cerr<<"c3 "<<c3.transpose() <<std::endl;
    std::cerr<<"d1 "<<d1.transpose() <<std::endl;
    std::cerr<<"d2 "<<d2.transpose() <<std::endl;
    std::cerr<<"d3 "<<d3.transpose() <<std::endl;
    std::cerr<<"e1 "<<e1.transpose() <<std::endl;
    std::cerr<<"e2 "<<e2.transpose() <<std::endl;
    std::cerr<<"e3 "<<e3.transpose() <<std::endl;
    std::cerr<<"f1 "<<f1.transpose() <<std::endl;
    std::cerr<<"f2 "<<f2.transpose() <<std::endl;
    std::cerr<<"f3 "<<f3.transpose() <<std::endl;
*/
}
//pre-computes the derivative matrices Jest, Hest
void NDTMatcherP2DL::computeDerivativesLocal(pcl::PointXYZ &pt,
        Eigen::Matrix<double,3,6> &_Jest,
        Eigen::Matrix<double,18,6> &_Hest,
        bool computeHessian)
{
	Eigen::Vector3d x;
	x(0)=pt.x;
	x(1)=pt.y;
	x(2)=pt.z;
    _Jest(0,4) = x(2);
    _Jest(0,5) = -x(1);
    _Jest(1,3) = -x(2);
    _Jest(1,5) = x(0);
    _Jest(2,3) = x(1);
    _Jest(2,4) = -x(0);

    Eigen::Matrix3d myBlock;
    if(computeHessian)
    {
	Eigen::Vector3d a,b,c,d,e,f;
	a<<0,-x(1),-x(2);
	b<<0,x(0),0;
	c<<0,0,x(0);
	d<<-x(0),0,-x(2);
	e<<0,0,x(1);
	f<<-x(0),-x(1),0;
        //Hest
        _Hest.block<3,1>(9,3) = a;
        _Hest.block<3,1>(12,3) = b;
        _Hest.block<3,1>(15,3) = c;
        _Hest.block<3,1>(9,4) = b;
        _Hest.block<3,1>(12,4) = d;
        _Hest.block<3,1>(15,4) = e;
        _Hest.block<3,1>(9,5) = c;
        _Hest.block<3,1>(12,5) = e;
        _Hest.block<3,1>(15,5) = f;
	}
}

double NDTMatcherP2DL::scorePointCloud(pcl::PointCloud<pcl::PointXYZ> *source,
        NDTMap **targetNDT)
{
    double score_here = 0;
    NUMBER_OF_ACTIVE_CELLS = 0;
	for(int jN=0;jN<NumInputs;jN++)
	{

		#pragma omp parallel num_threads(8)
		{
        #pragma omp for
		for(unsigned int i=0; i<source[jN].points.size(); i++)
		{
			Eigen::Vector3d point;
			point<<source[jN].points[i].x,source[jN].points[i].y,source[jN].points[i].z;
			double score_loc=0;
			unsigned int nact_loc=0;

			std::vector<NDTCell*> cells = targetNDT[jN]->getCellsForPoint(source[jN].points[i],current_resolution);
			for(unsigned int j=0; j<cells.size(); j++)
			{
				NDTCell *cell;
				cell = cells[j];
				Eigen::Matrix3d icov;
				Eigen::Vector3d mean;

				//{
				//    if(!targetNDT.getCellForPoint(source.points[i],cell)) {
				//        continue;
				//	}
				if(cell == NULL)
				{
					continue;
				}
				icov = cell->getInverseCov();
				mean = cell->getMean();
				double l = (point-mean).dot(icov*(point-mean));
				if(l*0 != 0) continue;

				if(l > 120) continue;

				score_loc+= (lfd1*exp(-lfd2*l/2));
				nact_loc+= 1;
			}
#pragma omp critical
			{
			score_here+=score_loc;
			NUMBER_OF_ACTIVE_CELLS+=nact_loc;
			}

		}
		}
    }
//    cout<<"here: "<<score_here<<" native: "<<score_native<<endl;
//    score_here /= NUMBER_OF_POINTS;
    return score_here;
}

//compute the score gradient of a point cloud + transformation to an NDT
void NDTMatcherP2DL::derivativesPointCloud(pcl::PointCloud<pcl::PointXYZ> *source,
        NDTMap **targetNDT,
        Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> &transform,
//	Eigen::Vector3d &eulerAngles,
        Eigen::Matrix<double,6,1> &score_gradient,
        Eigen::Matrix<double,6,6> &Hessian,
        bool computeHessian)
{

    //Eigen::Vector3d eulerAngles = transform.rotation().eulerAngles(0,1,2);
#define n_threads 6
    Eigen::MatrixXd score_gradient_omp;
    Eigen::MatrixXd Hessian_omp;

    int n_dimensions = score_gradient.rows();
    score_gradient_omp.resize(n_dimensions,n_threads);
    Hessian_omp.resize(n_dimensions,n_dimensions*n_threads);

    score_gradient_omp.setZero();
    Hessian_omp.setZero();
	unsigned long numPoints=0;
	for(int j=0;j<NumInputs;j++)
	{
		#pragma omp parallel num_threads(n_threads)
		{
        #pragma omp for
		for(unsigned int i=0; i<source[j].points.size(); i++)
		{
			NDTCell *cell;
			Eigen::Vector3d transformed;
			Eigen::Matrix3d Cinv;
			transformed<<source[j].points[i].x,source[j].points[i].y,source[j].points[i].z;
			if(!targetNDT[j]->getCellForPoint(source[j].points[i],cell)||cell==NULL)continue;

			transformed -=cell->getMean();
			Cinv = cell->getInverseCov();

            Eigen::Matrix<double,3,6> _Jest;
            Eigen::Matrix<double,18,6> _Hest;
            Eigen::Matrix<double,6,1> score_gradient_omp_loc;
            Eigen::Matrix<double,6,6> Hessian_omp_loc;
            double score_here_loc=0;
            int thread_id = omp_get_thread_num();

            score_gradient_omp_loc.setZero();
            Hessian_omp_loc.setZero();
            _Jest.setZero();
            _Jest.block<3,3>(0,0).setIdentity();
            _Hest.setZero();


			//compute Jest and Hest
			computeDerivativesLocal(source[j].points[i],_Jest,_Hest,computeHessian);

			//update score gradient
			if(!update_score_gradient(score_gradient_omp_loc, transformed, Cinv,_Jest))
			{
				continue;
			}

			//update hessian matrix
			if(computeHessian)
			{
				update_hessian(Hessian_omp_loc, transformed, Cinv,_Hest,_Jest);
			}
			cell = NULL;
            score_gradient_omp.col(thread_id) += score_gradient_omp_loc;
            Hessian_omp.block(0,n_dimensions*thread_id,n_dimensions,n_dimensions) += Hessian_omp_loc;
		}    
		}
		numPoints+=source[j].points.size();
	}
	score_gradient = score_gradient_omp.rowwise().sum();
	for(int i=0; i<n_threads; ++i)
		Hessian += Hessian_omp.block(0,n_dimensions*i,n_dimensions,n_dimensions);
	Hessian = -Hessian * (1.0 / numPoints);
	score_gradient= -score_gradient* (1.0 / numPoints);
}

//perform line search to find the best descent rate (More&Thuente)
double NDTMatcherP2DL::lineSearchMT(  Eigen::Matrix<double,6,1> &score_gradient_init,
        Eigen::Matrix<double,6,1> &increment,
        pcl::PointCloud<pcl::PointXYZ> *sourceCloud,
        Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> &globalT,
        NDTMap **targetNDT)
{

    // default params
    double stp = 4.0; //default step
    double recoverystep = 0.05;
    double dginit = 0.0;
    double ftol = 0.0001; //epsilon 1
    double gtol = 0.9999; //epsilon 2
    double stpmax = 10.0;
    double stpmin = 0.0001;
    int maxfev = 10; //max function evaluations
    double xtol = 0.01; //window of uncertainty around the optimal step

    double direction = 1.0;
    //my temporary variables
    pcl::PointCloud<pcl::PointXYZ> cloudHere[NumInputs];
	for(int i=0;i<NumInputs;i++)cloudHere[i]= sourceCloud[i];
    //cloudHere = lslgeneric::transformPointCloud(globalT,cloud);
    double score_init = 0.0;

    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> ps,ps2;
    Eigen::Matrix<double,6,1> pincr, score_gradient_here;
    Eigen::Matrix<double,6,6> pseudoH;
    Eigen::Vector3d eulerAngles;
    /////

    int info = 0;			// return code
    int infoc = 1;		// return code for subroutine cstep

    // Compute the initial gradient in the search direction and check
    // that s is a descent direction.

    //we want to maximize s, so we should minimize -s
    score_init = scorePointCloud(cloudHere,targetNDT);

    //gradient directions are opposite for the negated function
    //score_gradient_init = -score_gradient_init;

//  cout<<"score_init "<<score_init<<endl;
//  cout<<"score_gradient_init "<<score_gradient_init.transpose()<<endl;
//  cout<<"increment "<<increment.transpose()<<endl;

    dginit = increment.dot(score_gradient_init);
//  cout<<"dginit "<<dginit<<endl;

    if (dginit >= 0.0)
    {
//    cout << "MoreThuente::cvsrch - wrong direction (dginit = " << dginit << ")" << endl;
        //return recoverystep; //TODO TSV -1; //

        increment = -increment;
        dginit = -dginit;
        direction = -1;

        if (dginit >= 0.0)
        {
//    cout << "MoreThuente::cvsrch - Non-descent direction (dginit = " << dginit << ")" << endl;
            //stp = recoverystep;
            //newgrp.computeX(oldgrp, dir, stp);
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
    double eta_original = -1.0;
    double eta = 0.0;
    eta = eta_original;

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

        ps = Eigen::Translation<double,3>(pincr(0),pincr(1),pincr(2))*
             Eigen::AngleAxisd(pincr(3),Eigen::Vector3d::UnitX())*
             Eigen::AngleAxisd(pincr(4),Eigen::Vector3d::UnitY())*
             Eigen::AngleAxisd(pincr(5),Eigen::Vector3d::UnitZ());

        //ps2 = ps*globalT;
        //eulerAngles = ps2.rotation().eulerAngles(0,1,2);

        //eulerAngles<<pincr(3),pincr(4),pincr(5);
		for(int i=0;i<NumInputs;i++)
			cloudHere[i] = lslgeneric::transformPointCloud(ps,sourceCloud[i]);

        double f = 0.0;
        f = scorePointCloud(cloudHere,targetNDT);
        score_gradient_here.setZero();

        ps2.setIdentity();
        derivativesPointCloud(cloudHere,targetNDT,ps2,score_gradient_here,pseudoH,false);

        //derivativesPointCloud(cloud,ndt,ps,score_gradient_here,pseudoH,false);

        //derivativesPointCloud(cloudHere,ndt,ps2,score_gradient_here,pseudoH,false);
        //negate score gradient
        //score_gradient_here = -score_gradient_here;

        //cout<<"incr " <<pincr.transpose()<<endl;
        //cout<<"scg  " <<score_gradient_here.transpose()<<endl;
        //cout<<"score (f) "<<f<<endl;

        //VALGRIND_CHECK_VALUE_IS_DEFINED(score_gradient_here);
        //VALGRIND_CHECK_VALUE_IS_DEFINED(increment);
        double dg = 0.0;
        dg = increment.dot(score_gradient_here);


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


int NDTMatcherP2DL::MoreThuente::cstep(double& stx, double& fx, double& dx,
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

double NDTMatcherP2DL::MoreThuente::min(double a, double b)
{
    return (a < b ? a : b);
}

double NDTMatcherP2DL::MoreThuente::max(double a, double b)
{
    return (a > b ? a : b);
}

double NDTMatcherP2DL::MoreThuente::absmax(double a, double b, double c)
{
    a = fabs(a);
    b = fabs(b);
    c = fabs(c);

    if (a > b)
        return (a > c) ? a : c;
    else
        return (b > c) ? b : c;
}

double NDTMatcherP2DL::normalizeAngle(double a)
{
    //set the angle between -M_PI and M_PI
    return atan2(sin(a), cos(a));

}

}
