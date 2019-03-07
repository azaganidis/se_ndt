#include "se_ndt/ndt_matcher_d2d_se.h"
#include "ndt_map/ndt_cell.h"
#include "ndt_map/lazy_grid.h"
#include "ndt_map/pointcloud_utils.h"

#include "Eigen/Eigen"
#include <fstream>
#include <omp.h>
#include <sys/time.h>
#define NUM_MAX 50
#define n_threads 8
namespace perception_oru
{

bool NDTMatcherD2D_SE::match( NDTMap **targetNDT,
        NDTMap **sourceNDT,
        Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T ,
        bool useInitialGuess
		)
{

    //locals
    bool convergence = false;
    //double score=0;
    score_best = INT_MAX;
    //double DELTA_SCORE = 0.0005;
    //double NORM_MAX = current_resolution, ROT_MAX = M_PI/10; //
    int itr_ctr = 0;
    //double alpha = 0.95;
    double step_size = 1;
    Eigen::Matrix<double,6,1>  pose_increment_v, scg;
    Eigen::MatrixXd Hessian(6,6), score_gradient(6,1); //column vectors, pose_increment_v(6,1)

    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> TR, Tbest;
    Eigen::Vector3d transformed_vec, mean;
    bool ret = true;
    if(!useInitialGuess)
    {
        T.setIdentity();
    }
    Tbest = T;


    Eigen::Array<double,6,1> weights;
	std::vector<NDTCell*> *nextNDT=new std::vector<NDTCell*>[NumInputs]();
	for (int i=0;i<NumInputs;i++)
    	nextNDT[i]= sourceNDT[i]->pseudoTransformNDT(T);

    while(!convergence)
    {
        TR.setIdentity();
        Hessian.setZero();
        score_gradient.setZero();

        double score_here = derivativesNDT(nextNDT,targetNDT,score_gradient,Hessian,true);
        scg = score_gradient;
        if(score_here < score_best) 
        {
            Tbest = T;
            score_best = score_here;
        }
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,6,6> > Sol (Hessian);
        Eigen::Matrix<double,6,1> evals = Sol.eigenvalues().real();
        double minCoeff = evals.minCoeff();
        double maxCoeff = evals.maxCoeff();
        if(minCoeff < 0)  //|| evals.minCoeff()) // < 10e-5*evals.maxCoeff()) 
        {
            if(regularize) 
            {
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
            } else 
            {
                if(score_here > score_best) 
                {
                    T = Tbest;
                }
                //de-alloc nextNDT
                for(int i=0;i<NumInputs;i++)while(nextNDT[i].size()){delete nextNDT[i].back();nextNDT[i].pop_back();}
                delete[] nextNDT;
                return true;
            }
        }
        if (score_gradient.norm()<= DELTA_SCORE)
        {
            if(score_here > score_best) 
            {
                T = Tbest;
            }
            //de-alloc nextNDT
            for(int i=0;i<NumInputs;i++)while(nextNDT[i].size()){delete nextNDT[i].back();nextNDT[i].pop_back();}
            delete[] nextNDT;
            return true;
        }

        pose_increment_v = -Hessian.ldlt().solve(score_gradient);
        double dginit = pose_increment_v.dot(scg);
        if(dginit > 0)
        {
            if(score_here > score_best) 
            {
                T = Tbest;
            }
            for(int i=0;i<NumInputs;i++)while(nextNDT[i].size()){delete nextNDT[i].back();nextNDT[i].pop_back();}
            delete[] nextNDT;
            return true;
        }
        //check direction here:
        if(step_control) 
            step_size = lineSearchMT(pose_increment_v,nextNDT,targetNDT);
        else
            step_size = 1;
        pose_increment_v = step_size*pose_increment_v;

        TR.setIdentity();
        TR =  Eigen::Translation<double,3>(pose_increment_v(0),pose_increment_v(1),pose_increment_v(2))*
              Eigen::AngleAxis<double>(pose_increment_v(3),Eigen::Vector3d::UnitX()) *
              Eigen::AngleAxis<double>(pose_increment_v(4),Eigen::Vector3d::UnitY()) *
              Eigen::AngleAxis<double>(pose_increment_v(5),Eigen::Vector3d::UnitZ()) ;

        //transform source NDT
        T = TR*T;

        #pragma omp parallel num_threads(n_threads)
        {
        #pragma omp for
		for(unsigned int nS=0;nS<NumInputs;nS++)
			for(unsigned int i=0; i<nextNDT[nS].size(); i++)
			{
				//TRANSFORM
				Eigen::Vector3d meanC = nextNDT[nS][i]->getMean();
				Eigen::Matrix3d covC = nextNDT[nS][i]->getCov();
				meanC = TR*meanC;
				covC = TR.rotation()*covC*TR.rotation().transpose();
				nextNDT[nS][i]->setMean(meanC);
				nextNDT[nS][i]->setCov(covC);
			}
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
    score_gradient.setZero();
    double score_here = derivativesNDT(nextNDT,targetNDT,score_gradient,Hessian,true);
    if(score_here > score_best) 
    {
        T = Tbest;
    }
	for(int i=0;i<NumInputs;i++)while(nextNDT[i].size()){delete nextNDT[i].back();nextNDT[i].pop_back();}
    delete[] nextNDT;
    return ret;
}


//compute the score gradient of a point cloud + transformation to an NDT
double NDTMatcherD2D_SE::derivativesNDT(
    const std::vector<NDTCell*>  *sourceNDTMany,
    const NDTMap * const * targetNDTMany,
    Eigen::MatrixXd &score_gradient,
    Eigen::MatrixXd &Hessian,
    bool computeHessian,
    bool init_hessian_gradient
)
{
Profiler pr;
pr.start();
//    struct timeval tv_start, tv_end;
    double score_here = 0;
    int n_dimensions = score_gradient.rows();

//    gettimeofday(&tv_start,NULL);
    NUMBER_OF_ACTIVE_CELLS = 0;
    if(init_hessian_gradient){
        score_gradient.setZero();
        Hessian.setZero();
    }

    Eigen::MatrixXd score_gradient_omp;
    Eigen::MatrixXd score_here_omp;
    Eigen::MatrixXd Hessian_omp;


    //n_threads = omp_get_num_threads();
    score_gradient_omp.resize(n_dimensions,n_threads);
    score_here_omp.resize(1,n_threads);
    Hessian_omp.resize(n_dimensions,n_dimensions*n_threads);

    score_gradient_omp.setZero();
    score_here_omp.setZero();
    Hessian_omp.setZero();
    //std::cout<<n_threads<<" "<<omp_get_thread_num()<<std::endl;

	for(int nS=0;nS<NumInputs;nS++)
	{
		std::vector<NDTCell*> sourceNDT = sourceNDTMany[nS];
    #pragma omp parallel num_threads(n_threads)
    {
        #pragma omp for
        for(unsigned int i=0; i<sourceNDT.size(); i++)
        {
	    if(sourceNDT[i] == NULL) continue;
	    if(!sourceNDT[i]->hasGaussian_) continue;
            pcl::PointXYZ point;
            Eigen::Vector3d transformed;
            Eigen::Vector3d meanMoving, meanFixed;
            Eigen::Matrix3d CMoving, CFixed, CSum, Cinv, R;
            Eigen::MatrixXd score_gradient_omp_loc(n_dimensions,1);
            Eigen::MatrixXd Hessian_omp_loc(n_dimensions,n_dimensions);
            Eigen::Matrix<double,3,6> _Jest;
            Eigen::Matrix<double,18,6> _Hest;
            Eigen::Matrix<double,3,18> _Zest;
            Eigen::Matrix<double,18,18> _ZHest;
            double score_here_loc=0;
            int thread_id = omp_get_thread_num();
            NDTCell *cell;
            bool exists = false;
            double det = 0;


            score_gradient_omp_loc.setZero();
            Hessian_omp_loc.setZero();
            _Jest.setZero();
            _Jest.block<3,3>(0,0).setIdentity();
            _Hest.setZero();
            _Zest.setZero();
            _ZHest.setZero();

            meanMoving = sourceNDT[i]->getMean();
            CMoving= sourceNDT[i]->getCov();
            computeDerivativesLocal(meanMoving, CMoving, _Jest, _Hest, _Zest, _ZHest, computeHessian);

            point.x = meanMoving(0);
            point.y = meanMoving(1);
            point.z = meanMoving(2);
            std::vector<NDTCell*> cells = targetNDTMany[nS]->getCellsForPoint(point,n_neighbours); //targetNDT.getAllCells(); //
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
                    if(!update_gradient_hessian_local(score_gradient_omp_loc,Hessian_omp_loc,transformed, Cinv, sh,
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
		}
	}
    //end pragma block
    //std::cout<<"sgomp: "<<score_gradient_omp<<std::endl;
    //std::cout<<"somp: "<<score_here_omp<<std::endl;

    score_gradient += score_gradient_omp.rowwise().sum();
    score_here = score_here_omp.sum();
    if(computeHessian)
    {
        //std::cout<<"Homp: "<<Hessian_omp<<std::endl;
        for(int i=0; i<n_threads; ++i)
        {
            Hessian += Hessian_omp.block(0,n_dimensions*i,n_dimensions,n_dimensions);
        }
		HessianF=Hessian;
		score_gradientF=score_gradient;
    }
pr.check(1);
    return score_here;
}


//perform line search to find the best descent rate (More&Thuente)
double NDTMatcherD2D_SE::lineSearchMT(
    Eigen::Matrix<double,6,1> &increment,
    std::vector<NDTCell*>  *sourceNDT,
    NDTMap **targetNDT
)
{

    // default params
    double stp = 1.0; //default step
    double recoverystep = 0.1;
    double dginit = 0.0;
    double ftol = 0.11111; //epsilon 1
    double gtol = 0.99999; //epsilon 2
    double stpmax = 4.0;
    double stpmin = 0.001;
    int maxfev = 40; //max function evaluations
    double xtol = 0.01; //window of uncertainty around the optimal step

    //my temporary variables
	std::vector<NDTCell*> sourceNDTHere[NumInputs];
    double score_init = 0.0;

    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> ps;
    ps.setIdentity();

    Eigen::Matrix<double,6,1> scg_here;
    Eigen::MatrixXd pincr(6,1), score_gradient_here(6,1);
    Eigen::MatrixXd pseudoH(6,6);
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
    score_init = derivativesNDT(sourceNDT,targetNDT,score_gradient_here,pseudoH,false);
    scg_here = score_gradient_here;
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
            #pragma omp parallel num_threads(n_threads)
            {
            #pragma omp for
            for(unsigned int j=0; j<NumInputs; j++)
			{
				for(unsigned int i=0; i<sourceNDTHere[j].size(); i++)
					if(sourceNDTHere[j][i]!=NULL)
						delete sourceNDTHere[j][i];
				sourceNDTHere[j].clear();
			}
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

        ps = Eigen::Translation<double,3>(pincr(0),pincr(1),pincr(2))*
             Eigen::AngleAxisd(pincr(3),Eigen::Vector3d::UnitX())*
             Eigen::AngleAxisd(pincr(4),Eigen::Vector3d::UnitY())*
             Eigen::AngleAxisd(pincr(5),Eigen::Vector3d::UnitZ());

        #pragma omp parallel num_threads(n_threads)
        {
        #pragma omp for
		for(unsigned int j=0; j<NumInputs; j++)
		{
			for(unsigned int i=0; i<sourceNDTHere[j].size(); i++)
			{
				if(sourceNDTHere[j][i]!=NULL)
					delete sourceNDTHere[j][i];
			}
			sourceNDTHere[j].clear();
			for(unsigned int i=0; i<sourceNDT[j].size(); i++)
			{
				NDTCell *cell = sourceNDT[j][i];
				if(cell!=NULL)
				{
					Eigen::Vector3d mean = cell->getMean();
					Eigen::Matrix3d cov = cell->getCov();
					mean = ps*mean;
					cov = ps.rotation()*cov*ps.rotation().transpose();
					NDTCell* nd = (NDTCell*)cell->copy();
					nd->setMean(mean);
					nd->setCov(cov);
					sourceNDTHere[j].push_back(nd);
				}
			}
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
        f = derivativesNDT(sourceNDTHere,targetNDT,score_gradient_here,pseudoH,false);
        //std::cout<<"scg2  " <<score_gradient_here.transpose()<<std::endl;


        //cout<<"incr " <<pincr.transpose()<<endl;
        //cout<<"score (f) "<<f<<endl;

        double dg = 0.0;
        scg_here = score_gradient_here;
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
            #pragma omp parallel num_threads(n_threads)
            {
            #pragma omp for
			for(unsigned int j=0; j<NumInputs; j++)
			{
				for(unsigned int i=0; i<sourceNDTHere[j].size(); i++)
				{
					if(sourceNDTHere[j][i]!=NULL)
						delete sourceNDTHere[j][i];
				}
				sourceNDTHere[j].clear();
			}
            }
//      std::cout<<"nfev = "<<nfev<<std::endl;
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
bool NDTMatcherD2D_SE::covariance( 
    NDTMap *** targetNDTMany,
    NDTMap *** sourceNDTMany,
        Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T,
		std::vector<float>& resolutions,
        Eigen::MatrixXd &cov,
        bool inverse=false
                                                       )
{
    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> TR;
    TR.setIdentity();
    int NM=0,NM_=0;

	std::vector<NDTCell*> *nextNDT[resolutions.size()];
	std::vector<NDTCell*> *targetNDT[resolutions.size()];
	for(int nR=0;nR<resolutions.size();nR++)
	{
		nextNDT[nR]=new std::vector<NDTCell*>[NumInputs]();
		targetNDT[nR]=new std::vector<NDTCell*>[NumInputs]();
		for (int nS=0;nS<NumInputs;nS++)
		{
			nextNDT[nR][nS]= sourceNDTMany[nR][nS]->pseudoTransformNDT(T);
			targetNDT[nR][nS]= targetNDTMany[nR][nS]->pseudoTransformNDT(TR);
			NM+= nextNDT[nR][nS].size() + targetNDT[nR][nS].size();
		}
	}
    //double sigmaS = (0.03)*(0.03);
    double sigmaS = (0.03);


    Eigen::MatrixXd scg(6,1); //column vectors

    Eigen::MatrixXd Jdpdz(NM,6);

    NDTCell *cell;
    Eigen::Vector3d transformed;
    Eigen::Vector3d meanMoving, meanFixed;
    Eigen::Matrix3d CMoving, CFixed, Cinv;
    bool exists = false;
    double det = 0;
    Eigen::Matrix<double,6,1> ones;
    ones<<1,1,1,1,1,1;
    //TODO

    Eigen::Matrix3d Q;
    Jdpdz.setZero();
    Q.setZero();

    pcl::PointXYZ point;
    //now compute Jdpdz
	for(int nR=0;nR<resolutions.size();nR++)
	{
		current_resolution=resolutions.at(nR);
		derivativesNDT(nextNDT[nR],targetNDTMany[nR],scg,cov,true, false);
		for(int nS=0;nS<NumInputs;nS++)
		{
			std::vector<NDTCell*> sourceNDTN = nextNDT[nR][nS];
			std::vector<NDTCell*> targetNDTN = targetNDT[nR][nS];

			for(int i=0; i<sourceNDTN.size(); i++)
			{
				meanMoving = sourceNDTN[i]->getMean();
				point.x = meanMoving(0);
				point.y = meanMoving(1);
				point.z = meanMoving(2);

				if(!targetNDTMany[nR][nS]->getCellForPoint(point,cell)) continue;
				if(cell == NULL) continue;
				if(cell->hasGaussian_)
				{

					meanFixed = cell->getMean();
					transformed = meanMoving-meanFixed;
					CFixed = cell->getCov();
					CMoving= sourceNDTN[i]->getCov();

#ifndef DONT_COMPUTE_ZEST
					computeDerivatives(meanMoving,CMoving,false);
#endif
					(CFixed+CMoving).computeInverseAndDetWithCheck(Cinv,det,exists);
					if(!exists) continue;

					//compute Jdpdz.col(i)
					double factor = (-transformed.dot(Cinv*transformed)/2);
					//these conditions were copied from martin's code
					if(factor < -120)
					{
						continue;
					}
					factor = exp(lfd2*factor)/2;
					if(factor > 1 || factor < 0 || factor*0 !=0)
					{
						continue;
					}

					Q = -sigmaS*Cinv*Cinv;

					Eigen::Matrix<double,6,1> G, xtQJ;

					G.setZero();
					for(int q=3; q<6; q++)
					{
						G(q) =  -transformed.transpose()*Q*Zest.block<3,3>(0,3*q)*Cinv*transformed;
						G(q) = G(q) -transformed.transpose()*Cinv*Zest.block<3,3>(0,3*q)*Q*transformed;
					}

					xtQJ = transformed.transpose()*Q*Jest;

					double f1 = (transformed.transpose()*Q*transformed);
					G = G + xtQJ + (-lfd2/2)*f1*ones;
					G = G*factor*lfd1*lfd2/2;

					Jdpdz.row(NM_+i) = G.transpose();

					for(int j=0; j<targetNDTN.size(); j++)
					{
						if(targetNDTN[j]->getMean() == meanFixed)
						{

							Jdpdz.row(NM_+j+sourceNDTN.size()) = Jdpdz.row(NM_+j+sourceNDTN.size())+G.transpose();
							continue;
						}
					}

					cell = NULL;
				}
			}
			NM_+=sourceNDTN.size()+targetNDTN.size();
			for(unsigned int q=0; q<sourceNDTN.size(); q++)
			{
				delete sourceNDTN[q];
			}
			sourceNDTN.clear();
			nextNDT[nR][nS].clear();
			for(unsigned int q=0; q<targetNDTN.size(); q++)
			{
				delete targetNDTN[q];
			}
			targetNDTN .clear();
			targetNDT[nR][nS].clear();
		}
		delete[] targetNDT[nR];
		delete[] nextNDT[nR];
	}
    //cout<<Jdpdz.transpose()<<endl;

    Eigen::MatrixXd JK(6,6);
    JK = sigmaS*Jdpdz.transpose()*Jdpdz;

    //cout<<"J*J'\n"<<JK<<endl;
    //cout<<"H\n"<<cov<<endl;

    if(!inverse)
        cov = cov.inverse()*JK*cov.inverse();
    else
        cov = cov*JK.inverse()*cov;
    //cov = cov.inverse();//*fabsf(scoreNDT(sourceNDTN,targetNDT)*2/3);
    //cout<<"cov\n"<<cov<<endl;

    
    return true;
}
}
