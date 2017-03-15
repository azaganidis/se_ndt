#include "ndt_map/ndt_cell.h"
#include "ndt_map/lazy_grid.h"
#include "ndt_map/pointcloud_utils.h"
#include "ndt_registration/ndt_matcher_d2dl.h"

#include "Eigen/Eigen"
#include <fstream>
#include <omp.h>
#include <sys/time.h>
#define NUM_MAX 50
namespace lslgeneric
{


bool NDTMatcherD2DL::match( pcl::PointCloud<pcl::PointXYZ> *target,
        pcl::PointCloud<pcl::PointXYZ> *source,
        Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T ,
        bool useInitialGuess
	)
{

    struct timeval tv_start, tv_end;
    struct timeval tv_start0, tv_end0;
    double time_load =0, time_match=0, time_combined=0;

    gettimeofday(&tv_start0,NULL);

    //initial guess
    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> Temp, Tinit;
    Tinit.setIdentity();
    if(useInitialGuess)
    {
		for(int i=0;i<NumInputs;i++)
			lslgeneric::transformPointCloudInPlace(T,source[i]);
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
			NDTMap *targetNDT[NUM_MAX];
			for(int i=0;i<NumInputs;i++)
			{
				targetNDT[i]= new NDTMap(&prototypeTarget);
				targetNDT[i]->loadPointCloud( target[i] );
				targetNDT[i]->computeNDTCells();
			}

			NDTMap *sourceNDT[NUM_MAX];
			for(int i=0;i<NumInputs;i++)
			{
				sourceNDT[i]= new NDTMap(&prototypeSource);
				sourceNDT[i]->loadPointCloud( source[i]);
				sourceNDT[i]->computeNDTCells();
			}
            gettimeofday(&tv_end,NULL);

            time_load += (tv_end.tv_sec-tv_start.tv_sec)*1000.+(tv_end.tv_usec-tv_start.tv_usec)/1000.;
            Temp.setIdentity();

            gettimeofday(&tv_start,NULL);
            ret = this->match( targetNDT, sourceNDT, Temp);
			for(int i=0;i<NumInputs;i++)
				lslgeneric::transformPointCloudInPlace(Temp,source[i]);
            gettimeofday(&tv_end,NULL);

            time_match += (tv_end.tv_sec-tv_start.tv_sec)*1000.+(tv_end.tv_usec-tv_start.tv_usec)/1000.;

            //transform moving
            T = Temp*T; //ORIGINAL
            //T = T*Temp;

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
    if(useInitialGuess)
    {
	T = T*Tinit;
    }
    gettimeofday(&tv_end0,NULL);
    time_combined = (tv_end0.tv_sec-tv_start0.tv_sec)*1000.+(tv_end0.tv_usec-tv_start0.tv_usec)/1000.;
    //std::cout<<"load: "<<time_load<<" match "<<time_match<<" combined "<<time_combined<<std::endl;
    return ret;
}

bool NDTMatcherD2DL::match( NDTMap **targetNDT,
        NDTMap **sourceNDT,
        Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T ,
        bool useInitialGuess
		)
{

    //locals
    bool convergence = false;
    //double score=0;
    double score_best = INT_MAX;
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

    //std::cout<<"pose(:,"<<1<<") = ["<<T.translation().transpose()<<" "<<T.rotation().eulerAngles(0,1,2).transpose()<<"]';\n";
    while(!convergence)
    {
        TR.setIdentity();
        Hessian.setZero();
        score_gradient.setZero();

        double score_here = derivativesNDT(nextNDT,targetNDT,score_gradient,Hessian,true);
        //derivativesNDT(nextNDT,targetNDT,score_gradient,Hessian,true);
        scg = score_gradient;
//	std::cout<<"itr "<<itr_ctr<<std::endl;
	if(score_here < score_best) 
	{
	    Tbest = T;
	    score_best = score_here;
//	    std::cout<<"best score "<<score_best<<" at "<<itr_ctr<<std::endl;
	}

//	Hessian = Hessian + score_gradient.norm()*Eigen::Matrix<double,6,6>::Identity();
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,6,6> > Sol (Hessian);
        Eigen::Matrix<double,6,1> evals = Sol.eigenvalues().real();
        double minCoeff = evals.minCoeff();
        double maxCoeff = evals.maxCoeff();
        if(minCoeff < 0)  //|| evals.minCoeff()) // < 10e-5*evals.maxCoeff()) 
        {
	    if(regularize) {
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
	    } else {
			if(score_here > score_best) 
			{
				T = Tbest;
			}
			//de-alloc nextNDT
			for(int i=0;i<NumInputs;i++)while(nextNDT[i].size()){delete nextNDT[i].back();nextNDT[i].pop_back();}
			return true;
	    }
//            std::cerr<<"regularizing\n";
        }
//	std::cout<<"s("<<itr_ctr+1<<") = "<<score_here<<";\n";
//	std::cout<<"H(:,:,"<<itr_ctr+1<<")  =  ["<< Hessian<<"];\n"<<std::endl;				  //
//	std::cout<<"grad (:,"<<itr_ctr+1<<")= ["<<score_gradient.transpose()<<"];"<<std::endl;         //
        if (score_gradient.norm()<= DELTA_SCORE)
        {
//	    std::cout<<"incr(:,"<<itr_ctr+1<<") = [0 0 0 0 0 0]';\n";
//            std::cout<<"\%gradient vanished\n";
	    if(score_here > score_best) 
	    {
//		std::cout<<"crap iterations, best was "<<score_best<<" last was "<<score_here<<std::endl;
		T = Tbest;
	    }
            //de-alloc nextNDT
			for(int i=0;i<NumInputs;i++)while(nextNDT[i].size()){delete nextNDT[i].back();nextNDT[i].pop_back();}
//	    std::cout<<"itr "<<itr_ctr<<" dScore "<< 0 <<std::endl;
            return true;
        }
//	pose_increment_v.block<3,1>(0,0) = -Hessian.block<3,3>(0,0).ldlt().solve(score_gradient.block<3,1>(0,0));
//	pose_increment_v.block<3,1>(3,0) = -Hessian.block<3,3>(3,3).ldlt().solve(score_gradient.block<3,1>(3,0));
        /*
        	Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,6,6> > Sol (Hessian);
        	Eigen::Matrix<double,6,1> evals = Sol.eigenvalues().real();
        	if(evals.minCoeff() < 0 ) {
        	    std::cout<<"\%negative hessian\n";
        	    std::cout<<"incr(:,"<<itr_ctr+1<<") = [0 0 0 0 0 0]';\n";
        	    //de-alloc nextNDT
        	    for(unsigned int i=0; i<nextNDT.size(); i++) {
        		if(nextNDT[i]!=NULL)
        		    delete nextNDT[i];
        	    }
        	    return true;
        	} else {
        	}
        */

        pose_increment_v = -Hessian.ldlt().solve(score_gradient);
        //pose_increment_v =-0.0001*scg; //
        /*	if(itr_ctr%2) {
        	    pose_increment_v(0) = 0;
        	    pose_increment_v(1) = 0;
        	    pose_increment_v(2) = 0;
        	} else {
        	pose_increment_v(3) = 0;
        	pose_increment_v(4) = 0;
        	pose_increment_v(5) = 0;
        	}
        	*/
        /*
        weights = scg.array().abs();
        weights.block(0,0,3,1) /= weights.block(0,0,3,1).maxCoeff();
        weights.block(3,0,3,1) /= weights.block(3,0,3,1).maxCoeff();
        std::cout<<"w = ["<<weights<<"];\n";
        std::cout<<"pose_increment_v= ["<<pose_increment_v.transpose()<<"];"<<std::endl;       //
        pose_increment_v = weights*pose_increment_v.array();
        //if(pose_increment_v.norm() > 1)
        //    pose_increment_v.normalize();
        //scg.normalize();
        //pose_increment_v = -0.00001*scg;
         */
        /*
        pose_increment_v /= pv;
        */

        //pose_increment_v = Hessian.jacobiSvd(Eigen::ComputeFullU|Eigen::ComputeFullV).solve(-score_gradient);
        double dginit = pose_increment_v.dot(scg);
        if(dginit > 0)
        {
//	    std::cout<<"incr(:,"<<itr_ctr+1<<") = ["<<pose_increment_v.transpose()<<"]';\n";
//	    std::cout<<"\%dg  =  "<<dginit<<std::endl;     //
//            std::cout<<"\%can't decrease in this direction any more, done \n";
            //de-alloc nextNDT
	    if(score_here > score_best) 
	    {
//		std::cout<<"crap iterations, best was "<<score_best<<" last was "<<score_here<<std::endl;
		T = Tbest;
	    }
			for(int i=0;i<NumInputs;i++)while(nextNDT[i].size()){delete nextNDT[i].back();nextNDT[i].pop_back();}
//	    std::cout<<"itr "<<itr_ctr<<" dScore "<< 0 <<std::endl;
            return true;
        }
//	std::cout<<"score("<<itr_ctr+1<<") = "<<score_here<<";\n";
        /*
           	*/


//	if(dginit > 0) {
//	    std::cout<<"pose_increment_v= ["<<pose_increment_v.transpose()<<"]"<<std::endl;       //
//	    std::cout<<"dg  =  "<<dginit<<std::endl;     //
//	}
        /*
        */
        /*
        	std::cout<<"pose_increment_v= ["<<pose_increment_v.transpose()<<"]"<<std::endl;       //
        	pose_increment_v = Hessian.jacobiSvd(Eigen::ComputeFullU|Eigen::ComputeFullV).solve(-score_gradient);
        	std::cout<<"dg    "<<pose_increment_v.dot(score_gradient)<<std::endl;     //
        	//make the initial increment reasonable...
        	//cout<<"incr_init = ["<<pose_increment_v.transpose()<<"]"<<endl;
        	double pnorm = sqrt(pose_increment_v(0)*pose_increment_v(0) + pose_increment_v(1)*pose_increment_v(1)
        			    +pose_increment_v(2)*pose_increment_v(2));
        	if(pnorm > NORM_MAX) {
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

        */
        //check direction here:

	if(step_control) {
	    step_size = lineSearchMT(pose_increment_v,nextNDT,targetNDT);
	} else {
	    step_size = 1;
	}
        pose_increment_v = step_size*pose_increment_v;
        //std::cout<<"\%iteration "<<itr_ctr<<" pose norm "<<(pose_increment_v.norm())<<" score "<<score_here<<" step "<<step_size<<std::endl;
        /*
        pose_increment_v(2) = 0;
        pose_increment_v(3) = 0;
        pose_increment_v(4) = 0;
        pose_increment_v(5) = 0;
        */


        TR.setIdentity();
        TR =  Eigen::Translation<double,3>(pose_increment_v(0),pose_increment_v(1),pose_increment_v(2))*
              Eigen::AngleAxis<double>(pose_increment_v(3),Eigen::Vector3d::UnitX()) *
              Eigen::AngleAxis<double>(pose_increment_v(4),Eigen::Vector3d::UnitY()) *
              Eigen::AngleAxis<double>(pose_increment_v(5),Eigen::Vector3d::UnitZ()) ;

        //std::cout<<"incr= ["<<pose_increment_v.transpose()<<"]"<<std::endl;
        //transform source NDT
        T = TR*T;
//	std::cout<<"incr(:,"<<itr_ctr+1<<") = ["<<pose_increment_v.transpose()<<"]';\n";
//	std::cout<<"pose(:,"<<itr_ctr+2<<") = ["<<T.translation().transpose()<<" "<<T.rotation().eulerAngles(0,1,2).transpose()<<"]';\n";

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

        if(itr_ctr>0)
        {
            convergence = ((pose_increment_v.norm()) < DELTA_SCORE);
            //convergence = ((score_gradient.norm()) < DELTA_SCORE);
        }
        if(itr_ctr>ITR_MAX)
        {
            convergence = true;
            ret = false;
        }
        itr_ctr++;
        //step_size *= alpha;
        //std::cout<<"step size "<<step_size<<std::endl;
    }
    
//    std::cout<<"itr "<<itr_ctr<<" dScore "<< pose_increment_v.norm()<<std::endl;
    //std::vector<NDTCell<PointSource>*> nextNDT = sourceNDT.pseudoTransformNDT(T);
    score_gradient.setZero();
    double score_here = derivativesNDT(nextNDT,targetNDT,score_gradient,Hessian,false);
    if(score_here > score_best) 
    {
//	std::cout<<"crap iterations, best was "<<score_best<<" last was "<<score_here<<std::endl;
	T = Tbest;
    }
	for(int i=0;i<NumInputs;i++)while(nextNDT[i].size()){delete nextNDT[i].back();nextNDT[i].pop_back();}

//    std::cout<<"incr(:,"<<itr_ctr+1<<") = [0 0 0 0 0 0]';\n";
//    std::cout<<"grad(:,"<<itr_ctr+1<<") = [0 0 0 0 0 0]';\n";
    //std::cout<<"res "<<current_resolution<<" itr "<<itr_ctr<<std::endl;

//    this->finalscore = score/NUMBER_OF_ACTIVE_CELLS;

    return ret;
}

#define USE_OMP

//compute the score gradient of a point cloud + transformation to an NDT
double NDTMatcherD2DL::derivativesNDT(
    const std::vector<NDTCell*>  *sourceNDTMany,
    const NDTMap * const * targetNDTMany,
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

#ifdef USE_OMP
    Eigen::MatrixXd score_gradient_omp;
    Eigen::MatrixXd score_here_omp;
    Eigen::MatrixXd Hessian_omp;

#define n_threads 6

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

    score_gradient = score_gradient_omp.rowwise().sum();
    score_here = score_here_omp.sum();
    if(computeHessian)
    {
        //std::cout<<"Homp: "<<Hessian_omp<<std::endl;
        for(int i=0; i<n_threads; ++i)
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

	for(int nS=0;nS<NumInputs;nS++)
	{
		std::vector<NDTCell*> sourceNDT = sourceNDTMany[nN];
		NDTMap targetNDT = targetNDTMany[nN];
    for(unsigned int i=0; i<NumInputs; i++)
    {
        meanMoving = sourceNDT[i]->getMean();
        CMoving= sourceNDT[i]->getCov();
        computeDerivatives(meanMoving, CMoving, computeHessian);

        point.x = meanMoving(0);
        point.y = meanMoving(1);
        point.z = meanMoving(2);
        std::vector<NDTCell*> cells = targetNDT.getCellsForPoint(point,n_neighbours); //targetNDT.getAllCells(); //
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
                if(!update_gradient_hessian(score_gradient,Hessian,transformed, Cinv, sh, computeHessian))
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
	}
#endif
    /*
    if(computeHessian) {

        //regularization 0.01*NUMBER_OF_ACTIVE_CELLS*
        //Hessian = Hessian + 0.1*NUMBER_OF_ACTIVE_CELLS*Eigen::Matrix<double,6,6>::Identity();
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,6,6> > Sol (Hessian);
        Eigen::Matrix<double,6,1> evals = Sol.eigenvalues().real();

        double minCoeff = evals.minCoeff();
        double maxCoeff = evals.maxCoeff();

        if(minCoeff < 0 ) { //evals.minCoeff() < 0.01*evals.maxCoeff()) {
    	Eigen::Matrix<double,6,6> evecs = Sol.eigenvectors().real();
    	double regularizer = 0.01*maxCoeff - minCoeff;
    	Eigen::Matrix<double,6,1> reg;
    	//ugly
    	reg<<regularizer,regularizer,regularizer,regularizer,regularizer,regularizer;
    	evals += reg;
    	Eigen::Matrix<double,6,6> Lam;
    	Lam = evals.asDiagonal();
    	Hessian = evecs*Lam*(evecs.transpose());
    	std::cerr<<"regularizing\n";
        }
        */
    /*
    if(minCoeff < 0.001*maxCoeff) {
    Eigen::Matrix<double,6,6> evecs = Sol.eigenvectors().real();
    for(int q=0;q<6;++q) {
        if(evals(q) < 0.001*maxCoeff) {
    	evals(q)=0.001*maxCoeff;
        } else{
    	break;
        }
    }
    Eigen::Matrix<double,6,6> Lam;
    Lam = evals.asDiagonal();
    Hessian = evecs*Lam*(evecs.transpose());
    std::cerr<<"BAD_HESSIAN\n";
    }
    }
    */

//    gettimeofday(&tv_end,NULL);
//    double time_load = (tv_end.tv_sec-tv_start.tv_sec)*1000.+(tv_end.tv_usec-tv_start.tv_usec)/1000.;
//    std::cout<<"3D::: time derivatives took is: "<<time_load<<std::endl;
    return score_here;
}


//perform line search to find the best descent rate (More&Thuente)
double NDTMatcherD2DL::lineSearchMT(
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
            for(unsigned int j=0; j<NumInputs; j++)
			{
				for(unsigned int i=0; i<sourceNDTHere[j].size(); i++)
					if(sourceNDTHere[j][i]!=NULL)
						delete sourceNDTHere[j][i];
				sourceNDTHere[j].clear();
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

		for(unsigned int j=0; j<NumInputs; j++)
		{
			for(unsigned int i=0; i<sourceNDTHere[j].size(); i++)
			{
				if(sourceNDTHere[j][i]!=NULL)
					delete sourceNDTHere[j][i];
			}
			sourceNDTHere[j].clear();
		}
		for(unsigned int j=0;j<NumInputs;j++)
		{
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
			for(unsigned int j=0; j<NumInputs; j++)
			{
				for(unsigned int i=0; i<sourceNDTHere[j].size(); i++)
				{
					if(sourceNDTHere[j][i]!=NULL)
						delete sourceNDTHere[j][i];
				}
				sourceNDTHere[j].clear();
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
}
