#include "ndt_map/ndt_cell.h"
#include "ndt_map/lazy_grid.h"
#include "ndt_map/pointcloud_utils.h"
#include "ndt_registration/ndt_matcher_d2d.h"

#include "Eigen/Eigen"
#include <fstream>
#include <omp.h>
#include <sys/time.h>
namespace lslgeneric
{

//#define DO_DEBUG_PROC

void NDTMatcherD2D::init(bool _isIrregularGrid,
        bool useDefaultGridResolutions, std::vector<double> _resolutions)
{
    Jest.setZero();
    Jest.block<3,3>(0,0).setIdentity();
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
    ITR_MAX = 30;
    DELTA_SCORE = 10e-3*current_resolution;
    step_control = true;
    //should we try to regularize the hessian or just give up?
    regularize = true;
    //how many neighbours to use in the objective
    n_neighbours =2;

}

bool NDTMatcherD2D::match( pcl::PointCloud<pcl::PointXYZ>& target,
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

bool NDTMatcherD2D::match( NDTMap& targetNDT,
        NDTMap& sourceNDT,
        Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T ,
        bool useInitialGuess)
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
    std::vector<NDTCell*> nextNDT = sourceNDT.pseudoTransformNDT(T);

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
		for(unsigned int i=0; i<nextNDT.size(); i++)
		{
		    if(nextNDT[i]!=NULL)
			delete nextNDT[i];
		}
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
            for(unsigned int i=0; i<nextNDT.size(); i++)
            {
                if(nextNDT[i]!=NULL)
                    delete nextNDT[i];
            }
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
            for(unsigned int i=0; i<nextNDT.size(); i++)
            {
                if(nextNDT[i]!=NULL)
                    delete nextNDT[i];
            }
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
    for(unsigned int i=0; i<nextNDT.size(); i++)
    {
	if(nextNDT[i]!=NULL)
	    delete nextNDT[i];
    }

//    std::cout<<"incr(:,"<<itr_ctr+1<<") = [0 0 0 0 0 0]';\n";
//    std::cout<<"grad(:,"<<itr_ctr+1<<") = [0 0 0 0 0 0]';\n";
    //std::cout<<"res "<<current_resolution<<" itr "<<itr_ctr<<std::endl;

//    this->finalscore = score/NUMBER_OF_ACTIVE_CELLS;

    return ret;
}

//iteratively update the score gradient and hessian
bool NDTMatcherD2D::update_gradient_hessian_local(
    Eigen::MatrixXd &score_gradient,
    Eigen::MatrixXd &Hessian,
    const Eigen::Vector3d & x,
    const Eigen::Matrix3d & B,
    const double &likelihood,
    const Eigen::Matrix<double,3,6> &_Jest,
    const Eigen::Matrix<double,18,6> &_Hest,
    const Eigen::Matrix<double,3,18> &_Zest,
    const Eigen::Matrix<double,18,18> &_ZHest,
    bool computeHessian)
{


    //vars for gradient
    Eigen::Matrix<double,6,1> _xtBJ, _xtBZBx, _Q;
    //vars for hessian
    Eigen::Matrix<double,6,6> _xtBZBJ, _xtBH, _xtBZBZBx, _xtBZhBx;
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

    for(unsigned int i=0; i<6; i++)
    {
        _TMP1 = _xtB*_Zest.block<3,3>(0,3*i)*B;
        _xtBZBx(i) = _TMP1*x;
        if(computeHessian)
        {
            _xtBZBJ.col(i) = (_TMP1*_Jest).transpose(); //-
            for(unsigned int j=0; j<6; j++)
            {
                _xtBH(i,j) = _xtB*_Hest.block<3,1>(3*i,j);
                _xtBZBZBx(i,j) = _TMP1*_Zest.block<3,3>(0,3*j)*B*x;
                _xtBZhBx(i,j) = _xtB*_ZHest.block<3,3>(3*i,3*j)*B*x;
            }
        }
    }
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

bool NDTMatcherD2D::update_gradient_hessian(
    Eigen::MatrixXd &score_gradient,
    Eigen::MatrixXd &Hessian,
    const Eigen::Vector3d & x,
    const Eigen::Matrix3d & B,
    const double &likelihood,
    bool computeHessian)
{

    /*
        double lmax = 1;
        if(fabsf(likelihood) > lmax) {
    	return false;
        }
         */

    Eigen::MatrixXd Hloc(6,6), sg(6,1);

    Hloc.setZero();
    sg.setZero();
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

    for(unsigned int i=0; i<6; i++)
    {
        TMP1 = xtB*Zest.block<3,3>(0,3*i)*B;
        xtBZBx(i) = TMP1*x;
        if(computeHessian)
        {
            xtBZBJ.col(i) = (TMP1*Jest).transpose(); //-
            for(unsigned int j=0; j<6; j++)
            {
                xtBH(i,j) = xtB*Hest.block<3,1>(3*i,j);
                xtBZBZBx(i,j) = TMP1*Zest.block<3,3>(0,3*j)*B*x;
                xtBZhBx(i,j) = xtB*ZHest.block<3,3>(3*i,3*j)*B*x;
            }
        }
    }
    Q = 2*xtBJ-xtBZBx;
    /*
    std::cout<<"Zest = ["<<Zest<<"];\n";
    std::cout<<"Jest = ["<<Jest<<"];\n";
    std::cout<<"Q = ["<<Q.transpose()<<"]';\n";
    std::cout<<"xtBJ = ["<<xtBJ.transpose()<<"]';\n";
    std::cout<<"xtBZBx = ["<<xtBZBx.transpose()<<"]';\n";
    std::cout<<"Hest = ["<<Hest<<"];\n";
    std::cout<<"ZHest = ["<<ZHest<<"];\n";
    std::cout<<"xtBZBJ = ["<<xtBZBJ<<"];\n";
    std::cout<<"xtBH = ["<<xtBH<<"];\n";
    std::cout<<"xtBZBZBx = ["<<xtBZBZBx<<"];\n";
    std::cout<<"xtBZhBx = ["<<xtBZhBx<<"];\n";
    */
    //double factor = -(lfd2/2)*pow(1 - pow(likelihood/lmax,2),2);
    double factor = -(lfd2/2)*likelihood;
    sg = Q*factor;
//    double weight = pow(1 - pow(likelihood,2),2)
    score_gradient += sg;

    if(computeHessian)
    {
        Hloc= factor*(2*Jest.transpose()*B*Jest+2*xtBH -xtBZhBx -2*xtBZBJ.transpose()
                      -2*xtBZBJ +xtBZBZBx +xtBZBZBx.transpose() -lfd2*Q*Q.transpose()/2 ); // + Eigen::Matrix<double,6,6>::Identity();
        Hessian += Hloc;// + Eigen::Matrix<double,6,6>::Identity();
        //Hessian += Jest.transpose()*Jest;
    }

    /*
     std::cout<<"B = ["<<B<<"];\n x = ["<<x.transpose()<<"]';\n";
     std::cout<<"l = "<<likelihood<<";\nscg = ["
    <<sg.transpose()<<"]';\n H = ["
    << Hloc<<"];\n H2 = ["<<Jest.transpose()*Jest<<"];\n H3 = ["
    << sg*sg.transpose() <<"];\n";
     */
    return true;

}

//pre-computes the derivative matrices Jest, Hest, Zest, ZHest
void NDTMatcherD2D::computeDerivativesLocal(Eigen::Vector3d &x, Eigen::Matrix3d C1,
        Eigen::Matrix<double,3,6> &_Jest,
        Eigen::Matrix<double,18,6> &_Hest,
        Eigen::Matrix<double,3,18> &_Zest,
        Eigen::Matrix<double,18,18> &_ZHest,
        bool computeHessian)
{

    _Jest(0,4) = x(2);
    _Jest(0,5) = -x(1);
    _Jest(1,3) = -x(2);
    _Jest(1,5) = x(0);
    _Jest(2,3) = x(1);
    _Jest(2,4) = -x(0);

    Eigen::Matrix3d myBlock;
    //_Zest
    myBlock<<
           0,       -C1(0,2),      C1(0,1),
                    -C1(0,2),     -2*C1(1,2), -C1(2,2) + C1(1,1),
                    C1(0,1), -C1(2,2) + C1(1,1),    2*C1(1,2);
    _Zest.block<3,3>(0,9) = myBlock;
    myBlock<<
           2*C1(0,2), C1(1,2), -C1(0,0) + C1(2,2),
             C1(1,2),    0,       -C1(0,1),
             -C1(0,0) + C1(2,2),  -C1(0,1),     -2*C1(0,2);
    _Zest.block<3,3>(0,12) = myBlock;
    myBlock<<
           -2*C1(0,1), -C1(1,1) + C1(0,0),  -C1(1,2),
           -C1(1,1) + C1(0,0),    2*C1(0,1), C1(0,2),
           -C1(1,2),      C1(0,2),    0;
    _Zest.block<3,3>(0,15) = myBlock;

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

        //_ZHest
        myBlock<<
               0,          -C1(0,1),          -C1(0,2),
                           -C1(0,1), 2*C1(2,2) - 2*C1(1,1),        -4*C1(1,2),
                           -C1(0,2),        -4*C1(1,2), 2*C1(1,1) - 2*C1(2,2);
        _ZHest.block<3,3>(9,9) =   myBlock;

        myBlock<<
               0, C1(0,0) - C1(2,2),    C1(1,2),
                  C1(0,0) - C1(2,2),     2*C1(0,1),  2*C1(0,2),
                  C1(1,2),     2*C1(0,2), -2*C1(0,1);
        _ZHest.block<3,3>(9,12) = myBlock;

        myBlock<<
               0,    C1(1,2), C1(0,0) - C1(1,1),
                     C1(1,2), -2*C1(0,2),     2*C1(0,1),
                     C1(0,0) - C1(1,1),  2*C1(0,1),     2*C1(0,2);
        _ZHest.block<3,3>(9,15) = myBlock;

        myBlock<<
               2*C1(2,2) - 2*C1(0,0), -C1(0,1),        -4*C1(0,2),
                 -C1(0,1),    0,          -C1(1,2),
                 -4*C1(0,2), -C1(1,2), 2*C1(0,0) - 2*C1(2,2);
        _ZHest.block<3,3>(12,12) = myBlock;

        myBlock<<
               -2*C1(1,2),       C1(0,2),     2*C1(0,1),
               C1(0,2),         0, C1(1,1) - C1(0,0),
               2*C1(0,1), C1(1,1) - C1(0,0),     2*C1(1,2);
        _ZHest.block<3,3>(12,15) = myBlock;

        myBlock<<
               2*C1(1,1) - 2*C1(0,0),        -4*C1(0,1), -C1(0,2),
                 -4*C1(0,1), 2*C1(0,0) - 2*C1(1,1), -C1(1,2),
                 -C1(0,2),          -C1(1,2),    0;
        _ZHest.block<3,3>(15,15)= myBlock;

        _ZHest.block<3,3>(12,9) =    _ZHest.block<3,3>(9,12);
        _ZHest.block<3,3>(15,9) =    _ZHest.block<3,3>(9,15);
        _ZHest.block<3,3>(15,12)=    _ZHest.block<3,3>(12,15);
    }
}

void NDTMatcherD2D::computeDerivatives(Eigen::Vector3d &x, Eigen::Matrix3d C1, bool computeHessian)
{

    Jest(0,4) = x(2);
    Jest(0,5) = -x(1);
    Jest(1,3) = -x(2);
    Jest(1,5) = x(0);
    Jest(2,3) = x(1);
    Jest(2,4) = -x(0);

    Eigen::Matrix3d myBlock;
    //Zest
    myBlock<<
           0,       -C1(0,2),      C1(0,1),
                    -C1(0,2),     -2*C1(1,2), -C1(2,2) + C1(1,1),
                    C1(0,1), -C1(2,2) + C1(1,1),    2*C1(1,2);
    Zest.block<3,3>(0,9) = myBlock;
    myBlock<<
           2*C1(0,2), C1(1,2), -C1(0,0) + C1(2,2),
             C1(1,2),    0,       -C1(0,1),
             -C1(0,0) + C1(2,2),  -C1(0,1),     -2*C1(0,2);
    Zest.block<3,3>(0,12) = myBlock;
    myBlock<<
           -2*C1(0,1), -C1(1,1) + C1(0,0),  -C1(1,2),
           -C1(1,1) + C1(0,0),    2*C1(0,1), C1(0,2),
           -C1(1,2),      C1(0,2),    0;
    Zest.block<3,3>(0,15) = myBlock;

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
        Hest.block<3,1>(9,3) = a;
        Hest.block<3,1>(12,3) = b;
        Hest.block<3,1>(15,3) = c;
        Hest.block<3,1>(9,4) = b;
        Hest.block<3,1>(12,4) = d;
        Hest.block<3,1>(15,4) = e;
        Hest.block<3,1>(9,5) = c;
        Hest.block<3,1>(12,5) = e;
        Hest.block<3,1>(15,5) = f;

        //ZHest
        myBlock<<
               0,          -C1(0,1),          -C1(0,2),
                           -C1(0,1), 2*C1(2,2) - 2*C1(1,1),        -4*C1(1,2),
                           -C1(0,2),        -4*C1(1,2), 2*C1(1,1) - 2*C1(2,2);
        ZHest.block<3,3>(9,9) =   myBlock;

        myBlock<<
               0, C1(0,0) - C1(2,2),    C1(1,2),
                  C1(0,0) - C1(2,2),     2*C1(0,1),  2*C1(0,2),
                  C1(1,2),     2*C1(0,2), -2*C1(0,1);
        ZHest.block<3,3>(9,12) = myBlock;

        myBlock<<
               0,    C1(1,2), C1(0,0) - C1(1,1),
                     C1(1,2), -2*C1(0,2),     2*C1(0,1),
                     C1(0,0) - C1(1,1),  2*C1(0,1),     2*C1(0,2);
        ZHest.block<3,3>(9,15) = myBlock;

        myBlock<<
               2*C1(2,2) - 2*C1(0,0), -C1(0,1),        -4*C1(0,2),
                 -C1(0,1),    0,          -C1(1,2),
                 -4*C1(0,2), -C1(1,2), 2*C1(0,0) - 2*C1(2,2);
        ZHest.block<3,3>(12,12) = myBlock;

        myBlock<<
               -2*C1(1,2),       C1(0,2),     2*C1(0,1),
               C1(0,2),         0, C1(1,1) - C1(0,0),
               2*C1(0,1), C1(1,1) - C1(0,0),     2*C1(1,2);
        ZHest.block<3,3>(12,15) = myBlock;

        myBlock<<
               2*C1(1,1) - 2*C1(0,0),        -4*C1(0,1), -C1(0,2),
                 -4*C1(0,1), 2*C1(0,0) - 2*C1(1,1), -C1(1,2),
                 -C1(0,2),          -C1(1,2),    0;
        ZHest.block<3,3>(15,15)= myBlock;

        ZHest.block<3,3>(12,9) =    ZHest.block<3,3>(9,12);
        ZHest.block<3,3>(15,9) =    ZHest.block<3,3>(9,15);
        ZHest.block<3,3>(15,12)=    ZHest.block<3,3>(12,15);
    }
}

double NDTMatcherD2D::scoreNDT(std::vector<NDTCell*> &sourceNDT, NDTMap &targetNDT)
{
    NUMBER_OF_ACTIVE_CELLS = 0;
    NDTCell *cell;
    Eigen::Vector3d transformed, eps;
    Eigen::Vector3d meanMoving, meanFixed;
    Eigen::Matrix3d CMoving, CFixed, CSum, Cinv, R;
//	eps<<0.01,0.01,0.01;
    bool exists = false;
    double det = 0;
    double score_here = 0;

    pcl::PointXYZ point;
    for(unsigned int i=0; i<sourceNDT.size(); i++)
    {
        meanMoving = sourceNDT[i]->getMean();
        CMoving= sourceNDT[i]->getCov();
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
                double l = (transformed).dot(Cinv*(transformed));// + (eps).dot(Cinv*(eps));
                if(l*0 != 0)
                {
                    continue;
                }
                //if(l > 120) continue;
                double sh = -lfd1*(exp(-lfd2*l/2));
                NUMBER_OF_ACTIVE_CELLS++;
                score_here += sh;
                cell = NULL;
            }
            else
            {
                //delete cell;
                //std::cout << "had not gaussian!" << std::endl;
            }
        }
    }

    return score_here;
}

double NDTMatcherD2D::scoreNDT_OM(NDTMap &sourceNDT, NDTMap &targetNDT)
{
    NUMBER_OF_ACTIVE_CELLS = 0;
    Eigen::Vector3d transformed;
    Eigen::Vector3d meanMoving, meanFixed;
    Eigen::Matrix3d CMoving, CFixed, CSum, Cinv, R;
    bool exists = false;
    double det = 0;
    double score_here = 0;
    double importance_free = 0.1;

    pcl::PointXYZ point;
    LazyGrid* lz = dynamic_cast<LazyGrid*> (sourceNDT.getMyIndex());
    if(lz==NULL) return INT_MAX;
    std::vector<NDTCell*>::iterator it = lz->begin();

    //std::vector<NDTCell<PointSource>*> source = sourceNDT.getAllInitializedCells();

    //for(unsigned int i=0; i<source.size(); i++) {
    while(it!=lz->end()) {
        if(*it == NULL) {
            it++;
            continue;
        }
        NDTCell* source =(*it);
        if(source==NULL) {
            it++;
            continue;
        }
        source->consistency_score = 0;
        NDTCell* cell=NULL;
        point = source->getCenter();
        //SWITCHME
        std::vector<NDTCell*> all_cells = targetNDT.getCellsForPoint(point,2,false);

        for(unsigned int j=0; j<all_cells.size(); j++) {
            cell = all_cells[j];
            if(cell == NULL) {
                it++;
                continue;
            }
            double o1,o2;
            o1 = source->getOccupancyRescaled();
            o2 = cell->getOccupancyRescaled();
            //std::cout<<"o1 "<<o1<<" o2 "<<o2<<std::endl;

            if(source->hasGaussian_) {
                meanMoving = source->getMean();
                CMoving= source->getCov();
                //point.x = meanMoving(0); point.y = meanMoving(1); point.z = meanMoving(2);
                //targetNDT.getCellForPoint(point,cell,false);
                if(cell->hasGaussian_) {
                    //both have Gaussians
                    transformed = meanMoving - cell->getMean();
                    CFixed = cell->getCov();
                    CSum = (CFixed+CMoving);
                    CSum.computeInverseAndDetWithCheck(Cinv,det,exists);
                    if(!exists)
                    {
                        it++;
                        continue;
                    }
                    double l = (transformed).dot(Cinv*(transformed));
                    if(l*0 != 0) {
                        it++;
                        continue;
                    }
                    double sh = -lfd1*(exp(-lfd2*l/2));
                    //std::cout<<"sh: "<<sh<<std::endl;
                    score_here += o1*o2*sh;
                    //SWITCHME
                    //source->consistency_score += o1*o2*sh;
                    source->consistency_score += -o1*o2*sh;
                }
            }
            //add the part of score updates based on occupancy alone
            pcl::PointXYZ cen1 = cell->getCenter();
            pcl::PointXYZ cen2 = source->getCenter();
            double cell_dist = std::sqrt(pow(cen1.x-cen2.x,2)+pow(cen1.y-cen2.y,2)+pow(cen1.z-cen2.z,2));
            double lambda = 1 - cell_dist/(cell->getDiagonal());
            lambda = lambda < 0 ? 0 : importance_free*lambda;
            //std::cout<<"lam "<<lambda<<" comp "<<-lambda*((1-o1)*(1-o2) - o1*(1-o2) -o2*(1-o1))<<std::endl;
            score_here += -lambda*((1-o1)*(1-o2) - o1*(1-o2) -o2*(1-o1));

            //SWITCHME
            //source->consistency_score += -lambda*((1-o1)*(1-o2) - o1*(1-o2) -o2*(1-o1));
            double ll = cell->isInside(cen2) ? 1 : 0;
            //std::cout<<"j "<<j<<" sc "<<source->consistency_score<<" ll "<<ll<<" incr "<<ll*((1-o1)*(1-o2) - o1*(1-o2) -o2*(1-o1))<<std::endl;
            source->consistency_score += ll*((1-o1)*(1-o2) - o1*(1-o2) -o2*(1-o1));
        }
        it++;
    }
//cleanup
    /*
    for(unsigned int i=0; i<source.size(); i++) {
        if(source[i]==NULL) {
    	continue;
        }
        delete source[i];
        source[i] = NULL;
    }
    */
    return score_here;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Added by Jari - Uses 1-p as a score instead of -p
///
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
double NDTMatcherD2D::scoreNDTPositive(std::vector<NDTCell*> &sourceNDT, NDTMap &targetNDT,
        Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T)
{

    NUMBER_OF_ACTIVE_CELLS = 0;
    double score_here = 0;
    double det = 0;
    bool exists = false;
    NDTCell *cell;
    Eigen::Matrix3d covCombined, icov, R;
    Eigen::Vector3d meanFixed;
    Eigen::Vector3d meanMoving;
    pcl::PointXYZ point;

    R = T.rotation();
    for(unsigned int i=0; i<sourceNDT.size(); i++)
    {
        meanMoving = T*sourceNDT[i]->getMean();
        point.x = meanMoving(0);
        point.y = meanMoving(1);
        point.z = meanMoving(2);

        if(!targetNDT.getCellForPoint(point,cell))
        {
            score_here += 0.1;
            continue;
        }

        if(cell == NULL)
        {
            score_here += 0.1;
            continue;
        }
        if(cell->hasGaussian_)
        {
            meanFixed = cell->getMean();
            covCombined = cell->getCov() + R.transpose()*sourceNDT[i]->getCov()*R;
            covCombined.computeInverseAndDetWithCheck(icov,det,exists);
            if(!exists)
            {
                score_here+=0.1;
                continue;
            }
            double l = (meanMoving-meanFixed).dot(icov*(meanMoving-meanFixed));
            if(l*0 != 0)
            {
                score_here+=0.1;
                continue;
            }
            if(l > 120)
            {
                score_here+=0.1;
                continue;
            }

            double sh = lfd1*(exp(-lfd2*l/2));

            if(fabsf(sh) > 1e-10)
            {
                NUMBER_OF_ACTIVE_CELLS++;
            }
            score_here += (1.0 - sh);
            //fprintf(stderr," '[i]=%lf' ", (1.0-sh));
        }
        else
        {
            score_here+=0.1;
        }

    }
    return score_here;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



#define USE_OMP

//compute the score gradient of a point cloud + transformation to an NDT
double NDTMatcherD2D::derivativesNDT(
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
            std::vector<NDTCell*> cells = targetNDT.getCellsForPoint(point,n_neighbours); //targetNDT.getAllCells(); //
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
    } //end pragma block
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
    for(unsigned int i=0; i<sourceNDT.size(); i++)
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
double NDTMatcherD2D::lineSearchMT(
    Eigen::Matrix<double,6,1> &increment,
    std::vector<NDTCell*> &sourceNDT,
    NDTMap &targetNDT
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
    std::vector<NDTCell*> sourceNDTHere;
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

        ps = Eigen::Translation<double,3>(pincr(0),pincr(1),pincr(2))*
             Eigen::AngleAxisd(pincr(3),Eigen::Vector3d::UnitX())*
             Eigen::AngleAxisd(pincr(4),Eigen::Vector3d::UnitY())*
             Eigen::AngleAxisd(pincr(5),Eigen::Vector3d::UnitZ());

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
            for(unsigned int i=0; i<sourceNDTHere.size(); i++)
            {
                if(sourceNDTHere[i]!=NULL)
                    delete sourceNDTHere[i];
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

int NDTMatcherD2D::MoreThuente::cstep(double& stx, double& fx, double& dx,
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

double NDTMatcherD2D::MoreThuente::min(double a, double b)
{
    return (a < b ? a : b);
}

double NDTMatcherD2D::MoreThuente::max(double a, double b)
{
    return (a > b ? a : b);
}

double NDTMatcherD2D::MoreThuente::absmax(double a, double b, double c)
{
    a = fabs(a);
    b = fabs(b);
    c = fabs(c);

    if (a > b)
        return (a > c) ? a : c;
    else
        return (b > c) ? b : c;
}

double NDTMatcherD2D::normalizeAngle(double a)
{
    //set the angle between -M_PI and M_PI
    return atan2(sin(a), cos(a));

}

/*
void NDTMatcherD2D<PointSource,PointTarget>::generateScoreDebug(const char* out, pcl::PointCloud<pcl::PointXYZ>& fixed, pcl::PointCloud<pcl::PointXYZ>& moving) {

    std::ofstream lg(out,std::ios_base::out);
    int N_LINEAR = 100;
    int N_ROT	 = 100;
    //lfd1 = 1;
    //lfd2 = 1;

    cout<<"generating scores...\n";
    for(current_resolution = 4; current_resolution>=0.5; current_resolution/=2) {
	cout<<"res "<<current_resolution<<endl;
	double lfc1,lfc2,lfd3;
	double integral, outlier_ratio, support_size;
	integral = 0.1;
	outlier_ratio = 0.3;
	support_size = current_resolution;
	lfc1 = (1-outlier_ratio)/integral;
	lfc2 = outlier_ratio/pow(support_size,3);
	lfd3 = -log(lfc2);
	lfd1 = -(-log( lfc1 + lfc2 ) - lfd3);
	lfd2 = -log((-log( lfc1 * exp( -0.5 ) + lfc2 ) - lfd3 ) / -lfd1);

	lfd1 = 1;//lfd1;
	lfd2 = 0.05;//0.8*lfd2;

	double lmin=-2, lmax=2, rmin=-M_PI/2, rmax=M_PI/2;
	double lstep = (lmax-lmin)/(N_LINEAR-1);
	double rstep = (rmax-rmin)/(N_ROT-1);
	Eigen::MatrixXd S(6,N_LINEAR);
	Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> T;
	std::vector<NDTCell*> nextNDT;

	LazzyGrid prototype(current_resolution);
	NDTMap ndt( &prototype );
	ndt.loadPointCloud( fixed );
	ndt.computeNDTCells();

	double z=0, r=0, pitch=0, yaw=0;
	int k=0;
	for(double x=lmin; x<lmax; x+=lstep) {
	    T = Eigen::Translation<double,3>(x,0,0);
	    pcl::PointCloud<pcl::PointXYZ> cloud = moving;
	    lslgeneric::transformPointCloudInPlace(T,cloud);
	    NDTMap mov( &prototype );
	    mov.loadPointCloud( cloud );
	    mov.computeNDTCells();
	    T.setIdentity();
	    nextNDT = mov.pseudoTransformNDT(T);

	    S(0,k) = scoreNDT(nextNDT,ndt);
	    for(unsigned int i=0; i<nextNDT.size(); i++) {
		if(nextNDT[i]!=NULL) delete nextNDT[i];
	    }
	    k++;
	}
	k=0;
	for(double x=lmin; x<lmax; x+=lstep) {
	    T = Eigen::Translation<double,3>(0,x,0);
	    pcl::PointCloud<pcl::PointXYZ> cloud = moving;
	    lslgeneric::transformPointCloudInPlace(T,cloud);
	    NDTMap mov( &prototype );
	    mov.loadPointCloud( cloud );
	    mov.computeNDTCells();
	    T.setIdentity();
	    nextNDT = mov.pseudoTransformNDT(T);

	    S(1,k) = scoreNDT(nextNDT,ndt);
	    for(unsigned int i=0; i<nextNDT.size(); i++) {
		if(nextNDT[i]!=NULL) delete nextNDT[i];
	    }
	    k++;
	}
	k=0;
	for(double x=lmin; x<lmax; x+=lstep) {
	    T = Eigen::Translation<double,3>(0,0,x);
	    pcl::PointCloud<pcl::PointXYZ> cloud = moving;
	    lslgeneric::transformPointCloudInPlace(T,cloud);
	    NDTMap mov( &prototype );
	    mov.loadPointCloud( cloud );
	    mov.computeNDTCells();
	    T.setIdentity();
	    nextNDT = mov.pseudoTransformNDT(T);

	    S(2,k) = scoreNDT(nextNDT,ndt);
	    for(unsigned int i=0; i<nextNDT.size(); i++) {
		if(nextNDT[i]!=NULL) delete nextNDT[i];
	    }
	    k++;
	}

	k=0;
	for(double r=rmin; r<rmax; r+=rstep) {
	    T = Eigen::AngleAxis<double>(r,Eigen::Vector3d::UnitX()) *
		Eigen::AngleAxis<double>(0,Eigen::Vector3d::UnitY()) *
		Eigen::AngleAxis<double>(0,Eigen::Vector3d::UnitZ()) ;
	    pcl::PointCloud<pcl::PointXYZ> cloud = moving;
	    lslgeneric::transformPointCloudInPlace(T,cloud);
	    NDTMap mov( &prototype );
	    mov.loadPointCloud( cloud );
	    mov.computeNDTCells();
	    T.setIdentity();
	    nextNDT = mov.pseudoTransformNDT(T);

	    S(3,k) = scoreNDT(nextNDT,ndt);
	    for(unsigned int i=0; i<nextNDT.size(); i++) {
		if(nextNDT[i]!=NULL) delete nextNDT[i];
	    }
	    k++;
	}
	k=0;
	for(double r=rmin; r<rmax; r+=rstep) {
	    T = Eigen::AngleAxis<double>(0,Eigen::Vector3d::UnitX()) *
		Eigen::AngleAxis<double>(r,Eigen::Vector3d::UnitY()) *
		Eigen::AngleAxis<double>(0,Eigen::Vector3d::UnitZ()) ;
	    pcl::PointCloud<pcl::PointXYZ> cloud = moving;
	    lslgeneric::transformPointCloudInPlace(T,cloud);
	    NDTMap mov( &prototype );
	    mov.loadPointCloud( cloud );
	    mov.computeNDTCells();
	    T.setIdentity();
	    nextNDT = mov.pseudoTransformNDT(T);

	    S(4,k) = scoreNDT(nextNDT,ndt);
	    for(unsigned int i=0; i<nextNDT.size(); i++) {
		if(nextNDT[i]!=NULL) delete nextNDT[i];
	    }
	    k++;
	}
	k=0;
	for(double r=rmin; r<rmax; r+=rstep) {
	    T = Eigen::AngleAxis<double>(0,Eigen::Vector3d::UnitX()) *
		Eigen::AngleAxis<double>(0,Eigen::Vector3d::UnitY()) *
		Eigen::AngleAxis<double>(r,Eigen::Vector3d::UnitZ()) ;
	    pcl::PointCloud<pcl::PointXYZ> cloud = moving;
	    lslgeneric::transformPointCloudInPlace(T,cloud);
	    NDTMap mov( &prototype );
	    mov.loadPointCloud( cloud );
	    mov.computeNDTCells();
	    T.setIdentity();
	    nextNDT = mov.pseudoTransformNDT(T);

	    S(5,k) = scoreNDT(nextNDT,ndt);
	    for(unsigned int i=0; i<nextNDT.size(); i++) {
		if(nextNDT[i]!=NULL) delete nextNDT[i];
	    }
	    k++;
	}

	lg<<"Sf2f"<<(int)current_resolution<<" = ["<<S<<"];\n";
    }
    lg.close();

}*/

bool NDTMatcherD2D::covariance( NDTMap& targetNDT,
        NDTMap& sourceNDT,
        Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T,
        Eigen::MatrixXd &cov
                                                       )
{

    double sigmaS = (0.03)*(0.03);
    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> TR;
    TR.setIdentity();

    std::vector<NDTCell*> sourceNDTN = sourceNDT.pseudoTransformNDT(T);
    std::vector<NDTCell*> targetNDTN = targetNDT.pseudoTransformNDT(T);

    Eigen::MatrixXd scg(6,1); //column vectors
    int NM = sourceNDTN.size() + targetNDTN.size();

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
    derivativesNDT(sourceNDTN,targetNDT,scg,cov,true);

    Eigen::Matrix3d Q;
    Jdpdz.setZero();
    Q.setZero();

    pcl::PointXYZ point;
    //now compute Jdpdz
    for(int i=0; i<sourceNDTN.size(); i++)
    {
        meanMoving = sourceNDTN[i]->getMean();
        point.x = meanMoving(0);
        point.y = meanMoving(1);
        point.z = meanMoving(2);

        if(!targetNDT.getCellForPoint(point,cell))
        {
            continue;
        }
        if(cell == NULL)
        {
            continue;
        }
        if(cell->hasGaussian_)
        {

            meanFixed = cell->getMean();
            transformed = meanMoving-meanFixed;
            CFixed = cell->getCov();
            CMoving= sourceNDTN[i]->getCov();

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

            Jdpdz.row(i) = G.transpose();

            for(int j=0; j<targetNDTN.size(); j++)
            {
                if(targetNDTN[j]->getMean() == meanFixed)
                {

                    Jdpdz.row(j+sourceNDTN.size()) = Jdpdz.row(j+sourceNDTN.size())+G.transpose();
                    continue;
                }
            }

            cell = NULL;
        }
    }


    //cout<<Jdpdz.transpose()<<endl;

    Eigen::MatrixXd JK(6,6);
    JK = sigmaS*Jdpdz.transpose()*Jdpdz;

    //cout<<"J*J'\n"<<JK<<endl;
    //cout<<"H\n"<<cov<<endl;

    cov = cov.inverse()*JK*cov.inverse();
    //cov = cov.inverse();//*fabsf(scoreNDT(sourceNDTN,targetNDT)*2/3);
    //cout<<"cov\n"<<cov<<endl;

    for(unsigned int q=0; q<sourceNDTN.size(); q++)
    {
        delete sourceNDTN[q];
    }
    sourceNDTN.clear();

    return true;
}
bool NDTMatcherD2D::covariance( pcl::PointCloud<pcl::PointXYZ>& target,
        pcl::PointCloud<pcl::PointXYZ>& source,
        Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T,
        Eigen::MatrixXd &cov
                                                       )
{

    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> TR;
    TR.setIdentity();

    pcl::PointCloud<pcl::PointXYZ> sourceCloud = lslgeneric::transformPointCloud(T,source);

    LazyGrid prototypeSource(resolutions.front());
    LazyGrid prototypeTarget(resolutions.front());

    NDTMap targetNDT( &prototypeTarget );
    targetNDT.loadPointCloud( target );
    targetNDT.computeNDTCells();

    NDTMap sourceNDT( &prototypeSource );
    sourceNDT.loadPointCloud( sourceCloud );
    sourceNDT.computeNDTCells();

    this->covariance(targetNDT,sourceNDT,TR,cov);

    return true;
}


}
