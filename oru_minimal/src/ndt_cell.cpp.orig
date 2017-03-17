#include <ndt_map/ndt_cell.h>
#include <cstring>
#include <cstdio>

namespace lslgeneric
{


bool NDTCell::parametersSet_ = false;
double NDTCell::EVAL_ROUGH_THR;
double NDTCell::EVEC_INCLINED_THR;
double NDTCell::EVAL_FACTOR;

void NDTCell::setParameters(double _EVAL_ROUGH_THR   ,
                                    double _EVEC_INCLINED_THR,
                                    double _EVAL_FACTOR
                                   )
{

    NDTCell::EVAL_ROUGH_THR    = _EVAL_ROUGH_THR   ;
    NDTCell::EVEC_INCLINED_THR = cos(_EVEC_INCLINED_THR);
    NDTCell::EVAL_FACTOR       = _EVAL_FACTOR      ;
    parametersSet_ = true;
}

/**
* produces a new Cell* of the same type
*/
NDTCell* NDTCell::clone() const
{
    NDTCell *ret = new NDTCell();
    ret->setDimensions(this->xsize_,this->ysize_,this->zsize_);
    ret->setCenter(this->center_);
    ret->setRGB(this->R,this->G,this->B);
    ret->setOccupancy(occ);
    ret->setEmptyval(emptyval);
    ret->setEventData(edata);
    ret->setN(this->N);
    ret->isEmpty = this->isEmpty;
    ret->hasGaussian_ = this->hasGaussian_;
    return ret;
}
/** produces a new Cell* of the same type and sets it to have the same
  * parameters as this cell.
  */
NDTCell* NDTCell::copy() const
{
    NDTCell *ret = new NDTCell();

    ret->setDimensions(this->xsize_,this->ysize_,this->zsize_);
    ret->setCenter(this->center_);

    for(unsigned int i=0; i<this->points_.size(); i++)
    {
        ret->points_.push_back(this->points_[i]);
    }

    ret->setMean(this->getMean());
    ret->setCov(this->getCov());
    ret->setRGB(this->R,this->G,this->B);
    ret->setOccupancy(this->occ);
    ret->setEmptyval(emptyval);
    ret->setEventData(edata);
    ret->setN(this->N);
    ret->isEmpty = this->isEmpty;
    ret->hasGaussian_ = this->hasGaussian_;
    ret->consistency_score = this->consistency_score;
    ret->cost = this->cost;
    return ret;
}
/**
* Updates the current Sample mean and covariance based on
* give new sample mean @m2 and covariance @cov2,
* which have been computed from @numpointsindistribution number of points
*/
void NDTCell::updateSampleVariance(const Eigen::Matrix3d &cov2,const Eigen::Vector3d &m2, unsigned int numpointsindistribution, 
	bool updateOccupancyFlag, float max_occu, unsigned int maxnumpoints)
{

    if(numpointsindistribution<=2){
	fprintf(stderr,"updateSampleVariance:: INVALID NUMBER OF POINTS\n");
	return;
    }
    if(this->hasGaussian_)
    {
	Eigen::Vector3d msum1 = mean_ * (double) N;
	Eigen::Vector3d msum2 = m2 * (double) numpointsindistribution;

	Eigen::Matrix3d csum1 = cov_ * (double) (N-1);
	Eigen::Matrix3d csum2 = cov2 * (double) (numpointsindistribution-1);

	if( fabsf(N) < 1e-5) {
	    fprintf(stderr,"Divider error (%u %u)!\n",N,numpointsindistribution);
	    hasGaussian_ = false;
	    return;
	}
	double divider = (double)numpointsindistribution+(double)N;
	if(fabs(divider) < 1e-5)
	{
	    fprintf(stderr,"Divider error (%u %u)!\n",N,numpointsindistribution);
	    return;
	}
	mean_ = (msum1 + msum2) / (divider);

	double w1 =  ((double) N / (double)(numpointsindistribution*(N+numpointsindistribution)));
	double w2 = (double) (numpointsindistribution)/(double) N;
	Eigen::Matrix3d	csum3 = csum1 + csum2 + w1 * (w2 * msum1 - msum2) * ( w2 * msum1 - msum2).transpose();
	N = N + numpointsindistribution;
	cov_ = 1.0/((double)N-1.0)  * csum3;
	if(updateOccupancyFlag){
	    double likoccval = 0.6;
	    double logoddlikoccu = numpointsindistribution * log((likoccval)/(1.0-likoccval));
	    updateOccupancy(logoddlikoccu, max_occu); 
	}
    }
    else
    {
	mean_ = m2;
	cov_ = cov2;
	N = numpointsindistribution;
	hasGaussian_ = true;
	if(updateOccupancyFlag){
	    double likoccval = 0.6;
	    double logoddlikoccu = numpointsindistribution * log((likoccval)/(1.0-likoccval));
	    updateOccupancy(logoddlikoccu,max_occu); 
	}
    }
    if(N>maxnumpoints) N=maxnumpoints;
    if(this->occ<0){
	this->hasGaussian_ = false;
	return;
    }
    rescaleCovariance();
}

///Just computes the normal distribution parameters from the points and leaves the points into the map 
void NDTCell::computeGaussianSimple(){
        Eigen::Vector3d meanSum_;
        Eigen::Matrix3d covSum_;
			
				if(points_.size()<6){
					points_.clear();
					return;
				}
				
        mean_<<0,0,0;
        for(unsigned int i=0; i< points_.size(); i++)
        {
            Eigen::Vector3d tmp;
            tmp<<points_[i].x,points_[i].y,points_[i].z;
            mean_ += tmp;
        }
        meanSum_ = mean_;
        mean_ /= (points_.size());
        Eigen::MatrixXd mp;
        mp.resize(points_.size(),3);
        for(unsigned int i=0; i< points_.size(); i++)
        {
            mp(i,0) = points_[i].x - mean_(0);
            mp(i,1) = points_[i].y - mean_(1);
            mp(i,2) = points_[i].z - mean_(2);
        }
        covSum_ = mp.transpose()*mp;
        cov_ = covSum_/(points_.size()-1);
        this->rescaleCovariance();
        N = points_.size();	
				R = 0; G = 0; B = 0;
				updateColorInformation();
}
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
/// Robust estimation using Student-T 
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
void NDTCell::studentT(){
	Eigen::Vector3d meanSum_, meantmp_;
	Eigen::Matrix3d covSum_, covTmp_;

	double nu = 5; // degrees of freedom of the t-distribution
	unsigned int maxIter = 10; // maximum number of iterations
	//mean_<<0,0,0;
	// weights for each sample point. Initialize each weight to 1
	std::vector<double> lambda;
	unsigned int pnts = points_.size();
	lambda.reserve(pnts);
	for(unsigned int i=0;i<pnts;i++) lambda[i] = 1.0; 
	
		
	for (unsigned int j=0; j<maxIter; j++)
	{
			// update mean
			double lambdaSum=0;
			meantmp_<<0,0,0;
			for(unsigned int i=0; i< pnts; i++)
			{
					Eigen::Vector3d tmp;
					tmp<<points_[i].x,points_[i].y,points_[i].z;
					meantmp_ += lambda[i]*tmp;
					lambdaSum += lambda[i];
			
			}
			
			meanSum_ = meantmp_;
			meantmp_ /= lambdaSum;
			
			// update scalematrix
			Eigen::MatrixXd mp;
			mp.resize(points_.size(),3);
			for(unsigned int i=0; i< pnts; i++)
			{
					double sqrtLambda = sqrt(lambda[i]);
					mp(i,0) = sqrtLambda*(points_[i].x - meantmp_(0));
					mp(i,1) = sqrtLambda*(points_[i].y - meantmp_(1));
					mp(i,2) = sqrtLambda*(points_[i].z - meantmp_(2));
			}
			covSum_ = mp.transpose()*mp;
			
			covTmp_ = covSum_/(points_.size());
			
			// compute inverse scalematrix
			Eigen::Matrix3d invCov;
			double det=0;
			bool exists=false;
			covTmp_.computeInverseAndDetWithCheck(invCov,det,exists);
			if(!exists){
				///The inverse does not exist -- exit
				return;
			}
			
			Eigen::Vector3d tempVec;
			// update the weights
			for (unsigned int i=0; i< points_.size(); i++){
					tempVec(0) = points_[i].x - meantmp_(0);
					tempVec(1) = points_[i].y - meantmp_(1);
					tempVec(2) = points_[i].z - meantmp_(2);
					double temp = nu;
					temp += squareSum(invCov, tempVec);
					lambda[i] = (nu+3)/(temp);
			}
	}
	double temp;
	temp = nu/(nu-2.0);
	covTmp_ = temp*covTmp_;
	
	if(!hasGaussian_){
		mean_ = meantmp_;
		cov_ = covTmp_;
		N = pnts;
		this->rescaleCovariance();
	}else{
		updateSampleVariance(covTmp_, meantmp_, pnts, false);
	}
	
	points_.clear();
}

/**
		Attempts to fit a gaussian in the cell.
		computes covariance and mean using observations
*/
void NDTCell::computeGaussian(int mode, unsigned int maxnumpoints, float occupancy_limit, Eigen::Vector3d origin, double sensor_noise)
{
    ///Occupancy update part
    ///This part infers the given "inconsistency" check done on update phase and
    ///updates the values accordingly
	
    ///Continous occupancy update
    double likoccval = 0.6;
    
    
    double logoddlikoccu = points_.size() * log((likoccval)/(1.0-likoccval));
//		if(mode != CELL_UPDATE_MODE_ERROR_REFINEMENT){
			if(logoddlikoccu > 0.4)  ///We have to trust that it is occupied if we have measurements from the cell!
			{
					updateOccupancy(logoddlikoccu, occupancy_limit);
			}
			else
			{
	#ifdef ICRA_2013_NDT_OM_SIMPLE_MODE
				double logoddlikempty = emptyval * log((0.49)/(1.0-0.49));
				updateOccupancy(logoddlikempty, occupancy_limit);	
	#else
				//updateOccupancy(logoddlikoccu+emptylik, occupancy_limit); << MUST BE ENABLED FOR OMG-NR
				updateOccupancy(logoddlikoccu, occupancy_limit);
				
	#endif
			}
//		}
    //occ++;
		isEmpty = -1;
		emptyval = 0;
		emptylik = 0;
		emptydist = 0;
		
		if(occ<=0){
      hasGaussian_ = false;
      return;
		}
		
/**
    if(points_.size()>=3)
    {
        //updateOccupancy(1.0, occupancy_limit);
        edata.updateSimple(EVENTMAP_OCCU); ///The cell is observed occupied
        isEmpty = -1;
        emptyval = 0;
        emptylik = 0;
        emptydist = 0;
    }
    else
    {
				points_.clear();
        if(emptyval>=6)
        {
            //updateOccupancy(-1.0, occupancy_limit);
            edata.updateSimple(EVENTMAP_FREE); ///The cell is observed empty
            isEmpty = 1;
            if(occ<=0)
            {
                hasGaussian_ = false;
            }
            //points_.clear();
            emptyval = 0;
            emptylik = 0;
            emptydist = 0;
            return;
        }
        else   ///No proper observation
        {
            isEmpty = 0;
            points_.clear();
            emptyval = 0;
            emptylik = 0;
            emptydist = 0;
            if(occ<0) hasGaussian_ = false;
            // 				//R = 0;
            return; //end
        }
    }
    ***/
    
		if((hasGaussian_==false && points_.size()< 3) || points_.size()==0){
			points_.clear();
			return; ///< not enough to compute the gaussian
		}
		
		if(mode==CELL_UPDATE_MODE_STUDENT_T){
			studentT();
			return;
		}


    if(!hasGaussian_)
    {
        Eigen::Vector3d meanSum_;
        Eigen::Matrix3d covSum_;

        mean_<<0,0,0;
        for(unsigned int i=0; i< points_.size(); i++)
        {
            Eigen::Vector3d tmp;
            tmp<<points_[i].x,points_[i].y,points_[i].z;
            mean_ += tmp;
        }
        meanSum_ = mean_;
        mean_ /= (points_.size());
        Eigen::MatrixXd mp;
        mp.resize(points_.size(),3);
        for(unsigned int i=0; i< points_.size(); i++)
        {
            mp(i,0) = points_[i].x - mean_(0);
            mp(i,1) = points_[i].y - mean_(1);
            mp(i,2) = points_[i].z - mean_(2);
        }
        covSum_ = mp.transpose()*mp;
        cov_ = covSum_/(points_.size()-1);
        this->rescaleCovariance();
        N = points_.size();
        

        /////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////
    }
    else if(mode == CELL_UPDATE_MODE_COVARIANCE_INTERSECTION)   ///Update using new information
    {
        /////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////
        Eigen::Vector3d m2;
        m2<<0,0,0;
        for(unsigned int i=0; i< points_.size(); i++)
        {
            Eigen::Vector3d tmp;
            tmp<<points_[i].x,points_[i].y,points_[i].z;
            m2 += tmp;
        }
        m2 /= (points_.size());


        Eigen::MatrixXd mp;
        mp.resize(points_.size(),3);
        for(unsigned int i=0; i< points_.size(); i++)
        {
            mp(i,0) = points_[i].x - m2(0);
            mp(i,1) = points_[i].y - m2(1);
            mp(i,2) = points_[i].z - m2(2);
        }
        Eigen::Matrix3d c2;
        c2 = mp.transpose()*mp/(points_.size()-1);
        double w1 = 0.98;
        Eigen::Matrix3d c3,icov2,icov3;
        bool exists = false;
        double det=0;
        exists = rescaleCovariance(c2, icov2);

        if(exists)
        {
            c3 = w1 * icov_ + (1.0-w1) * icov2;
            c3.computeInverseAndDetWithCheck(icov3,det,exists);
            if(exists)
            {
                cov_ = icov3;
                mean_ = icov3 * (w1*icov_*mean_ + (1.0-w1)*icov2*m2);

                this->rescaleCovariance();
                N += points_.size();
            }
            else
            {
                fprintf(stderr,"Covariance Intersection failed - Inverse does not exist (2)\n");
            }
        }
        else
        {
            points_.clear();
            fprintf(stderr,"Covariance Intersection failed - Inverse does not exist (1)\n");
        }
        /////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////
    }
    else if(mode == CELL_UPDATE_MODE_SAMPLE_VARIANCE)
    {
        /////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////
        Eigen::Vector3d meanSum_ = mean_ * N;
        Eigen::Matrix3d covSum_ = cov_ * (N-1);

        Eigen::Vector3d m2;
        m2<<0,0,0;
        for(unsigned int i=0; i< points_.size(); i++)
        {
            Eigen::Vector3d tmp;
            tmp<<points_[i].x,points_[i].y,points_[i].z;
            m2 += tmp;
        }
        Eigen::Vector3d T2 = m2;

        m2 /= (points_.size());

        Eigen::MatrixXd mp;
        mp.resize(points_.size(),3);
        for(unsigned int i=0; i< points_.size(); i++)
        {
            mp(i,0) = points_[i].x - m2(0);
            mp(i,1) = points_[i].y - m2(1);
            mp(i,2) = points_[i].z - m2(2);
        }

        Eigen::Matrix3d c2;
        c2 = mp.transpose()*mp;
        Eigen::Matrix3d c3;

        double w1 =  ((double) N / (double)(points_.size()*(N+points_.size())));
        double w2 = (double) (points_.size())/(double) N;

        c3 = covSum_ + c2 + w1 * (w2 * meanSum_ - (T2)) * ( w2 * meanSum_ - (T2)).transpose();

        meanSum_ = meanSum_ + T2;
        covSum_ = c3;
        N = N + points_.size();

        /**
         * If the timescaling is set.
         * This does "Sliding Average" for the sample mean and covariance
         */
        if(maxnumpoints > 0)
        {
            if(maxnumpoints < N)
            {
                meanSum_ = meanSum_ * ((double)maxnumpoints / (double)N);
                covSum_ = covSum_ * ((double)(maxnumpoints-1) / (double)(N-1));
                N = maxnumpoints;
            }
        }
        mean_ = meanSum_ / (double) N;
        cov_ = covSum_ / (double) (N-1);
        this->rescaleCovariance();

        /////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////
    }
    else if(mode == CELL_UPDATE_MODE_ERROR_REFINEMENT)
    {
        /////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////
        
        double e_min = 5.0e-4;
        double e_max = 5.0e-2;
        double o_max = 10.0;
        if(occ>o_max) occ=o_max;
        if(occ<-o_max) occ = -o_max;

        double epsilon = ((e_min - e_max) / (2.0*o_max)) * (occ+o_max)+e_max;

        for(unsigned int i=0; i< points_.size(); i++)
        {
            Eigen::Vector3d tmp;
            tmp<<points_[i].x,points_[i].y,points_[i].z;
            mean_ = mean_ + epsilon * (tmp - mean_);

            cov_ = cov_+epsilon * ((tmp-mean_) * (tmp-mean_).transpose() - cov_);
            //occ += 1.0;
            //if(occ>o_max) occ=o_max;
						//if(occ<-o_max) occ = -o_max;
        }
        hasGaussian_ = true;
        this->rescaleCovariance();

        /////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////
    }
    else if(mode == CELL_UPDATE_MODE_SAMPLE_VARIANCE_SURFACE_ESTIMATION)
    {
        /////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////
        Eigen::Vector3d meanSum_ = mean_ * N;
        Eigen::Matrix3d covSum_ = cov_ * (N-1);

        Eigen::Vector3d m2;
        m2<<0,0,0;
        for(unsigned int i=0; i< points_.size(); i++)
        {
            Eigen::Vector3d tmp;
            tmp<<points_[i].x,points_[i].y,points_[i].z;
            m2 += tmp;
        }
        Eigen::Vector3d T2 = m2;

        m2 /= (points_.size());

	pcl::PointXYZ orig;
        orig.x = origin[0];
        orig.y=origin[1];
        orig.z = origin[2];

        Eigen::Vector3d Xm;
        Eigen::Matrix3d c2=Eigen::Matrix3d::Zero();



        for(unsigned int i=0; i< points_.size(); i++)
        {
            Eigen::Vector3d X;
            X<<points_[i].x,points_[i].y,points_[i].z;

            computeMaximumLikelihoodAlongLine(orig, points_[i], Xm);
            double dist = (origin-X).norm();
            //double sigma_dist = 0.5 * ((dist*dist)/12.0); ///test for distance based sensor noise
            double sigma_dist = 0.5 * ((dist)/20.0); ///test for distance based sensor noise
            double snoise = sigma_dist + sensor_noise;
            double model_trust_given_meas = 0.49 + 0.5*exp(-0.5 * ((Xm-X).norm()*(Xm-X).norm())/(snoise*snoise));

            //double loglik_m = log(model_trust_given_meas / (1-model_trust_given_meas));
            double oval = occ/10.0;
            if(oval > 5.0) oval = 5.0;

            double likweight = oval;

            double weight = (1.0 - 1.0/(1.0+exp(likweight)))*model_trust_given_meas;
            if(N<100) weight = 0;
            c2 = c2 + (1.0-weight) * (X-m2)*(X-m2).transpose() + weight * (Xm-mean_)*(Xm-mean_).transpose();
        }


        Eigen::Matrix3d c3;

        double w1 =  ((double) N / (double)(points_.size()*(N+points_.size())));
        double w2 = (double) (points_.size())/(double) N;

        c3 = covSum_ + c2 + w1 * (w2 * meanSum_ - (T2)) * ( w2 * meanSum_ - (T2)).transpose();

        meanSum_ = meanSum_ + T2;
        covSum_ = c3;
        N = N + points_.size();

        /**
         * If the timescaling is set.
         * This does "Sliding Average" for the sample mean and covariance
         */
        if(maxnumpoints > 0)
        {
            if(maxnumpoints < N)
            {
                meanSum_ = meanSum_ * ((double)maxnumpoints / (double)N);
                covSum_ = covSum_ * ((double)(maxnumpoints-1) / (double)(N-1));
                N = maxnumpoints;
            }
        }
        mean_ = meanSum_ / (double) N;
        cov_ = covSum_ / (double) (N-1);
        this->rescaleCovariance();
        /////////////////////////////////////////////////////////////////////////////
        /////////////////////////////////////////////////////////////////////////////
    }///End of update
    /////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////
    updateColorInformation();
		
		points_.clear();
		/*
    if( !(hasGaussian_ && occ < 0))
    {
        points_.clear();
    }*/


}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/// Update the color information (template specialization)
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
///Default - Do nothing
void NDTCell::updateColorInformation()
{
}

//FIXME fix the template specialization here, is it still necessary?
#if 0
///RGB
inline void NDTCell::updateColorInformation()
{
    double r=0,g=0,b=0;
    for(unsigned int i=0; i< points_.size(); i++)
    {
        r+= ((double)points_[i].r)/255.0;
        g+= ((double)points_[i].g)/255.0;
        b+= ((double)points_[i].b)/255.0;
    }
    if(R == 0 && G == 0 && B == 0)
    {
        this->setRGB(r/(double)points_.size(), g/(double)points_.size(),b/(double)points_.size());
    }
    else
    {
        r /=  points_.size();
        g /=  points_.size();
        b /=  points_.size();


        R = 0.9*R + 0.1*r;
        G = 0.9*G + 0.1*g;
        B = 0.9*B + 0.1*b;
    }
}
///Intensity
inline void NDTCell<pcl::PointXYZI>::updateColorInformation()
{
    double intensity=0;
    for(unsigned int i=0; i< points_.size(); i++)
    {
        intensity+= ((double)points_[i].intensity)/255.0;
    }
    this->setRGB(intensity/(double)points_.size(), intensity/(double)points_.size(),intensity/(double)points_.size());
}
#endif


////////////////////////////////////////////////////////////////////////////////////
/**
* Rescales the covariance to protect against near sigularities
* and computes the inverse - This does not change class member values
*/
bool NDTCell::rescaleCovariance(Eigen::Matrix3d &cov, Eigen::Matrix3d &invCov)
{
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> Sol (cov);

    Eigen::Matrix3d evecs;
    Eigen::Vector3d evals;

    evecs = Sol.eigenvectors().real();
    evals = Sol.eigenvalues().real();

    if(evals(0) <= 0 || evals(1) <= 0 || evals(2) <= 0)
    {
        //fprintf(stderr,"Negative Eigenvalues!\n");
        hasGaussian_ = false;
        return false;
    }
    else
    {
        bool recalc = false;
        //guard against near singular matrices::
        int idMax;
        //double minEval = evals.minCoeff(&idMin);
        double maxEval = evals.maxCoeff(&idMax);
        if(maxEval > evals(0)*EVAL_FACTOR)
        {
            evals(0) = evals(idMax)/EVAL_FACTOR;
            recalc = true;
        }
        if(maxEval > evals(1)*EVAL_FACTOR)
        {
            evals(1) = evals(idMax)/EVAL_FACTOR;
            recalc = true;
        }
        if(maxEval > evals(2)*EVAL_FACTOR)
        {
            evals(2) = evals(idMax)/EVAL_FACTOR;
            recalc = true;
        }

        if(recalc)
        {
            Eigen::Matrix3d Lam;
            Lam = evals.asDiagonal();
            cov = evecs*Lam*(evecs.transpose());
        }

        //compute inverse covariance
        Eigen::Matrix3d Lam;
        Lam = evals.asDiagonal();
        invCov = evecs*(Lam.inverse())*(evecs.transpose());
    }
    return true;
}




///////////////////////////////////////////////////////////////////////////////////
void NDTCell::rescaleCovariance()
{
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> Sol (cov_);

    evecs_ = Sol.eigenvectors().real();
    evals_ = Sol.eigenvalues().real();

    //std::cout<<"evals = [ "<<evals_.transpose()<<"]';\n";
    if(evals_(0) <= 0 || evals_(1) <= 0|| evals_(2) <= 0)
    {
        //fprintf(stderr,"evals <=0\n");
        hasGaussian_ = false;
    }
    else
    {
        hasGaussian_ = true;
        
        bool recalc = false;
        //guard against near singular matrices::
        int idMax;
        //double minEval = evals.minCoeff(&idMin);
        double maxEval = evals_.maxCoeff(&idMax);
        if(maxEval > evals_(0)*EVAL_FACTOR)
        {
            evals_(0) = evals_(idMax)/EVAL_FACTOR;
            recalc = true;
        }
        if(maxEval > evals_(1)*EVAL_FACTOR)
        {
            evals_(1) = evals_(idMax)/EVAL_FACTOR;
            recalc = true;
        }
        if(maxEval > evals_(2)*EVAL_FACTOR)
        {
            evals_(2) = evals_(idMax)/EVAL_FACTOR;
            recalc = true;
        }

        if(recalc)
        {
            Eigen::Matrix3d Lam;
            Lam = evals_.asDiagonal();
            cov_ = evecs_*Lam*(evecs_.transpose());
        }
        classify();
        //compute inverse covariance
        Eigen::Matrix3d Lam;
        Lam = (evals_).asDiagonal();
        icov_ = evecs_*(Lam.inverse())*(evecs_.transpose());
    }
}


/** helper function for writeToJFF()
  */
void NDTCell::writeJFFMatrix(FILE * jffout, Eigen::Matrix3d &mat)
{

    double dtemp[6];

    dtemp[0] = mat.coeff(0,0);
    dtemp[1] = mat.coeff(1,0);
    dtemp[2] = mat.coeff(2,0);
    dtemp[3] = mat.coeff(1,1);
    dtemp[4] = mat.coeff(2,1);
    dtemp[5] = mat.coeff(2,2);

    fwrite(dtemp, sizeof(double), 6, jffout);

}

/** another helper function for writeToJFF()
  */
void NDTCell::writeJFFVector(FILE * jffout, Eigen::Vector3d &vec)
{

    double dtemp[3];

    for(int i=0; i<3; i++)
    {
        dtemp[i] = vec.coeff(i);
    }

    fwrite(dtemp, sizeof(double), 3, jffout);

}

/** yet another helper function for writeToJFF()
  */
void NDTCell::writeJFFEventData(FILE * jffout, TEventData &evdata)
{

    float    ftemp[4];
    uint8_t  ocval[1] = {evdata.occval};
    uint64_t evnts[1] = {evdata.events};

    ftemp[0] = evdata.a_exit_event;
    ftemp[1] = evdata.b_exit_event;
    ftemp[2] = evdata.a_entry_event;
    ftemp[3] = evdata.b_entry_event;

    fwrite(ocval, sizeof(uint8_t),  1, jffout);
    fwrite(ftemp, sizeof( float ),  4, jffout);
    fwrite(evnts, sizeof(uint64_t), 1, jffout);

}

/** output method to save the ndt cell as part of a jff v0.5 file
  */
int NDTCell::writeToJFF(FILE * jffout)
{

    pcl::PointXYZ * center = &(this->center_);
    fwrite(center, sizeof(pcl::PointXYZ), 1, jffout);
    double cell_size[3] = {this->xsize_, this->ysize_, this->zsize_};
    fwrite(cell_size, sizeof(double), 3, jffout);

    writeJFFMatrix(jffout, cov_);
    //writeJFFMatrix(jffout, covSum_);
    writeJFFVector(jffout, mean_);
    //writeJFFVector(jffout, meanSum_);

    // Temporary arrays to write all cell data
    double dtemp[2] = {d1_, d2_};
    //int    itemp[2] = {N, emptyval};
    int    itemp[3] = {N, emptyval,(int) hasGaussian_};
    float  ftemp[5] = {R, G, B, occ};

    fwrite(dtemp, sizeof(double), 2, jffout);
    fwrite(itemp, sizeof( int ),  3, jffout);
    fwrite(ftemp, sizeof(float),  4, jffout);

    writeJFFEventData(jffout, edata);

    return 0;

}

/** helper function for loadFromJFF()
  */
int NDTCell::loadJFFMatrix(FILE * jffin, Eigen::Matrix3d &mat)
{

    double dtemp[6];

    if(fread(&dtemp, sizeof(double), 6, jffin) <= 0)
        return -1;

    mat(0,0) = dtemp[0];
    mat(1,0) = dtemp[1];
    mat(2,0) = dtemp[2];
    mat(1,1) = dtemp[3];
    mat(2,1) = dtemp[4];
    mat(2,2) = dtemp[5];
    mat(0,1) = dtemp[1];
    mat(0,2) = dtemp[2];
    mat(1,2) = dtemp[4];

    return 0;

}

/** another helper function for loadFromJFF()
  */
int NDTCell::loadJFFVector(FILE * jffin, Eigen::Vector3d &vec)
{

    double dtemp[3];

    if(fread(&dtemp, sizeof(double), 3, jffin) <= 0)
        return -1;

    vec << dtemp[0], dtemp[1], dtemp[2];

    return 0;

}

/** yet another helper function for loadFromJFF()
  */
int NDTCell::loadJFFEventData(FILE * jffin, TEventData &evdata)
{

    float    ftemp[4];
    uint8_t  ocval;
    uint64_t evnts;

    if(fread(&ocval, sizeof(uint8_t),  1, jffin) <= 0)
        return -1;
    if(fread(&ftemp, sizeof( float ),  4, jffin) <= 0)
        return -1;
    if(fread(&evnts, sizeof(uint64_t), 1, jffin) <= 0)
        return -1;

    evdata.a_exit_event  = ftemp[0];
    evdata.b_exit_event  = ftemp[1];
    evdata.a_entry_event = ftemp[2];
    evdata.b_entry_event = ftemp[3];
    evdata.occval        = ocval;
    evdata.events        = evnts;

    return 0;
}

/** input method to load the ndt cell from a jff v0.5 file
  */
int NDTCell::loadFromJFF(FILE * jffin)
{

    pcl::PointXYZ center;
    if(fread(&center, sizeof(pcl::PointXYZ), 1, jffin) <= 0)
        return -1;
    this->setCenter(center);

    double dimensions[3];
    if(fread(&dimensions, sizeof(double), 3, jffin) <= 0)
        return -1;
    this->setDimensions(dimensions[0], dimensions[1], dimensions[2]);

    Eigen::Matrix3d temp_matrix;
    Eigen::Vector3d temp_vector;
    if(loadJFFMatrix(jffin, temp_matrix) < 0)
        return -1;
    this->setCov(temp_matrix);

    if(loadJFFVector(jffin, temp_vector) < 0)
        return -1;
    this->setMean(temp_vector);

    // Temporary arrays to load all cell data to
    double dtemp[2];// = {d1_, d2_};
    int    itemp[3];// = {N, emptyval, hasGaussian_};
    float  ftemp[5];// = {R, G, B, occ, cellConfidence};

    if(fread(&dtemp, sizeof(double), 2, jffin) <= 0)
        return -1;
    if(fread(&itemp, sizeof( int ),  3, jffin) <= 0)
        return -1;
    if(fread(&ftemp, sizeof(float),  4, jffin) <= 0)
        return -1;

    this->d1_ = dtemp[0];
    this->d2_ = dtemp[1];

    this->setN(itemp[0]);
    this->setEmptyval(itemp[1]);
    this->hasGaussian_ = (bool) itemp[2];

    this->setRGB(ftemp[0], ftemp[1], ftemp[2]);
    this->setOccupancy(ftemp[3]);

    TEventData edata;
    loadJFFEventData(jffin, edata);
    this->setEventData(edata);

    return 0;

}


/** classifies the cell according to the covariance matrix properties
    if smallest eigenval is bigger then roughness thershold, it's a rough cell
    evaluate inclination of the corresponding evector and classify as vertica, horizontal or inclined
  */
void NDTCell::classify()
{

    cl_ = UNKNOWN;

    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> tr;
    tr = tr.rotate(evecs_);

    int index=-1;
    double minEval = evals_.minCoeff(&index);
    if(index<0 || index > 2) return;

    if(minEval > EVAL_ROUGH_THR )
    {
        cl_ = ROUGH;
    }
    //if 1 eval << other 2 -> planar distr
    else
    {
        //default case -> sloping surface
        cl_ = INCLINED;

        //check the orientation of the vanishing axis
        Eigen::Vector3d e3;
        e3<<0,0,1;

        Eigen::Vector3d minorAxis = evecs_.col(index); //tr*e3;

        //dot product with vertical axis gives us the angle
        double d = minorAxis.dot(e3);
        double l = minorAxis.norm();
        double ac = d/l;
        if(fabsf(ac) < EVEC_INCLINED_THR)
        {
            //angle is nearly perpendicular => vertical surface
            cl_ = VERTICAL;
        }

        if(fabsf(ac) > 1-EVEC_INCLINED_THR)
        {
            //angle is nearly 0 => horizontal surface
            cl_ = HORIZONTAL;
        }
    }
}

double NDTCell::getLikelihood(const pcl::PointXYZ &pt) const
{
    //compute likelihood
    if(!hasGaussian_) return -1;
    Eigen::Vector3d vec (pt.x,pt.y,pt.z);
    vec = vec-mean_;
    double likelihood = vec.dot(icov_*vec);
    if(std::isnan(likelihood)) return -1;

    return exp(-likelihood/2);
    //return -d1*exp(-d2*likelihood/2);
}


/** setter for covariance
  */
void NDTCell::setCov(const Eigen::Matrix3d &_cov)
{
    cov_ = _cov;
    this->rescaleCovariance();
}
/**
* Computes the maximum likelihood that a point moving along a line
* defined by two points p1 and p2, gets measured agains the normaldistribution that
* is within this cell.
* This is used in raytracing to check if a measurement ray passes through a previously
* mapped object (thus provides evidence of inconsistency)
*
* @param p1 One point along the ray
* @param p2 second point along the ray (it must hold that p1 != p2);
*/
double NDTCell::computeMaximumLikelihoodAlongLine(const pcl::PointXYZ &p1, const pcl::PointXYZ &p2, Eigen::Vector3d &out)
{
    Eigen::Vector3d v1,v2;
    v1 << p1.x,p1.y,p1.z;
    v2 << p2.x,p2.y,p2.z;

    Eigen::Vector3d L = (v2-v1)/(v2-v1).norm();
    Eigen::Vector3d A = icov_ * L;
    Eigen::Vector3d B = (v2 - mean_);

    double sigma = A(0)*L(0)+A(1)*L(1)+A(2)*L(2);
    if(sigma == 0) return 1.0;

    double t = -(A(0)*B(0)+A(1)*B(1)+A(2)*B(2))/sigma;

    Eigen::Vector3d X = L*t+v2; ///X marks the spot

    pcl::PointXYZ p;
    p.x = X(0);
    p.y = X(1);
    p.z = X(2);
    out = X;
    return getLikelihood(p);

}




}; // end namespace
