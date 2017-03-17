#include <Eigen/Eigen>
#include <ndt_registration/ndt_matcher_d2d_feature.h>
#include <boost/bind.hpp>
#include <sys/time.h>

namespace lslgeneric
{

bool
NDTMatcherFeatureD2D::covariance( NDTMap& targetNDT,
        NDTMap& sourceNDT,
        Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T,
        Eigen::Matrix<double,6,6> &cov)
{
    assert(false);
};

//DEPRECATED???
double
NDTMatcherFeatureD2D::scoreNDT(std::vector<NDTCell*> &sourceNDT, NDTMap &targetNDT,
        Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor>& T)
{
    NUMBER_OF_ACTIVE_CELLS = 0;
    double score_here = 0;
    double det = 0;
    bool exists = false;
    NDTCell *cell;
    Eigen::Matrix3d covCombined, icov;
    Eigen::Vector3d meanFixed;
    Eigen::Vector3d meanMoving;
    Eigen::Matrix3d R = T.rotation();
    std::vector<std::pair<unsigned int, double> > scores;
    for(unsigned int j=0; j<_corr.size(); j++)
    {
        unsigned int i = _corr[j].second;
        if (_corr[j].second >= (int)sourceNDT.size())
        {
            std::cout << "second correspondance : " << _corr[j].second << ", " << sourceNDT.size() << std::endl;
        }
        if (sourceNDT[i] == NULL)
        {
            std::cout << "sourceNDT[i] == NULL!" << std::endl;
        }
        meanMoving = T*sourceNDT[i]->getMean();

        cell = targetNDT.getCellIdx(_corr[j].first);
        {

            if(cell == NULL)
            {
                std::cout << "cell== NULL!!!" << std::endl;
            }
            else
            {
                if(cell->hasGaussian_)
                {
                    meanFixed = cell->getMean();
                    covCombined = cell->getCov() + R.transpose()*sourceNDT[i]->getCov()*R;
                    covCombined.computeInverseAndDetWithCheck(icov,det,exists);
                    if(!exists) continue;
                    double l = (meanMoving-meanFixed).dot(icov*(meanMoving-meanFixed));
                    if(l*0 != 0) continue;
                    if(l > 120) continue;

                    double sh = -lfd1*(exp(-lfd2*l/2));

                    if(fabsf(sh) > 1e-10)
                    {
                        NUMBER_OF_ACTIVE_CELLS++;
                    }
                    scores.push_back(std::pair<unsigned int, double>(j, sh));
                    score_here += sh;
                    //score_here += l;
                }
            }
        }
    }

    if (_trimFactor == 1.)
    {
        return score_here;
    }
    else
    {
        // Determine the score value
        if (scores.empty()) // TODO, this happens(!), why??!??
            return score_here;

        score_here = 0.;
        unsigned int index = static_cast<unsigned int>(_trimFactor * (scores.size() - 1));
        //	std::nth_element (scores.begin(), scores.begin()+index, scores.end(), sort_scores()); //boost::bind(&std::pair<unsigned int, double>::second, _1) < boost::bind(&std::pair<unsigned int, double>::second, _2));
        std::nth_element (scores.begin(), scores.begin()+index, scores.end(), boost::bind(&std::pair<unsigned int, double>::second, _1) < boost::bind(&std::pair<unsigned int, double>::second, _2));
        std::fill(_goodCorr.begin(), _goodCorr.end(), false);
        //	std::cout << "_goodCorr.size() : " << _goodCorr.size() << " scores.size() : " << scores.size() << " index : " << index << std::endl;
        for (unsigned int i = 0; i < _goodCorr.size(); i++)
        {
            if (i <= index)
            {
                score_here += scores[i].second;
                _goodCorr[scores[i].first] = true;
            }
        }
        return score_here;
    }
}

double
NDTMatcherFeatureD2D::derivativesNDT( 
    const std::vector<NDTCell*> &sourceNDT,
    const NDTMap &targetNDT,
    Eigen::MatrixXd &score_gradient,
    Eigen::MatrixXd &Hessian,
    bool computeHessian
)
{

    struct timeval tv_start, tv_end;
    double score_here = 0;

    gettimeofday(&tv_start,NULL);
    NUMBER_OF_ACTIVE_CELLS = 0;
    score_gradient.setZero();
    Hessian.setZero();

    pcl::PointXYZ point;
    Eigen::Vector3d transformed;
    Eigen::Vector3d meanMoving, meanFixed;
    Eigen::Matrix3d CMoving, CFixed, CSum, Cinv, R;
    NDTCell *cell;
    bool exists = false;
    double det = 0;
    for (unsigned int j = 0; j < _corr.size(); j++)
    {
        if (!_goodCorr[j])
            continue;

        unsigned int i = _corr[j].second;
        if (i >= sourceNDT.size())
        {
            std::cout << "sourceNDT.size() : " << sourceNDT.size() << ", i: " << i << std::endl;
        }
        assert(i < sourceNDT.size());
        
	meanMoving = sourceNDT[i]->getMean();
        CMoving= sourceNDT[i]->getCov();
        
	this->computeDerivatives(meanMoving, CMoving, computeHessian);

        point.x = meanMoving(0);
        point.y = meanMoving(1);
        point.z = meanMoving(2);
       
        cell = targetNDT.getCellIdx(_corr[j].first);
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
	    double sh = -lfd1*(exp(-lfd2*l/2));
	    //std::cout<<"m1 = ["<<meanMoving.transpose()<<"]';\n m2 = ["<<cell->getMean().transpose()<<"]';\n";
	    //std::cout<<"C1 = ["<<CMoving<<"];\n C2 = ["<<CFixed<<"];\n";
	    //update score gradient
	    if(!this->update_gradient_hessian(score_gradient,Hessian,transformed, Cinv, sh, computeHessian))
	    {
		continue;
	    }
	    score_here += sh;
	    cell = NULL;
	}
    }
    gettimeofday(&tv_end,NULL);

    //double time_load = (tv_end.tv_sec-tv_start.tv_sec)*1000.+(tv_end.tv_usec-tv_start.tv_usec)/1000.;
    //std::cout<<"time derivatives took is: "<<time_load<<std::endl;
    return score_here;

}

}
