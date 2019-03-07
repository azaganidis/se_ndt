#include <ndt_map/pointcloud_utils.h>
#include <pcl/registration/transformation_estimation_svd.h>
#include <se_ndt/ndt_histogram.h>
#include <cmath>
#define n_threads 8
#define CROSS_ENTROPY

namespace perception_oru{

  void NDTHistogram::init         (int linear_classes,
                              int flat_classes,
                              int spherical_classes,
                              int n_classes,
	    		      double _D1,
			      double _D2) {
    N_LINE_BINS = linear_classes;
    N_FLAT_BINS = flat_classes;
    N_SPHERE_BINS = spherical_classes;
    N_CLASSES=n_classes;

    histogramBinsLine = Eigen::MatrixXi::Zero(N_LINE_BINS,N_CLASSES);
    histogramBinsFlat = Eigen::MatrixXi::Zero(N_FLAT_BINS,N_CLASSES);
    histogramBinsSphere = Eigen::MatrixXi::Zero(N_SPHERE_BINS,N_CLASSES);

    for(int i=0; i<3; i++)
      {
        dist_histogramBinsLine[i] = Eigen::MatrixXi::Zero(N_LINE_BINS,N_CLASSES);
        dist_histogramBinsFlat[i] = Eigen::MatrixXi::Zero(N_FLAT_BINS,N_CLASSES);
        dist_histogramBinsSphere[i] = Eigen::MatrixXi::Zero(N_SPHERE_BINS,N_CLASSES);
      }

    D1 = _D1;
    D2 = _D2;

    averageDirections = std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >(N_FLAT_BINS,Eigen::Vector3d(0,0,0));
    computeDirections();

    topThree.reserve(3);
    for(int r=0; r<3; r++)
      {
        topThree[r].setIdentity();
        topThreeS[r] = INT_MAX;
      }
    inited = true;
  }
  NDTHistogram::NDTHistogram (int linear_classes,
                              int flat_classes,
                              int spherical_classes,
                              int n_classes,
	    		      double _D1,
			      double _D2) {
    init(linear_classes, flat_classes, spherical_classes,n_classes, _D1, _D2);
  }
  NDTHistogram::NDTHistogram (NDTMap **map,
                              int linear_classes,
                              int flat_classes,
                              int spherical_classes,
                              int n_classes,
	    		      double _D1,
			      double _D2) {
    init(linear_classes, flat_classes, spherical_classes,n_classes, _D1, _D2);

    constructHistogram(map);
  }


  NDTHistogram::NDTHistogram (const NDTHistogram& other){

    N_LINE_BINS = other.N_LINE_BINS;
    N_FLAT_BINS = other.N_FLAT_BINS;
    N_SPHERE_BINS = other.N_SPHERE_BINS;
    N_CLASSES=other.N_CLASSES;

    histogramBinsLine =   other.histogramBinsLine;
    histogramBinsFlat =   other.histogramBinsFlat;
    histogramBinsSphere = other.histogramBinsSphere;

    for(int i=0; i<3; i++)
      {
        dist_histogramBinsLine[i] =   other.dist_histogramBinsLine[i];
        dist_histogramBinsFlat[i] =   other.dist_histogramBinsFlat[i];
        dist_histogramBinsSphere[i] = other.dist_histogramBinsSphere[i];
      }

    D1 = other.D1;
    D2 = other.D2;
    ENTROPY=other.ENTROPY;

    averageDirections = other.averageDirections;
    directions = other.directions;

    topThree.reserve(3);
    for(int r=0; r<3; r++){
      topThree[r].setIdentity();
      topThreeS[r] = INT_MAX;
    }
    inited = true;
  }


  void NDTHistogram::computeDirections(){

    double dlong = M_PI*(3-sqrt(5));  /* ~2.39996323 */
    double dz    = 2.0/N_FLAT_BINS;
    double longitude = 0;
    double z    = 1 - dz/2;

    for (int k = 0; k<N_FLAT_BINS; k++){
      double r    = sqrt(1-z*z);
      Eigen::Vector3d v;
      v<<cos(longitude)*r, sin(longitude)*r, z;
      directions.push_back(v);
      z    = z - dz;
      longitude = longitude + dlong;
    }
  }

  void NDTHistogram::constructHistogram(NDTMap **map,
                                        double linear_factor,
                                        double flat_factor
                                        ){

      omp_init_lock(&writelock);
    #pragma omp parallel num_threads(n_threads)
    {
        #pragma omp for
    for(int cl_ind=0;cl_ind<N_CLASSES; cl_ind++)
    {
        SpatialIndex *si = map[cl_ind]->getMyIndex();
        if(si==NULL) continue;

        // double LINEAR_FACTOR = 50;
        // double FLAT_FACTOR = 50;

        typename std::set<NDTCell*>::iterator it = si->begin();
        if(dynamic_cast<LazyGrid *>(si)->size()<1)
        {
            std::cerr<<"Empty map in histogram. (Semantic No "<<cl_ind<<")"<<std::endl;
            continue;
        }
        volatile int *sensor_pose = dynamic_cast<LazyGrid *>(si)->sensor_pose;
        double cellSizeX=dynamic_cast<LazyGrid *>(si)->cellSizeX;
        double cellSizeY=dynamic_cast<LazyGrid *>(si)->cellSizeY;
        double cellSizeZ=dynamic_cast<LazyGrid *>(si)->cellSizeZ;
        Eigen::Vector3d sensor_poseV;
        sensor_poseV(0)=sensor_pose[0]*cellSizeX;
        sensor_poseV(1)=sensor_pose[1]*cellSizeX;
        sensor_poseV(2)=sensor_pose[2]*cellSizeX;
        ///std::cerr<<cl_ind<<"\t NCELLS "<< dynamic_cast<LazyGrid *>(si)->size()<<std::endl;
        while(it!=si->end())
          {

            //NDTCell *ndcell = dynamic_cast<NDTCell* >(*it);
            if(*it == NULL){
                it++;
                continue;
              }
            if(!(*it)->hasGaussian_){
                it++;
                continue;
              }
            Eigen::Matrix3d evecs = (*it)->getEvecs();
            Eigen::Vector3d evals = (*it)->getEvals();

            int idMin,idMax,idMid;
            double minEval = evals.minCoeff(&idMin);
            double maxEval = evals.maxCoeff(&idMax);
            double midEval = -1;
            idMid = -1;
            for(int j=0; j<3; j++){
              if(j!=idMin && j!=idMax){
                midEval = evals(j);
                idMid = j;
              }
            }
            double dist = ((*it)->getMean()-sensor_poseV).norm();
            //three cases:
            //maxEval >> midEval -> linear
            if(maxEval > midEval*linear_factor){
              incrementLineBin(dist, cl_ind);
              it++;
              continue;
            }
            //maxEval ~ midEval >> minEval -> planar
            if(midEval > minEval*flat_factor){
                Eigen::Vector3d normal = evecs.col(idMin);
                Eigen::Vector3d mean = (*it)->getMean();
                if(normal.dot(mean) < 0){
                  //		std::cout<<"switching normal direction\n";
                  normal = -normal;
                }
                incrementFlatBin(normal,dist, cl_ind);
                it++;
                continue;
              }

            //maxEval ~ midEval ~ minEval -> spherical
            incrementSphereBin(dist, cl_ind);

            it++;
          }

      }
    }
    omp_destroy_lock(&writelock);
      for(int i=0; i<averageDirections.size(); i++)
          averageDirections[i].normalize();
      //std::cerr<<"H\t"<<histogramBinsFlat.sum()<<" "<<histogramBinsSphere.sum()<<" "<<histogramBinsLine.sum()<<std::endl;
  }

  void NDTHistogram::incrementLineBin(double d, int c){
    histogramBinsLine(0,c) ++;
    if(d<D1) dist_histogramBinsLine[0](0,c) ++;
    else if(d>D2) dist_histogramBinsLine[2](0,c) ++;
    else dist_histogramBinsLine[1](0,c) ++;
  }

  void NDTHistogram::incrementFlatBin(Eigen::Vector3d &normal, double d, int c){
    //std::cout<<"n "<<normal.transpose()<<std::endl;
    normal.normalize();
    //bins are in 3D. go through directions, find smallest difference
    double mindist = INT_MAX;
    int idmin = -1;
    for(unsigned int i=0; i<directions.size(); i++)
      {
        double dist = (directions[i]-normal).norm();
        if(mindist > dist)
          {
            mindist = dist;
            idmin = i;
          }
      }
    //std::cout<<idmin<<std::endl;
    if(idmin >=0 && idmin < N_FLAT_BINS)
      {
        histogramBinsFlat(idmin,c) ++;
        omp_set_lock(&writelock);
        averageDirections[idmin] += normal;
        omp_unset_lock(&writelock);
        if(d<D1) dist_histogramBinsFlat[0](idmin,c) ++;
        else if(d>D2) dist_histogramBinsFlat[2](idmin,c) ++;
        else dist_histogramBinsFlat[1](idmin,c) ++;
      }
  }

  void NDTHistogram::incrementSphereBin(double d, int c){
    histogramBinsSphere(0,c) ++;
    if(d<D1){
      int id = floor(((double)d*N_SPHERE_BINS)/(double)D1);
      dist_histogramBinsSphere[0](id,c) ++;
    }
    else if(d>D2){
      dist_histogramBinsSphere[2](0,c) ++;
    }
    else{
      int id = floor(((double)(d-D1)*N_SPHERE_BINS)/(double)D2);
      dist_histogramBinsSphere[1](id,c) ++;
    }
  }


  pcl::PointCloud<pcl::PointXYZ> NDTHistogram::getDominantDirections(int nDirections){
    
    pcl::PointCloud<pcl::PointXYZ> ret;
    std::vector<bool> dominated (directions.size(),false);
    double NORM_MIN = 0.2;
    int MIN_SUPPORT = 3;
    
    for(int i=0; i<nDirections; i++){
      //get the next direction
      pcl::PointXYZ current;
      //find max in histogram, that is not dominated
      bool found = false;
      int maxBin, idMax;
      Eigen::MatrixXi histBinsFlatSum = histogramBinsFlat.rowwise().sum();
      while (!found){
        maxBin = -1;
        idMax = -1;
        for(int j=0; j<N_FLAT_BINS; j++){
          if(histBinsFlatSum(j) > maxBin && !dominated[j]){
            maxBin = histBinsFlatSum(j);
            idMax = j;
          }
        }

        found = !dominated[idMax];
        //check if any of the already found directions are "duals"
        /* for(int j=0; j<ret.points.size() && found; j++) {
           Eigen::Vector3d v(ret.points[j].x,ret.points[j].y,ret.points[j].z);
           v = v + directions[idMax];
           found = found && (v.norm() > NORM_MIN);
           }*/
        //suppress this max and neighbours --- independent of the value of found, this is necessarry
        dominated[idMax] = true;
        if(idMax-1 >=0) dominated[idMax-1] = true;
        if(idMax+1 <dominated.size()) dominated[idMax+1] = true;
        if(maxBin < MIN_SUPPORT) break;
      }

      if(maxBin < MIN_SUPPORT) break;
      //current.x = directions[idMax](0);
      //current.y = directions[idMax](1);
      //current.z = directions[idMax](2);
      current.x = averageDirections[idMax](0);
      current.y = averageDirections[idMax](1);
      current.z = averageDirections[idMax](2);
      //current.intensity = maxBin;
      //std::cout<<directions[idMax].transpose()<<" e "<<maxBin<<std::endl;
      //std::cout<<averageDirections[idMax].transpose()<<" e "<<maxBin<<std::endl;
      ret.points.push_back(current);
    }
    return ret;
  }

  void NDTHistogram::bestFitToHistogram(NDTHistogram &target, Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> &T, bool bound_transform)
  {

    //find the top N dominant directions
    int N = 3;
    //store vectors as pcl::points
    pcl::PointCloud<pcl::PointXYZ> dominantBinsMine, dominantBinsTarget;

    //    std::cout<<"d1 : \n";
    dominantBinsMine = this->getDominantDirections(N);
    //    std::cout<<"d2 : \n";
    dominantBinsTarget = target.getDominantDirections(N);

    /*    double N_THIS=0, N_OTHER=0;
          for(int i=0; i<histogramBinsFlat.size(); i++) {
          N_THIS += histogramBinsFlat[i];
          N_OTHER += target.histogramBinsFlat[i];
          }
    */    //estimate least-squares fit, assuming correspondence
    /*    pcl::registration::TransformationEstimationSVD<PointT,PointT> trEst;
          Eigen::Matrix4f TR;
          trEst.estimateRigidTransformation(dominantBinsMine, dominantBinsTarget, TR);
          T = TR.cast<double>();
    */
    //check for best fitting combination

    for(int r=0; r<3; r++)
      {
        topThree[r].setIdentity();
        topThreeS[r] = INT_MAX;
      }
    pcl::PointCloud<pcl::PointXYZ> mineNew, mineNew2;
    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> localT, localT2;

    //    double best = INT_MAX;
    //choose 3 out of N
    for(int a1 = 0; a1<dominantBinsMine.points.size(); a1++){
      //int a1 = 0; {
      for(int b1 = a1+1; b1<dominantBinsMine.points.size(); b1++){
        //if(a1 == b1) continue;
        //for(int c1 = 0; c1<dominantBinsMine.points.size(); c1++) {
        //if(b1 == c1 || a1 == c1) continue;
        //int a2 = 0; {
        for(int a2 = 0; a2<dominantBinsTarget.points.size(); a2++)
          {
            for(int b2 = 0; b2<dominantBinsTarget.points.size(); b2++)
              {
                if(a2 == b2) continue;
                //for(int c2 = 0; c2<dominantBinsTarget.points.size(); c2++) {
                //    if(b2 == c2 || a2 == c2) continue;
                pcl::PointCloud<pcl::PointXYZ> useMine, useTarget;
                useMine.points.push_back(dominantBinsMine[a1]);
                useMine.points.push_back(dominantBinsMine[b1]);
                //useMine.points.push_back(dominantBinsMine[c1]);
                useTarget.points.push_back(dominantBinsTarget[a2]);
                useTarget.points.push_back(dominantBinsTarget[b2]);
                //useTarget.points.push_back(dominantBinsTarget[c2]);

                closedFormSolution(useMine, useTarget, localT2);
                double good2 = this->getSimilarity(target,localT2);
                //	    std::cout<<good2<<std::endl;
                //compute score of fitness "good2"
                /*			    double good2 = 0, scale2 = 0;
                                for(int q = 0; q<averageDirections.size(); q++) {

                                Eigen::Vector3d tr = localT2*averageDirections[q];
                                if( this->histogramBinsFlat[q] == 0) {
                                tr = directions[q]; //fall back to bin direction
                                }
                                //find B = the bin in which tr falls
                                tr.normalize();
                                //std::cout<<"TR"<<tr.transpose()<<std::endl;
                                double mindist = INT_MAX;
                                int idmin = -1;
                                for(unsigned int i=0; i<directions.size(); i++) {
                                double dist = (directions[i]-tr).norm();
                                if(mindist > dist) {
                                mindist = dist;
                                idmin = i;
                                }
                                }
                                //std::cout<<idmin<<std::endl;
                                if(!(idmin >=0 && idmin < histogramBinsFlat.size())) {
                                continue;
                                }


                                //find the averageDirection other->average[B]
                                Eigen::Vector3d other_tr = target.averageDirections[idmin];
                                if (target.histogramBinsFlat[idmin] == 0) {
                                other_tr = target.directions[idmin]; //fall back to bin direction
                                }

                                //compute norm of difference, scale by cardinality
                                //double factor = fabsf((double)(this->histogramBinsFlat[q] - target.histogramBinsFlat[idmin]+1)/
                                //			(double)(this->histogramBinsFlat[q]+target.histogramBinsFlat[idmin]+1));
                                //double factor = this->histogramBinsFlat[q]+target.histogramBinsFlat[idmin];
                                //scale2 += factor;
                                //good2 += factor*(tr-other_tr).norm();
                                good2 += pow((double)this->histogramBinsFlat[q]/N_THIS - (double)target.histogramBinsFlat[idmin]/N_OTHER,2);

                                }
                                good2 = sqrt(good2);
                */

                /*
                //don't use them for finding but just for verifying
                useMine.points.push_back(dominantBinsMine[c1]);
                useTarget.points.push_back(dominantBinsTarget[c2]);
                //calculate goodness of fit
                //			    mineNew = transformPointCloud(localT,useMine);

                //now also considering last orientation....
                //			    closedFormSolution(useMine, useTarget, localT2);
                mineNew2 = transformPointCloud(localT2,useMine);

                //			    double good = 0, scale = 0;
                //			    for(int i=0; i<mineNew.points.size(); i++) {
                //				double factor = (mineNew.points[i].intensity+useTarget.points[i].intensity);
                //				good += factor*sqrt( pow(mineNew.points[i].x-useTarget.points[i].x,2) +
                //					pow(mineNew.points[i].y-useTarget.points[i].y,2) +
                //					pow(mineNew.points[i].z-useTarget.points[i].z,2));
                //				scale += factor;
                //			    }
                //			    good = good/scale;


                double good2 = 0, scale2 = 0;
                for(int i=0; i<mineNew2.points.size(); i++) {
                double factor = (mineNew2.points[i].intensity+useTarget.points[i].intensity);
                good2 += factor*sqrt( pow(mineNew2.points[i].x-useTarget.points[i].x,2) +
                pow(mineNew2.points[i].y-useTarget.points[i].y,2) +
                pow(mineNew2.points[i].z-useTarget.points[i].z,2));
                scale2 += factor;
                }
                good2 = good2/scale2;
                */

                //			    std::cout<<"combo "<<a1<<" "<<b1<<" "<<" -- "
                //				<<a2<<" "<<b2<<" "<<" fit = "<<good2<<std::endl;
                /*if(good < best) {
                  std::cout<<"local minimum at combo "<<a1<<" "<<b1<<" "<<c1<<" -- "
                  <<a2<<" "<<b2<<" "<<c2<<" fit = "<<good<<std::endl;
                  best = good;
                  T = localT;
                  }*/
                Eigen::Quaternion<double> id, errorQ;
                id.setIdentity();
                errorQ = localT2.rotation();
                double angle = fabsf(acos(id.dot(errorQ))/2);
                if(angle > M_PI/8 && bound_transform)
                  {
                    //transform is too big!
                    continue;
                  }
                if(good2 < topThreeS[2])
                  {
                    //				std::cout<<"local minimum at combo "<<a1<<" "<<b1<<" "<<" -- "
                    //				    <<a2<<" "<<b2<<" "<<" fit = "<<good2<<std::endl;
                    if(good2 < topThreeS[1])
                      {
                        if(good2 < topThreeS[0])
                          {
                            topThree[2] = topThree[1];
                            topThree[1] = topThree[0];
                            topThree[0] = localT2;
                            topThreeS[2] = topThreeS[1];
                            topThreeS[1] = topThreeS[0];
                            topThreeS[0] = good2;
                          }
                        else
                          {
                            topThree[2] = topThree[1];
                            topThree[1] = localT2;
                            topThreeS[2] = topThreeS[1];
                            topThreeS[1] = good2;
                          }
                      }
                    else
                      {
                        topThree[2] = localT2;
                        topThreeS[2] = good2;
                      }
                  }
                //			}
              }
          }
        //    }
      }
    }

    T = topThree[0];

    if(dominantBinsMine.points.size() < 2 || dominantBinsTarget.points.size() < 2)
      {
        T.setIdentity();
      }

  }

  void NDTHistogram::closedFormSolution(pcl::PointCloud<pcl::PointXYZ> &src, pcl::PointCloud<pcl::PointXYZ> &tgt, Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> &T){
    T.setIdentity ();
    Eigen::Matrix3d H; // = (cloud_src_demean * cloud_tgt_demean.transpose ()).topLeftCorner<3, 3>();
    Eigen::MatrixXd P1, P2;         //temporary points for sets 1 and 2
    unsigned int itr=0;             //loop counters
    size_t size = src.points.size();

    // Assemble the correlation matrix H = source * target'
    P1 = Eigen::MatrixXd(size,3);
    P2 = Eigen::MatrixXd(size,3);
    //compute values needed for the N matrix and for scale
    for(itr=0; itr<size; itr++)
      {
        P1(itr,0) = src.points[itr].x;
        P1(itr,1) = src.points[itr].y;
        P1(itr,2) = src.points[itr].z;
        P2(itr,0) = tgt.points[itr].x;
        P2(itr,1) = tgt.points[itr].y;
        P2(itr,2) = tgt.points[itr].z;
      }
    H = P1.transpose()*P2;

    // Compute the Singular Value Decomposition
    Eigen::JacobiSVD<Eigen::Matrix3d> svd (H, Eigen::ComputeFullU | Eigen::ComputeFullV);
    Eigen::Matrix3d u = svd.matrixU ();
    Eigen::Matrix3d v = svd.matrixV ();

    // Compute R = V * U'
    if (u.determinant () * v.determinant () < 0)
      {
        for (int x = 0; x < 3; ++x)
          v (x, 2) *= -1;
      }

    Eigen::Matrix3d R = v * u.transpose ();
    // Return the correct transformation
    T = R;

  }

  double NDTHistogram::getSimilarity(NDTHistogram &other)
  {

    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> T;

    this->bestFitToHistogram(other,T,false);

    return this->getSimilarity(other,T);
  }

double NDTHistogram::calculateEntropy()
{
    double sum=0;
    double entropy=0;
    for(int r=0;r<2;r++)
    {
        sum+= this->dist_histogramBinsFlat[r].sum();
        sum+= this->dist_histogramBinsSphere[r].sum();
        sum+= this->dist_histogramBinsLine[r].row(0).sum();
    }
    for(int r=0;r<2;r++)
    {
        Eigen::MatrixXd F_=(dist_histogramBinsFlat[r].cast<double>()/sum);
        entropy+=(F_.array()*F_.array().log().unaryExpr([](double v) { return std::isfinite(v)? v : 0.0; })).sum();

        Eigen::MatrixXd S_=(dist_histogramBinsSphere[r].cast<double>()/sum);
        entropy+=(S_.array()*S_.array().log().unaryExpr([](double v) { return std::isfinite(v)? v : 0.0; })).sum();

        Eigen::MatrixXd L_=(dist_histogramBinsLine[r].row(0).cast<double>()/sum);
        entropy+=(L_.array()*L_.array().log().unaryExpr([](double v) { return std::isfinite(v)? v : 0.0; })).sum();
    }
//    return -entropy;
    if(std::isnan(entropy))//FIXME sum=0
        std::cerr<<sum<<" "<<std::endl;
    assert(!std::isnan(entropy));
    double num_el=N_LINE_BINS*N_SPHERE_BINS*N_FLAT_BINS*N_CLASSES*2;
    double max_entropy=-log(1/num_el);
    double normalized_entropy = -entropy/max_entropy;
    ENTROPY=normalized_entropy;
    return normalized_entropy;
}
/* //FOR CROSS ENTROPY
double entropySafe(double p, double q){
    if(q<=0)
        return -p*1000;
    return p*log(q);
}
*/
//FOR KL DIVERGENCE
double entropySafe(double p, double q){
    if(p==0)
        return 0;
    return p*log(q/p);
}

  double NDTHistogram::getSimilarity(NDTHistogram &other, Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> &T){
    double score[3];
    double scoreT=0;
    double score_final[3]={0,0,0};
    Eigen::MatrixXd N_THIS=Eigen::MatrixXd::Zero(3, N_CLASSES);
    Eigen::MatrixXd N_OTHER=Eigen::MatrixXd::Zero(3, N_CLASSES);
    Eigen::MatrixXd scoreM=Eigen::MatrixXd::Zero(3, N_CLASSES);
    for(unsigned int r = 0; r<3; r++){
        N_THIS.row(r)= this->dist_histogramBinsFlat[r].colwise().sum().cast<double> ();
        N_OTHER.row(r)= other.dist_histogramBinsFlat[r].colwise().sum().cast<double> ();

        N_THIS.row(r)+= this->dist_histogramBinsSphere[r].colwise().sum().cast<double> ();
        N_OTHER.row(r)+= other.dist_histogramBinsSphere[r].colwise().sum().cast<double> ();

        N_THIS.row(r)+= this->dist_histogramBinsLine[r].row(0).cast<double> ();
        N_OTHER.row(r)+= other.dist_histogramBinsLine[r].row(0).cast<double> ();
    }
    for(int i_class=0;i_class<N_CLASSES;i_class++)
    {
        for(unsigned int r = 0; r<3; r++)
        {
            score[r]=0;
            //if((N_THIS(r, i_class)==0) != (N_OTHER(r,i_class)==0))
            //    score[r]=1;
        }
        for(int q = 0; q<averageDirections.size(); q++)
        {

            Eigen::Vector3d tr = T*averageDirections[q];
            if( this->histogramBinsFlat(q, i_class) == 0)
            {
                tr = directions[q]; //fall back to bin direction
            }

            //find B = the bin in which tr falls
            tr.normalize();
            double mindist = INT_MAX;
            int idmin = -1;
            for(unsigned int i=0; i<directions.size(); i++)
            {
                double dist = (directions[i]-tr).norm();
                if(mindist > dist)
                  {
                    mindist = dist;
                    idmin = i;
                  }
            }
            //std::cout<<idmin<<std::endl;
            if(!(idmin >=0 && idmin < N_FLAT_BINS))
            {
                continue;
            }

            for(unsigned int r = 0; r<3; r++)
            {
                double d1=0,d2=0;
                if(N_THIS(r, i_class)!=0)
                    d1 = ((double)this->dist_histogramBinsFlat[r](q, i_class))/N_THIS(r,i_class);
                if(N_OTHER(r, i_class)!=0)
                    d2 = ((double)other.dist_histogramBinsFlat[r](idmin, i_class))/N_OTHER(r,i_class);
#ifdef CROSS_ENTROPY
                d1 = ((double)this->dist_histogramBinsFlat[r](q, i_class));
                d2 = ((double)other.dist_histogramBinsFlat[r](idmin, i_class));
                if(d2==0)//d1 can be zero, with correct result.
                    d2=0.01;//Set a small number.
                d1=d1/N_THIS.sum();
                d2=d2/N_OTHER.sum();
                scoreT-=entropySafe(d1,d2);
#else
                score[r] += pow( d1-d2 ,2);
#endif
            }

        }
        for(unsigned int r = 0; r<3; r++)
          {
              //scoreM.row(r)=(this->dist_histogramBinsSphere[r].cast<double>().array()/N_THIS.row(r).array() - other.dist_histogramBinsSphere[r].cast<double>().array()/N_OTHER.row(r).array()).pow(2).colwise().sum();
            for(int i=0; i<N_SPHERE_BINS; i++)
              {
                double d1=0,d2=0;
                if(N_THIS(r, i_class)!=0)
                    d1=((double)this->dist_histogramBinsSphere[r](i, i_class))/N_THIS(r,i_class);
                if(N_OTHER(r, i_class)!=0)
                    d2=((double)other.dist_histogramBinsSphere[r](i, i_class))/N_OTHER(r,i_class);
#ifdef CROSS_ENTROPY
                d1 = ((double)this->dist_histogramBinsSphere[r](i, i_class));
                d2 = ((double)other.dist_histogramBinsSphere[r](i, i_class));
                if(d2==0)//d1 can be zero, with correct result.
                    d2=0.01;//Set a small number.
                d1=d1/N_THIS.sum();
                d2=d2/N_OTHER.sum();
                scoreT-=entropySafe(d1,d2);
#else
                score[r] += pow( d1-d2 ,2);
#endif
              }
            //score[r]+=scoreM(r,i_class);

            double d1=0,d2=0;
            if(N_THIS(r, i_class)!=0)
                d1=((double)this->dist_histogramBinsLine[r](0, i_class))/N_THIS(r,i_class);
            if(N_OTHER(r, i_class)!=0)
                d2=((double)other.dist_histogramBinsLine[r](0, i_class))/N_OTHER(r,i_class);
#ifdef CROSS_ENTROPY
            d1 = ((double)this->dist_histogramBinsLine[r](0, i_class));
            d2 = ((double)other.dist_histogramBinsLine[r](0, i_class));
            if(d2==0)//d1 can be zero, with correct result.
                d2=0.01;//Set a small number.
            d1=d1/N_THIS.sum();
            d2=d2/N_OTHER.sum();
            scoreT-=entropySafe(d1,d2);
#else
            score[r] += pow( d1-d2 ,2);
#endif
            double maxN, minN;
            maxN = (N_THIS(r,i_class) > N_OTHER(r,i_class)) ? N_THIS(r,i_class) : N_OTHER(r,i_class);
            minN = (N_THIS(r,i_class) < N_OTHER(r,i_class)) ? N_THIS(r,i_class) : N_OTHER(r,i_class);
            minN = (minN < 1) ? 1 : minN;
            //double scale_factor = -log(minN/maxN);//ANESTIS CONTRIB

#ifndef CROSS_ENTROPY
            score_final[r] += maxN*sqrt(score[r])/minN/N_CLASSES;
#endif
          }
    }

//    if(score_final[0]+score_final[1]==0)
//        std::cout<<"EDW"<<std::endl;

#ifdef CROSS_ENTROPY
    //return score_final[0]+score_final[1];
    return scoreT;

#else
    return (score_final[0]+score_final[1]+score_final[2])/3;
#endif
  }

}
