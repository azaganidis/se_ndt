#include <ndt_map/pointcloud_utils.h>
#include <pcl/registration/transformation_estimation_svd.h>
#include <ndt_map/ndt_histogram.h>

namespace lslgeneric{

  NDTHistogram::NDTHistogram(){
    N_LINE_BINS = 1;
    N_FLAT_BINS = 40;
    N_SPHERE_BINS = 10;

    histogramBinsLine = std::vector<int>(N_LINE_BINS,0);
    histogramBinsFlat = std::vector<int>(N_FLAT_BINS,0);
    histogramBinsSphere = std::vector<int>(N_SPHERE_BINS,0);

    for(int i=0; i<3; i++)
      {
        dist_histogramBinsLine[i] = std::vector<int>(N_LINE_BINS,0);
        dist_histogramBinsFlat[i] = std::vector<int>(N_FLAT_BINS,0);
        dist_histogramBinsSphere[i] = std::vector<int>(N_SPHERE_BINS,0);
      }

    D1 = 5;
    D2 = 10;

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

  NDTHistogram::NDTHistogram (const NDTHistogram& other){

    histogramBinsLine =   other.histogramBinsLine;
    histogramBinsFlat =   other.histogramBinsFlat;
    histogramBinsSphere = other.histogramBinsSphere;

    for(int i=0; i<3; i++)
      {
        dist_histogramBinsLine[i] =   other.dist_histogramBinsLine[i];
        dist_histogramBinsFlat[i] =   other.dist_histogramBinsFlat[i];
        dist_histogramBinsSphere[i] = other.dist_histogramBinsSphere[i];
      }

    D1 = 5;
    D2 = 10;

    averageDirections = other.averageDirections;
    directions = other.directions;

    topThree.reserve(3);
    for(int r=0; r<3; r++){
      topThree[r].setIdentity();
      topThreeS[r] = INT_MAX;
    }
    inited = true;
  }

  NDTHistogram::NDTHistogram (NDTMap &map){

    N_LINE_BINS = 1;
    N_FLAT_BINS = 40;
    N_SPHERE_BINS = 10;

    histogramBinsLine = std::vector<int>(N_LINE_BINS,0);
    histogramBinsFlat = std::vector<int>(N_FLAT_BINS,0);
    histogramBinsSphere = std::vector<int>(N_SPHERE_BINS,0);

    for(int i=0; i<3; i++){
      dist_histogramBinsLine[i] = std::vector<int>(N_LINE_BINS,0);
      dist_histogramBinsFlat[i] = std::vector<int>(N_FLAT_BINS,0);
      dist_histogramBinsSphere[i] = std::vector<int>(N_SPHERE_BINS,0);
    }

    D1 = 5;
    D2 = 10;

    averageDirections = std::vector<Eigen::Vector3d,Eigen::aligned_allocator<Eigen::Vector3d> >(N_FLAT_BINS,Eigen::Vector3d(0,0,0));
    //populate directions
    computeDirections();
    constructHistogram(map);

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


  void NDTHistogram::constructHistogram(NDTMap &map){

    SpatialIndex *si = map.getMyIndex();
    if(si==NULL) return;

    double LINEAR_FACTOR = 50;
    double FLAT_FACTOR = 50;

    typename std::vector<NDTCell*>::iterator it = si->begin();
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
        double dist = (*it)->getMean().norm();
        //three cases:
        //maxEval >> midEval -> linear
        if(maxEval > midEval*LINEAR_FACTOR){
          incrementLineBin(dist);
          it++;
          continue;
        }
        //maxEval ~ midEval >> minEval -> planar
        if(midEval > minEval*FLAT_FACTOR){
            Eigen::Vector3d normal = evecs.col(idMin);
            Eigen::Vector3d mean = (*it)->getMean();
            if(normal.dot(mean) < 0){
              //		std::cout<<"switching normal direction\n";
              normal = -normal;
            }
            incrementFlatBin(normal,dist);
            it++;
            continue;
          }

        //maxEval ~ midEval ~ minEval -> spherical
        incrementSphereBin(dist);

        it++;
      }

    for(int i=0; i<averageDirections.size(); i++)
      {
        averageDirections[i].normalize();
      }

  }

  void NDTHistogram::incrementLineBin(double d){
    histogramBinsLine[0] ++;
    if(d<D1) dist_histogramBinsLine[0][0] ++;
    else if(d>D2) dist_histogramBinsLine[2][0] ++;
    else dist_histogramBinsLine[1][0] ++;
  }

  void NDTHistogram::incrementFlatBin(Eigen::Vector3d &normal, double d){
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
    if(idmin >=0 && idmin < histogramBinsFlat.size())
      {
        histogramBinsFlat[idmin] ++;
        averageDirections[idmin] += normal;
        if(d<D1) dist_histogramBinsFlat[0][idmin] ++;
        else if(d>D2) dist_histogramBinsFlat[2][idmin] ++;
        else dist_histogramBinsFlat[1][idmin] ++;
      }
  }

  void NDTHistogram::incrementSphereBin(double d){
    histogramBinsSphere[0] ++;
    if(d<D1){
      int id = floor(((double)d*N_SPHERE_BINS)/(double)D1);
      dist_histogramBinsSphere[0][id] ++;
    }
    else if(d>D2){
      dist_histogramBinsSphere[2][0] ++;
    }
    else{
      int id = floor(((double)(d-D1)*N_SPHERE_BINS)/(double)D2);
      dist_histogramBinsSphere[1][id] ++;
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
      while (!found){
        maxBin = -1;
        idMax = -1;
        for(int j=0; j<histogramBinsFlat.size(); j++){
          if(histogramBinsFlat[j] > maxBin && !dominated[j]){
            maxBin = histogramBinsFlat[j];
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

  void NDTHistogram::printHistogram(bool bMatlab){
    if(bMatlab){
      //prints in a format suitable for matlab plotting
      std::cout<<"L=[ ";
      for(unsigned int i=0; i<histogramBinsLine.size(); i++){
        std::cout<<histogramBinsLine[i]<<" ";
      }
      std::cout<<"];\n";
      std::cout<<"F=[";
      for(unsigned int i=0; i<histogramBinsFlat.size(); i++)
        {
          std::cout<<histogramBinsFlat[i]<<" ";
        }
      std::cout<<"];\n";
      for(unsigned int q=0; q<3; q++)
        {
          std::cout<<"F"<<q<<" = [";
          for(unsigned int i=0; i<dist_histogramBinsFlat[q].size(); i++)
            {
              std::cout<<dist_histogramBinsFlat[q][i]<<" ";
            }
          std::cout<<"];\n";
        }

      std::cout<<"];\nS=[";
      for(unsigned int i=0; i<histogramBinsSphere.size(); i++)
        {
          std::cout<<histogramBinsSphere[i]<<" ";
        }
      std::cout<<"];\n";
      for(unsigned int q=0; q<3; q++)
        {
          std::cout<<"S"<<q<<" = [";
          for(unsigned int i=0; i<dist_histogramBinsSphere[q].size(); i++)
            {
              std::cout<<dist_histogramBinsSphere[q][i]<<" ";
            }
          std::cout<<"];\n";
        }
      /*	std::cout<<"];\nD=[";
        	for(unsigned int i=0; i<directions.size(); i++) {
        	std::cout<<directions[i].transpose()<<"; ";
        	}
      */

    }
    else
      {
        std::cout<<"L: ";
        for(unsigned int i=0; i<histogramBinsLine.size(); i++)
          {
            std::cout<<histogramBinsLine[i]<<" ";
          }

        std::cout<<"\nF: ";
        for(unsigned int i=0; i<histogramBinsFlat.size(); i++)
          {
            std::cout<<histogramBinsFlat[i]<<" ";
          }

        std::cout<<"\nS: ";
        for(unsigned int i=0; i<histogramBinsSphere.size(); i++)
          {
            std::cout<<histogramBinsSphere[i]<<" ";
          }
        std::cout<<"\n";
      }
  }


  double NDTHistogram::getSimilarity(NDTHistogram &other)
  {

    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> T;

    this->bestFitToHistogram(other,T,false);

    return this->getSimilarity(other,T);
  }


  double NDTHistogram::getSimilarity(NDTHistogram &other, Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> &T){
    double score[3];
    double N_THIS[3], N_OTHER[3];
    for(unsigned int r = 0; r<3; r++){
      N_THIS[r] = 0;
      N_OTHER[r] = 0;
      score[r] = 0;
      for(int i=0; i<histogramBinsFlat.size(); i++)
        {
          N_THIS[r] += dist_histogramBinsFlat[r][i];
          N_OTHER[r] += other.dist_histogramBinsFlat[r][i];
        }
      for(int i=0; i<histogramBinsSphere.size(); i++)
        {
          N_THIS[r] += dist_histogramBinsSphere[r][i] ;
          N_OTHER[r] += other.dist_histogramBinsSphere[r][i] ;
        }
      N_THIS[r] += dist_histogramBinsLine[r][0];
      N_OTHER[r]+= other.dist_histogramBinsLine[r][0];
      N_THIS[r] = N_THIS[r]==0 ? INT_MAX : N_THIS[r];
      N_OTHER[r] = N_OTHER[r]==0 ? INT_MAX : N_OTHER[r];
    }

    for(int q = 0; q<averageDirections.size(); q++)
      {

        Eigen::Vector3d tr = T*averageDirections[q];
        if( this->histogramBinsFlat[q] == 0)
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
        if(!(idmin >=0 && idmin < histogramBinsFlat.size()))
          {
            continue;
          }

        for(unsigned int r = 0; r<3; r++)
          {
            score[r] += pow((double)this->dist_histogramBinsFlat[r][q]/N_THIS[r] - (double)other.dist_histogramBinsFlat[r][idmin]/N_OTHER[r],2);
          }

      }
    for(unsigned int r = 0; r<3; r++)
      {
        for(int i=0; i<histogramBinsSphere.size(); i++)
          {
            score[r] += pow( (double)this->dist_histogramBinsSphere[r][i]/N_THIS[r] - (double)other.dist_histogramBinsSphere[r][i]/N_OTHER[r] ,2);
          }

        score[r] += pow( (double)this->dist_histogramBinsLine[r][0]/N_THIS[r] - (double)other.dist_histogramBinsLine[r][0]/N_OTHER[r] ,2);
        double maxN, minN;
        maxN = (N_THIS[r] > N_OTHER[r]) ? N_THIS[r] : N_OTHER[r];
        minN = (N_THIS[r] < N_OTHER[r]) ? N_THIS[r] : N_OTHER[r];
        minN = (minN < 1) ? 1 : minN;

        score[r] = maxN*sqrt(score[r])/minN;
      }


    return score[0]+score[1];
  }

}
