//general
#include <PointCloudUtils.hh>
#include <fstream>

//ndt
#include <NDTMap.hh>
#include <OctTree.hh>
#include <NDTMatcherF2F.hh>

using namespace std;

pcl::PointCloud<pcl::PointXYZ> generateNegatives(pcl::PointCloud<pcl::PointXYZ> &cloud)
{

    pcl::PointCloud<pcl::PointXYZ> negCloud;
    double min_offset = 0.1;
    double max_offset = 2;
    pcl::PointXYZ origin, neg;

    origin.x = origin.y = origin.z = 0;
    for(int i =0; i<cloud.points.size(); i++)
    {
        double rand_offset = min_offset + (max_offset-min_offset)*
                             (double)rand()/(double)RAND_MAX;


        neg.x = (cloud.points[i].x-origin.x);
        neg.y = (cloud.points[i].y-origin.y);
        neg.z = (cloud.points[i].z-origin.z);
        double len = sqrt(neg.x*neg.x + neg.y*neg.y + neg.z*neg.z);
        double factor = (len - rand_offset) / len;
        factor = (factor < 0) ? 0 : factor;
        neg.x = factor*(cloud.points[i].x-origin.x) + origin.x;
        neg.y = factor*(cloud.points[i].y-origin.y) + origin.y;
        neg.z = factor*(cloud.points[i].z-origin.z) + origin.z;
        negCloud.points.push_back(neg);
    }
    return negCloud;
}

int main(int argc, char **argv)
{


    char filename[200];
    FILE *fout;
    //read in config file
    //read in config file name
    if(argc != 3)
    {
        cout<<"usage: ./representationTester logDirectory N_CLOUDS\n";
        return -1;
    }

    bool doMatch = true;
    bool doGroundTruth = false;

    pcl::PointCloud<pcl::PointXYZ> cloudLaser, negLaser;
    pcl::PointCloud<pcl::PointXYZ> cloudKinect, negKinect;
    pcl::PointCloud<pcl::PointXYZ> cloudSR, negSR;
    pcl::PointCloud<pcl::PointXYZ> cloudFotonic, negFotonic;
    pcl::PointCloud<pcl::PointXYZ> gt, negGt;

    string dirName = argv[1];
    string outputName;

    int N_CLOUDS = atoi(argv[2]);
    double LTHRESH = 0.0005;

    int tpKinect[N_CLOUDS], fpKinect[N_CLOUDS], tnKinect[N_CLOUDS], fnKinect[N_CLOUDS];
    int tpSR[N_CLOUDS], fpSR[N_CLOUDS], tnSR[N_CLOUDS], fnSR[N_CLOUDS];
    int tpFotonic[N_CLOUDS], fpFotonic[N_CLOUDS], tnFotonic[N_CLOUDS], fnFotonic[N_CLOUDS];
    int tpLaser[N_CLOUDS], fpLaser[N_CLOUDS], tnLaser[N_CLOUDS], fnLaser[N_CLOUDS];

    outputName = dirName+"results.m";
    ofstream foutResults((outputName).c_str());

    lslgeneric::OctTree::BIG_CELL_SIZE = 6.4;
    lslgeneric::OctTree::SMALL_CELL_SIZE = 0.4;
    lslgeneric::OctTree prototype;
    int tp=0,fp=0,tn=0,fn=0;

    for(int i=0; i<N_CLOUDS; i++)
    {

        char cloudname[200];
        snprintf(cloudname,199,"%s/laser/pc%04d.wrl",dirName.c_str(),i);
        cloudLaser = lslgeneric::readVRML(cloudname);

        snprintf(cloudname,199,"%s/kinect/pc%04d.wrl",dirName.c_str(),i);
        cloudKinect = lslgeneric::readVRML(cloudname);

        snprintf(cloudname,199,"%s/sr/pc%04d.wrl",dirName.c_str(),i);
        cloudSR = lslgeneric::readVRML(cloudname);

        snprintf(cloudname,199,"%s/fotonic/pc%04d.wrl",dirName.c_str(),i);
        cloudFotonic = lslgeneric::readVRML(cloudname);


        lslgeneric::NDTMap ndt(&prototype);

        if(doGroundTruth)
        {
            gt = cloudLaser;
            snprintf(cloudname,199,"%s/gt.wrl",dirName.c_str());
            cloudLaser= lslgeneric::readVRML(cloudname);
            ndt.loadPointCloud(cloudLaser);
        }

        if(doMatch)
        {
            Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> T;
            //lslgeneric::NDTMatcherF2F matcher;
            double __res[] = {0.2, 0.4, 1, 2};
            std::vector<double> resolutions (__res, __res+sizeof(__res)/sizeof(double));
            lslgeneric::NDTMatcherF2F matcher(false, false, false, resolutions);

            pcl::PointCloud<pcl::PointXYZ> tmp;
            if(!doGroundTruth)
            {
                for(int q=0; q<cloudLaser.points.size(); q++)
                {
                    pcl::PointXYZ pt = cloudLaser.points[q];
                    double angle_x = atan2(pt.x,pt.y);
                    double angle_y = atan2(pt.x,pt.z);
                    double angle_z = atan2(pt.y,pt.z);
                    //cout<<angle_x<<" "<<angle_y<<" "<<angle_z<<endl;
                    double dist = sqrt(pt.x*pt.x + pt.y*pt.y + pt.z*pt.z);
                    if(angle_y > -0.33 && dist<3.5) tmp.points.push_back(pt);
                }
                cloudLaser = tmp;
                T = Eigen::Translation<double,3>(-0.15,0.17,0.44)*
                    Eigen::AngleAxis<double>(0.06,Eigen::Vector3d::UnitX()) *
                    Eigen::AngleAxis<double>(0.1,Eigen::Vector3d::UnitY()) *
                    Eigen::AngleAxis<double>(-0.15,Eigen::Vector3d::UnitZ()) ;
                cout<<"Laser: "<<T.translation().transpose()<<" "<<T.rotation().eulerAngles(0,1,2).transpose()<<endl;
                lslgeneric::transformPointCloudInPlace(T,cloudLaser);
                ndt.loadPointCloud(cloudLaser);
            }


            tmp.points.clear();
            for(int q=0; q<cloudKinect.points.size(); q++)
            {
                pcl::PointXYZ pt = cloudKinect.points[q];
                double angle_x = atan2(pt.x,pt.y);
                double angle_y = atan2(pt.x,pt.z);
                double angle_z = atan2(pt.y,pt.z);
                double dist = sqrt(pt.x*pt.x + pt.y*pt.y + pt.z*pt.z);
                //cout<<angle_x<<" "<<angle_y<<" "<<angle_z<<endl;
                if(angle_y > -0.33 && angle_y < 0.33 && angle_z > -0.25 && dist<3.5) tmp.points.push_back(pt);
            }
            cloudKinect = tmp;

            T = Eigen::Translation<double,3>(-0.121,0.2278,0.514)*
                Eigen::AngleAxis<double>(0.066,Eigen::Vector3d::UnitX()) *
                Eigen::AngleAxis<double>(0.1098,Eigen::Vector3d::UnitY()) *
                Eigen::AngleAxis<double>(-0.0707,Eigen::Vector3d::UnitZ()) ;
            /*T = Eigen::Translation<double,3>(-0.07,0.058,0.44)*
            Eigen::AngleAxis<double>(0.07,Eigen::Vector3d::UnitX()) *
            Eigen::AngleAxis<double>(0.09,Eigen::Vector3d::UnitY()) *
            Eigen::AngleAxis<double>(-0.08,Eigen::Vector3d::UnitZ()) ;*/
            lslgeneric::transformPointCloudInPlace(T,cloudKinect);
            matcher.match(cloudLaser, cloudKinect, T);
            cout<<"Kinect: "<<T.translation().transpose()<<" "<<T.rotation().eulerAngles(0,1,2).transpose()<<endl;
            lslgeneric::transformPointCloudInPlace(T,cloudKinect);
            snprintf(cloudname,199,"%s/ndt/rk%04d.wrl",dirName.c_str(),i);
            fout = fopen(cloudname,"w");
            fprintf(fout,"#VRML V2.0 utf8\n");
            lslgeneric::writeToVRML(fout,cloudLaser,Eigen::Vector3d(1,1,1));
            lslgeneric::writeToVRML(fout,cloudKinect,Eigen::Vector3d(1,0,0));
            fclose(fout);

            tmp.points.clear();
            for(int q=0; q<cloudSR.points.size(); q++)
            {
                pcl::PointXYZ pt = cloudSR.points[q];
                double angle_x = atan2(pt.x,pt.y);
                double angle_y = atan2(pt.x,pt.z);
                double angle_z = atan2(pt.y,pt.z);
                double dist = sqrt(pt.x*pt.x + pt.y*pt.y + pt.z*pt.z);
                //cout<<angle_x<<" "<<angle_y<<" "<<angle_z<<endl;
                if(angle_y > -0.33 && angle_y < 0.33 && angle_z > -0.25 && dist<3.5) tmp.points.push_back(pt);
            }
            cloudSR = tmp;
            T = Eigen::Translation<double,3>(0.052,0.036,0.42)*
                Eigen::AngleAxis<double>(0.067,Eigen::Vector3d::UnitX()) *
                Eigen::AngleAxis<double>(0.061,Eigen::Vector3d::UnitY()) *
                Eigen::AngleAxis<double>(-0.0859,Eigen::Vector3d::UnitZ()) ;
            lslgeneric::transformPointCloudInPlace(T,cloudSR);
            matcher.match(cloudLaser, cloudSR, T);
            cout<<"SR: "<<T.translation().transpose()<<" "<<T.rotation().eulerAngles(0,1,2).transpose()<<endl;
            lslgeneric::transformPointCloudInPlace(T,cloudSR);
            snprintf(cloudname,199,"%s/ndt/rs%04d.wrl",dirName.c_str(),i);
            fout = fopen(cloudname,"w");
            fprintf(fout,"#VRML V2.0 utf8\n");
            lslgeneric::writeToVRML(fout,cloudLaser,Eigen::Vector3d(1,1,1));
            lslgeneric::writeToVRML(fout,cloudSR,Eigen::Vector3d(1,0,0));
            fclose(fout);

            tmp.points.clear();
            for(int q=0; q<cloudFotonic.points.size(); q++)
            {
                pcl::PointXYZ pt = cloudFotonic.points[q];
                double angle_x = atan2(pt.x,pt.y);
                double angle_y = atan2(pt.x,pt.z);
                double angle_z = atan2(pt.y,pt.z);
                //cout<<angle_x<<" "<<angle_y<<" "<<angle_z<<endl;
                double dist = sqrt(pt.x*pt.x + pt.y*pt.y + pt.z*pt.z);
                if(angle_y > -0.33 && angle_y < 0.33 && angle_z > -0.25 && dist<3.5) tmp.points.push_back(pt);
            }
            cloudFotonic = tmp;
            T = Eigen::Translation<double,3>(-0.0988,0.045,0.452)*
                Eigen::AngleAxis<double>(0.0952,Eigen::Vector3d::UnitX()) *
                Eigen::AngleAxis<double>(0.0718,Eigen::Vector3d::UnitY()) *
                Eigen::AngleAxis<double>(-0.0758,Eigen::Vector3d::UnitZ()) ;
            lslgeneric::transformPointCloudInPlace(T,cloudFotonic);
            matcher.match(cloudLaser, cloudFotonic, T);
            cout<<"Fotonic: "<<T.translation().transpose()<<" "<<T.rotation().eulerAngles(0,1,2).transpose()<<endl;
            lslgeneric::transformPointCloudInPlace(T,cloudFotonic);
            snprintf(cloudname,199,"%s/ndt/rf%04d.wrl",dirName.c_str(),i);
            fout = fopen(cloudname,"w");
            fprintf(fout,"#VRML V2.0 utf8\n");
            lslgeneric::writeToVRML(fout,cloudLaser,Eigen::Vector3d(1,1,1));
            lslgeneric::writeToVRML(fout,cloudFotonic,Eigen::Vector3d(1,0,0));
            fclose(fout);


            if(doGroundTruth)
            {
                pcl::PointCloud<pcl::PointXYZ> gt2;
                for(int q=0; q<gt.points.size(); q++)
                {
                    pcl::PointXYZ pt = gt.points[q];
                    double angle_x = atan2(pt.x,pt.y);
                    double angle_y = atan2(pt.x,pt.z);
                    double angle_z = atan2(pt.y,pt.z);
                    //cout<<angle_x<<" "<<angle_y<<" "<<angle_z<<endl;
                    double dist = sqrt(pt.x*pt.x + pt.y*pt.y + pt.z*pt.z);
                    if(angle_y > -0.33 && dist <3.5) gt2.points.push_back(pt);
                }
                gt = gt2;
                T = Eigen::Translation<double,3>(-0.15,0.17,0.44)*
                    Eigen::AngleAxis<double>(0.06,Eigen::Vector3d::UnitX()) *
                    Eigen::AngleAxis<double>(0.1,Eigen::Vector3d::UnitY()) *
                    Eigen::AngleAxis<double>(-0.15,Eigen::Vector3d::UnitZ()) ;
                //lslgeneric::transformPointCloudInPlace(T,gt);
                //matcher.match(cloudLaser, gt, T);
                cout<<"Laser: "<<T.translation().transpose()<<" "<<T.rotation().eulerAngles(0,1,2).transpose()<<endl;
                lslgeneric::transformPointCloudInPlace(T,gt);
                snprintf(cloudname,199,"%s/ndt/rl%04d.wrl",dirName.c_str(),i);
                fout = fopen(cloudname,"w");
                fprintf(fout,"#VRML V2.0 utf8\n");
                lslgeneric::writeToVRML(fout,cloudLaser,Eigen::Vector3d(1,1,1));
                lslgeneric::writeToVRML(fout,gt,Eigen::Vector3d(1,0,0));
                fclose(fout);
                cloudLaser = gt;
            }
        }

        ndt.computeNDTCells();
        snprintf(filename,199,"%s/ndt/ndtTree%d.wrl",dirName.c_str(),i);
        ndt.writeToVRML(filename);

        negLaser = generateNegatives(cloudLaser);
        snprintf(cloudname,199,"%s/laser/negpc%04d.wrl",dirName.c_str(),i);
        lslgeneric::writeToVRML(cloudname,negLaser);
        negKinect = generateNegatives(cloudKinect);
        snprintf(cloudname,199,"%s/kinect/negpc%04d.wrl",dirName.c_str(),i);
        lslgeneric::writeToVRML(cloudname,negKinect);
        negSR = generateNegatives(cloudSR);
        snprintf(cloudname,199,"%s/sr/negpc%04d.wrl",dirName.c_str(),i);
        lslgeneric::writeToVRML(cloudname,negSR);
        negFotonic = generateNegatives(cloudFotonic);
        snprintf(cloudname,199,"%s/fotonic/negpc%04d.wrl",dirName.c_str(),i);
        lslgeneric::writeToVRML(cloudname,negFotonic);

        /*	tp=0;fp=0;tn=0;fn=0;
        	lslgeneric::NDTMap ndtLaser(&prototype);
        	if(doGroundTruth) {
        	    ndtLaser.loadPointCloud(gt);
        	} else {
        	    ndtLaser.loadPointCloud(cloudLaser);
        	}
        	ndtLaser.computeNDTCells();
        	for(int j=0; j<cloudLaser.points.size(); j++) {
        	    double prob = ndtLaser.getLikelihoodForPoint(cloudLaser.points[j]);
        	    double nprob = ndtLaser.getLikelihoodForPoint(negLaser.points[j]);
        	    (prob > LTHRESH) ? tp++ : fp++;
        	    (nprob > LTHRESH) ? fn++ : tn++;
        	}
        	cout<<"Laser at "<<i<<" : "<<tp<<" "<<fp<<" "<<tn<<" "<<fn<<endl;
        	tpLaser[i] = tp; fpLaser[i] = fp; tnLaser[i] = tn; fnLaser[i] = fn;

        	tp=0;fp=0;tn=0;fn=0;
        	lslgeneric::NDTMap ndtKinect(&prototype);
        	ndtKinect.loadPointCloud(cloudKinect);
        	ndtKinect.computeNDTCells();
        	for(int j=0; j<cloudLaser.points.size(); j++) {
        	    double prob = ndtKinect.getLikelihoodForPoint(cloudLaser.points[j]);
        	    double nprob = ndtKinect.getLikelihoodForPoint(negLaser.points[j]);
        	    (prob > LTHRESH) ? tp++ : fp++;
        	    (nprob > LTHRESH) ? fn++ : tn++;
        	}
        	cout<<"Kinect at "<<i<<" : "<<tp<<" "<<fp<<" "<<tn<<" "<<fn<<endl;
        	tpKinect[i] = tp; fpKinect[i] = fp; tnKinect[i] = tn; fnKinect[i] = fn;

        	tp=0;fp=0;tn=0;fn=0;
        	lslgeneric::NDTMap ndtSR(&prototype);
        	ndtSR.loadPointCloud(cloudSR);
        	ndtSR.computeNDTCells();
        	for(int j=0; j<cloudLaser.points.size(); j++) {
        	    double prob = ndtSR.getLikelihoodForPoint(cloudLaser.points[j]);
        	    double nprob = ndtSR.getLikelihoodForPoint(negLaser.points[j]);
        	    (prob > LTHRESH) ? tp++ : fp++;
        	    (nprob > LTHRESH) ? fn++ : tn++;
        	}
        	cout<<"SR at "<<i<<" : "<<tp<<" "<<fp<<" "<<tn<<" "<<fn<<endl;
        	tpSR[i] = tp; fpSR[i] = fp; tnSR[i] = tn; fnSR[i] = fn;
        	tp=0;fp=0;tn=0;fn=0;

        	lslgeneric::NDTMap ndtFotonic(&prototype);
        	ndtFotonic.loadPointCloud(cloudFotonic);
        	ndtFotonic.computeNDTCells();
        	for(int j=0; j<cloudLaser.points.size(); j++) {
        	    double prob = ndtFotonic.getLikelihoodForPoint(cloudLaser.points[j]);
        	    double nprob = ndtFotonic.getLikelihoodForPoint(negLaser.points[j]);
        	    (prob > LTHRESH) ? tp++ : fp++;
        	    (nprob > LTHRESH) ? fn++ : tn++;
        	}
        	cout<<"Fotonic at "<<i<<" : "<<tp<<" "<<fp<<" "<<tn<<" "<<fn<<endl;
        	tpFotonic[i] = tp; fpFotonic[i] = fp; tnFotonic[i] = tn; fnFotonic[i] = fn;
        */
        tp=0;
        fp=0;
        tn=0;
        fn=0;
        for(int j=0; j<cloudLaser.points.size(); j++)
        {
            double prob = ndt.getLikelihoodForPoint(cloudLaser.points[j]);
            double nprob = ndt.getLikelihoodForPoint(negLaser.points[j]);
            (prob > LTHRESH) ? tp++ : fp++;
            (nprob > LTHRESH) ? fn++ : tn++;
        }
        cout<<"Laser at "<<i<<" : "<<tp<<" "<<fp<<" "<<tn<<" "<<fn<<endl;
        tpLaser[i] = tp;
        fpLaser[i] = fp;
        tnLaser[i] = tn;
        fnLaser[i] = fn;


        tp=0;
        fp=0;
        tn=0;
        fn=0;
        for(int j=0; j<cloudKinect.points.size(); j++)
        {
            double prob = ndt.getLikelihoodForPoint(cloudKinect.points[j]);
            double nprob = ndt.getLikelihoodForPoint(negKinect.points[j]);
            (prob > LTHRESH) ? tp++ : fp++;
            (nprob > LTHRESH) ? fn++ : tn++;
        }
        cout<<"Kinect at "<<i<<" : "<<tp<<" "<<fp<<" "<<tn<<" "<<fn<<endl;
        tpKinect[i] = tp;
        fpKinect[i] = fp;
        tnKinect[i] = tn;
        fnKinect[i] = fn;


        tp=0;
        fp=0;
        tn=0;
        fn=0;
        for(int j=0; j<cloudSR.points.size(); j++)
        {
            double prob = ndt.getLikelihoodForPoint(cloudSR.points[j]);
            double nprob = ndt.getLikelihoodForPoint(negSR.points[j]);
            (prob > LTHRESH) ? tp++ : fp++;
            (nprob > LTHRESH) ? fn++ : tn++;
        }
        cout<<"SR at "<<i<<" : "<<tp<<" "<<fp<<" "<<tn<<" "<<fn<<endl;
        tpSR[i] = tp;
        fpSR[i] = fp;
        tnSR[i] = tn;
        fnSR[i] = fn;

        tp=0;
        fp=0;
        tn=0;
        fn=0;
        for(int j=0; j<cloudFotonic.points.size(); j++)
        {
            double prob = ndt.getLikelihoodForPoint(cloudFotonic.points[j]);
            double nprob = ndt.getLikelihoodForPoint(negFotonic.points[j]);
            (prob > LTHRESH) ? tp++ : fp++;
            (nprob > LTHRESH) ? fn++ : tn++;
        }
        cout<<"Fotonic at "<<i<<" : "<<tp<<" "<<fp<<" "<<tn<<" "<<fn<<endl;
        tpFotonic[i] = tp;
        fpFotonic[i] = fp;
        tnFotonic[i] = tn;
        fnFotonic[i] = fn;

    }

    cout<<"Done, now writing results for "<<N_CLOUDS<<" clouds at "<<outputName<<"...\n";
    //laser
    foutResults<<"tpLaser = [";
    for(int i=0; i<N_CLOUDS; i++)
    {
        foutResults<<(double)tpLaser[i]/(double)(tpLaser[i]+fnLaser[i])<<" ";
    }
    foutResults<<"];\nfpLaser = [";
    for(int i=0; i<N_CLOUDS; i++)
    {
        foutResults<<(double)fpLaser[i]/(double)(tnLaser[i]+fpLaser[i])<<" ";
    }
    foutResults<<"];\ntnLaser = [";
    for(int i=0; i<N_CLOUDS; i++)
    {
        foutResults<<(double)tnLaser[i]/(double)(tpLaser[i]+fnLaser[i])<<" ";
    }
    foutResults<<"];\nfnLaser = [";
    for(int i=0; i<N_CLOUDS; i++)
    {
        foutResults<<(double)fnLaser[i]/(double)(tnLaser[i]+fpLaser[i])<<" ";
    }
    foutResults<<"];\n";

    //kinect
    foutResults<<"tpKinect = [";
    for(int i=0; i<N_CLOUDS; i++)
    {
        foutResults<<(double)tpKinect[i]/(double)(tpKinect[i]+fnKinect[i])<<" ";
    }
    foutResults<<"];\nfpKinect = [";
    for(int i=0; i<N_CLOUDS; i++)
    {
        foutResults<<(double)fpKinect[i]/(double)(tnKinect[i]+fpKinect[i])<<" ";
    }
    foutResults<<"];\ntnKinect = [";
    for(int i=0; i<N_CLOUDS; i++)
    {
        foutResults<<(double)tnKinect[i]/(double)(tpKinect[i]+fnKinect[i])<<" ";
    }
    foutResults<<"];\nfnKinect = [";
    for(int i=0; i<N_CLOUDS; i++)
    {
        foutResults<<(double)fnKinect[i]/(double)(tnKinect[i]+fpKinect[i])<<" ";
    }
    foutResults<<"];\n";

    //sr
    foutResults<<"tpSR = [";
    for(int i=0; i<N_CLOUDS; i++)
    {
        foutResults<<(double)tpSR[i]/(double)(tpSR[i]+fnSR[i])<<" ";
    }
    foutResults<<"];\nfpSR = [";
    for(int i=0; i<N_CLOUDS; i++)
    {
        foutResults<<(double)fpSR[i]/(double)(tnSR[i]+fpSR[i])<<" ";
    }
    foutResults<<"];\ntnSR = [";
    for(int i=0; i<N_CLOUDS; i++)
    {
        foutResults<<(double)tnSR[i]/(double)(tnSR[i]+fpSR[i])<<" ";
    }
    foutResults<<"];\nfnSR = [";
    for(int i=0; i<N_CLOUDS; i++)
    {
        foutResults<<(double)fnSR[i]/(double)(tpSR[i]+fnSR[i])<<" ";
    }
    foutResults<<"];\n";

    //fotonic
    foutResults<<"tpFotonic = [";
    for(int i=0; i<N_CLOUDS; i++)
    {
        foutResults<<(double)tpFotonic[i]/(double)(tpFotonic[i]+fnFotonic[i])<<" ";
    }
    foutResults<<"];\nfpFotonic = [";
    for(int i=0; i<N_CLOUDS; i++)
    {
        foutResults<<(double)fpFotonic[i]/(double)(tnFotonic[i]+fpFotonic[i])<<" ";
    }
    foutResults<<"];\ntnFotonic = [";
    for(int i=0; i<N_CLOUDS; i++)
    {
        foutResults<<(double)tnFotonic[i]/(double)(tnFotonic[i]+fpFotonic[i])<<" ";
    }
    foutResults<<"];\nfnFotonic = [";
    for(int i=0; i<N_CLOUDS; i++)
    {
        foutResults<<(double)fnFotonic[i]/(double)(tpFotonic[i]+fnFotonic[i])<<" ";
    }
    foutResults<<"];\n";

    return 0;

}


