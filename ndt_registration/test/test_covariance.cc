#include <NDTMatcher.hh>
#include <NDTMatcherF2F.hh>
#include <NDTMap.hh>
#include <OctTree.hh>
#include <AdaptiveOctTree.hh>
#include "ros/ros.h"
#include "pcl/point_cloud.h"
#include "sensor_msgs/PointCloud2.h"
#include "pcl/io/pcd_io.h"
#include "pcl/features/feature.h"
#include "pcl/registration/icp.h"
#include "pcl/filters/voxel_grid.h"

#include <cstdio>
#include <Eigen/Eigen>
#include <PointCloudUtils.hh>
#include <fstream>

#define DEBUG_COVARIANCE 0

using namespace std;

float ranf()           /* ranf() is uniform in 0..1 */
{
    return (float)rand()/(float)RAND_MAX;
}

float box_muller(float m, float s)  /* normal random variate generator */
{
    /* mean m, standard deviation s */
    float x1, x2, w, y1;
    static float y2;
    static int use_last = 0;

    if (use_last)	        /* use value from previous call */
    {
        y1 = y2;
        use_last = 0;
    }
    else
    {
        do
        {
            x1 = 2.0 * ranf() - 1.0;
            x2 = 2.0 * ranf() - 1.0;
            w = x1 * x1 + x2 * x2;
        }
        while ( w >= 1.0 );

        w = sqrt( (-2.0 * log( w ) ) / w );
        y1 = x1 * w;
        y2 = x2 * w;
        use_last = 1;
    }

    return( m + y1 * s );
}


bool matchICP(pcl::PointCloud<pcl::PointXYZ> &fixed,  pcl::PointCloud<pcl::PointXYZ> &moving,
              Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> &Tout)
{

    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_in (new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_out (new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::ConstPtr f (new pcl::PointCloud<pcl::PointXYZ>(fixed) );
    pcl::PointCloud<pcl::PointXYZ>::ConstPtr m (new pcl::PointCloud<pcl::PointXYZ>(moving) );

    pcl::VoxelGrid<pcl::PointXYZ> gr1,gr2;
    gr1.setLeafSize(0.1,0.1,0.1);
    gr2.setLeafSize(0.1,0.1,0.1);

    gr1.setInputCloud(m);
    gr2.setInputCloud(f);

    cloud_in->height = 1;
    cloud_in->width = cloud_in->points.size();
    cloud_out->height = 1;
    cloud_out->width = cloud_out->points.size();
    cloud_in->is_dense = false;
    cloud_out->is_dense = false;

    gr1.filter(*cloud_in);
    gr2.filter(*cloud_out);

    pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;

    icp.setMaximumIterations(1000);
    //cout<<"max itr are "<<icp.getMaximumIterations()<<endl;
    icp.setInputCloud(cloud_in);
    icp.setInputTarget(cloud_out);

    icp.setRANSACOutlierRejectionThreshold (2);
    icp.setMaxCorrespondenceDistance(10);
    icp.setTransformationEpsilon(0.00001);
//    cout<<"ransac outlier thersh   : "<<icp.getRANSACOutlierRejectionThreshold ()<<endl;
//    cout<<"correspondance max dist : "<<icp.getMaxCorrespondenceDistance() << endl;
//    cout<<"epsilon : "<<icp.getTransformationEpsilon() << endl;
    pcl::PointCloud<pcl::PointXYZ> Final;
    icp.align(Final);


//    std::cout << "has converged:" << icp.hasConverged() << " score: " <<
//	icp.getFitnessScore() << std::endl;
//    std::cout << icp.getFinalTransformation() << std::endl;

    //Eigen::Transform<float,3,Eigen::Affine,Eigen::ColMajor> tTemp;
    Tout = (icp.getFinalTransformation()).cast<double>();

    /*    char fname[50];
        snprintf(fname,49,"/home/tsv/ndt_tmp/c2_offset.wrl");
        FILE *fout = fopen(fname,"w");
        fprintf(fout,"#VRML V2.0 utf8\n");
        lslgeneric::writeToVRML(fout,*cloud_out,Eigen::Vector3d(0,1,0));
        lslgeneric::writeToVRML(fout,Final,Eigen::Vector3d(1,0,0));
        lslgeneric::writeToVRML(fout,*cloud_in,Eigen::Vector3d(1,1,1));
        fclose(fout);
    */
    return icp.hasConverged();

}

int
main (int argc, char** argv)
{
    std::ofstream logger ("/home/tsv/ndt_tmp/covariance.m");
    cout.precision(15);

    pcl::PointCloud<pcl::PointXYZ> prev, curr, currHere, currHere2;
    char prevName[500], currName[500];
    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> prevPose, currPose;

    double __res[] = {0.2, 0.4, 1, 2};
    std::vector<double> resolutions (__res, __res+sizeof(__res)/sizeof(double));
    lslgeneric::NDTMatcherF2F matcherF2F(false, false, false, resolutions);
    lslgeneric::NDTMatcher matcherP2F;

    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> Pr, R;
    struct timeval tv_now;

    gettimeofday(&tv_now, NULL);
    srand(tv_now.tv_usec);

    int N_SAMPLES = 100;

    if(argc != 2)
    {
        std::cout<<"Usage: "<<argv[0]<<" configFile\n";
        return -1;
    }

    FILE *fin = fopen(argv[1],"r");
    double xd,yd,zd, x,y,z,w, ts;
    string prefix;
    char *line = NULL;
    size_t len;
    bool first = true;
    double randX, randY, randZ, randRoll, randPitch, randYaw;
    double dev_x = 0.3 ,dev_y=0.3 ,dev_z =0.3 ,dev_roll =0.1,dev_pitch =0.1,dev_yaw=0.1;

    int ctr = 0;

    //get first line
    int n = getline(&line,&len,fin);
    if(n <=0 ) return -1;
    prefix = line;
    *(prefix.rbegin()) = '\0';

    while(getline(&line,&len,fin) > 0)
    {

        int n = sscanf(line,"%lf %lf %lf %lf %lf %lf %lf %lf",
                       &ts,&xd,&yd,&zd,&x,&y,&z,&w);
        if(n != 8)
        {
            cout<<"wrong format of pose at : "<<line<<endl;
            break;
        }

        if(first)
        {
            //set previous pose to curent pose
            prevPose =Eigen::Translation<double,3>(xd,yd,zd)*
                      Eigen::Quaternion<double>(w,x,y,z);
            first = false;
            ctr++;
            continue;
        }


        currPose =  Eigen::Translation<double,3>(xd,yd,zd)*
                    Eigen::Quaternion<double>(w,x,y,z);

        /*
        r = M_PI*r/180;
        p = M_PI*p/180;
            //due to the way we collected data, offset by 90, shouldn't matter
        y = y+90;
        y = M_PI*y/180;

        if(first) {
            //set previous pose to curent pose
            prevPose =Eigen::Translation<double,3>(xd,yd,zd)*
        	      Eigen::AngleAxis<double>(r,Eigen::Vector3d::UnitX())*
        	      Eigen::AngleAxis<double>(p,Eigen::Vector3d::UnitY())*
        	      Eigen::AngleAxis<double>(y,Eigen::Vector3d::UnitZ());
            first = false;
            ctr++;
            continue;
        }


        currPose =Eigen::Translation<double,3>(xd,yd,zd)*
            Eigen::AngleAxis<double>(r,Eigen::Vector3d::UnitX())*
            Eigen::AngleAxis<double>(p,Eigen::Vector3d::UnitY())*
            Eigen::AngleAxis<double>(y,Eigen::Vector3d::UnitZ());
        */


        //compute ground truth relative pose
        Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> P = prevPose.inverse()*currPose;

        cout<<"testing scan "<<ctr-1<<" to "<<ctr<<endl;
#if DEBUG_COVARIANCE
        cout<<"prev pose    : "<<prevPose.translation().transpose()<<" rpy "<<prevPose.rotation().eulerAngles(0,1,2).transpose()<<endl;
        cout<<"curr pose    : "<<currPose.translation().transpose()<<" rpy "<<r<<" "<<p<<" "<<y<<endl;
        cout<<"gt difference: "<<P.translation().transpose()<<" rpy "<<P.rotation().eulerAngles(0,1,2).transpose()<<endl;
#endif

        for( int i=0; i<N_SAMPLES; i++)
        {
            cout<<"itr: "<<i<<endl;
            Pr.setIdentity();
            snprintf(prevName,499,"%s%03d_%03d.wrl",prefix.c_str(),ctr-1,i);
            snprintf(currName,499,"%s%03d_%03d.wrl",prefix.c_str(),ctr,i);
            prev = lslgeneric::readVRML(prevName);
            curr = lslgeneric::readVRML(currName);

            //draw random numbers
            randX = box_muller(0,dev_x);
            randY = box_muller(0,dev_y);
            randZ = box_muller(0,dev_z);
            randRoll = box_muller(0,dev_roll);
            randPitch = box_muller(0,dev_pitch);
            randYaw = box_muller(0,dev_yaw);

            currPose =  Eigen::Translation<double,3>(xd+randX,yd+randY,zd+randZ)*
                        Eigen::Quaternion<double>(w,x,y,z)*
                        Eigen::AngleAxis<double>(randRoll,Eigen::Vector3d::UnitX())*
                        Eigen::AngleAxis<double>(randPitch,Eigen::Vector3d::UnitY())*
                        Eigen::AngleAxis<double>(randYaw,Eigen::Vector3d::UnitZ());
            //compute current pose
            /*currPose =Eigen::Translation<double,3>(xd+randX,yd+randY,zd+randZ)*
            Eigen::AngleAxis<double>(r+randRoll,Eigen::Vector3d::UnitX())*
            Eigen::AngleAxis<double>(p+randPitch,Eigen::Vector3d::UnitY())*
            Eigen::AngleAxis<double>(y+randYaw,Eigen::Vector3d::UnitZ());
            */

            //compute relative pose
            Pr = prevPose.inverse()*currPose;

#if DEBUG_COVARIANCE
            cout<<"curr pose: "<<currPose.translation().transpose()<<" "<<r+randRoll<<" "<<p+randPitch<<" "<<y+randYaw<<endl;
            cout<<"diff: "<<Pr.translation().transpose()<<" rpy "<<Pr.rotation().eulerAngles(0,1,2).transpose()<<endl;
#endif
            bool ret;
            Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> offset;
            /*
            	    //point2NDT
            	    //transform curr scan by Pr
            	    currHere = lslgeneric::transformPointCloud(Pr,curr);
            	    ret = matcherP2F.match(prev,currHere,R);
            	    //compute P[-](Pr[+]R)
            	    offset = P.inverse()*(R*Pr);
            	    //that's the error! log it.
            	    logger<<"pnt2ndt("<<ctr<<","<<i+1<<",:) = ["<<offset.translation().transpose()<<" "<<offset.rotation().eulerAngles(0,1,2).transpose()<<"];\n";

            	    //ICP
            	    //transform curr scan by Pr
            	    currHere = lslgeneric::transformPointCloud(Pr,curr);
            	    ret = matchICP(prev,currHere,R);
            	    //compute P[-](Pr[+]R)
            	    offset = P.inverse()*(R*Pr);
            	    //that's the error! log it.
            	    logger<<"icp("<<ctr<<","<<i+1<<",:) = ["<<offset.translation().transpose()<<" "<<offset.rotation().eulerAngles(0,1,2).transpose()<<"];\n";
            */
            //NDT2NDT
            //transform curr scan by Pr
            currHere = lslgeneric::transformPointCloud(Pr,curr);
            ret = matcherF2F.match(prev,currHere,R);
#if DEBUG_COVARIANCE
            currHere2 = lslgeneric::transformPointCloud(R,currHere);
            char fname[50];
            snprintf(fname,49,"/home/tsv/ndt_tmp/c_%03d.wrl",ctr);
            FILE *fout = fopen(fname,"w");
            fprintf(fout,"#VRML V2.0 utf8\n");
            lslgeneric::writeToVRML(fout,currHere2,Eigen::Vector3d(0,1,0));
            lslgeneric::writeToVRML(fout,currHere,Eigen::Vector3d(1,0,0));
            lslgeneric::writeToVRML(fout,prev,Eigen::Vector3d(1,1,1));
            fclose(fout);
            cout<<"registration: "<<R.translation().transpose()<<" rpy "<<R.rotation().eulerAngles(0,1,2).transpose()<<endl;
#endif
            //compute P[-](Pr[+]R)
            offset = P.inverse()*(R*Pr);
            //that's the error! log it.
            logger<<"ndt2ndt("<<ctr<<","<<i+1<<",:) = ["<<offset.translation().transpose()<<" "<<offset.rotation().eulerAngles(0,1,2).transpose()<<"];\n";
//	    cout<<offset.translation().transpose()<<" "<<offset.rotation().eulerAngles(0,1,2).transpose()<<endl;

        }

        Eigen::Matrix<double,6,6> cov;
//	matcherP2F.covariance(prev,curr,P,cov);
//	logger<<"COVpnt2ndt("<<ctr<<",:,:) = ["<<cov.row(0)<<";\n"<<cov.row(1)<<";\n"<<cov.row(2)<<";\n"<<cov.row(3)<<";\n"<<cov.row(4)<<";\n"<<cov.row(5)<<"];\n";
        matcherF2F.covariance(prev,curr,P,cov);
        logger<<"COVndt2ndt("<<ctr<<",:,:) = ["<<cov.row(0)<<";\n"<<cov.row(1)<<";\n"<<cov.row(2)<<";\n"<<cov.row(3)<<";\n"<<cov.row(4)<<";\n"<<cov.row(5)<<"];\n";

        prevPose =Eigen::Translation<double,3>(xd,yd,zd)*
                  Eigen::Quaternion<double>(w,x,y,z);
        ctr++;
    }

}


