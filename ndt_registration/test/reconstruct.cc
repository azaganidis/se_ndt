#include <NDTMatcherF2F.hh>
#include <NDTMatcher.hh>
#include <NDTMap.hh>
#include <OctTree.hh>
#include <AdaptiveOctTree.hh>
#include "ros/ros.h"
#include "pcl/point_cloud.h"
#include "sensor_msgs/PointCloud2.h"

#include "pcl/io/pcd_io.h"
#include <pcl/filters/radius_outlier_removal.h>

#include "pcl/features/feature.h"
#include <cstdio>
#include <Eigen/Eigen>
#include <PointCloudUtils.hh>
#include <fstream>
#include <NDTMap.hh>
#include <LazzyGrid.hh>

using namespace std;
using namespace lslgeneric;

void filterDensity(pcl::PointCloud<pcl::PointXYZ> &rawCloud, pcl::PointCloud<pcl::PointXYZ> &pc)
{
    pcl::RadiusOutlierRemoval<pcl::PointXYZ> filter;
    filter.setRadiusSearch(0.2);
    filter.setMinNeighborsInRadius(3);
    pcl::PointCloud<pcl::PointXYZ>::Ptr pcptr(new pcl::PointCloud<pcl::PointXYZ>());
    *pcptr = rawCloud;
    filter.setInputCloud(pcptr);
    filter.filter(pc);
}

int
main (int argc, char** argv)
{

    pcl::PointCloud<pcl::PointXYZ> prev, thisone1, thisone, cloudFinal;
    double []__res = {0.2, 0.4, 1, 2};
    std::vector<double> resolutions (__res, __res+sizeof(__res)/sizeof(double));
    lslgeneric::NDTMatcherF2F matcher(false, false, false, resolutions);
    //lslgeneric::NDTMatcher matcher;
    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> Tprev, Tloc, Treg, Tglob;

    pcl::PCDReader reader;
    Tprev.setIdentity();
    Tloc.setIdentity();
    Treg.setIdentity();
    Tglob.setIdentity();

    //we do a single scan to scan registration
    if(argc != 2)
    {
        cout<<"usage: ./reconstruct configFile\n";
        return -1;
    }

    FILE *fin = fopen(argv[1],"r");
    double xd,yd,zd, q1,q2,q3,q4;
    double xoffset=0, yoffset=0, zoffset=-0.01;
    string pcName;
    char *line = NULL;
    size_t len;
    bool first = true;

    while(getline(&line,&len,fin) > 0)
    {

        int n = sscanf(line,"%lf,%lf,%lf %lf",
                       &q1,&q2,&q3,&q4);
        if(n != 4)
        {
            cout<<"wrong format of pose at : "<<line<<endl;
            break;
        }
        if(!getline(&line,&len,fin) > 0) break;
        n = sscanf(line,"%lf,%lf,%lf",
                   &xd,&yd,&zd);
        if(n != 3)
        {
            cout<<"wrong format of pose at : "<<line<<endl;
            break;
        }

        if(!getline(&line,&len,fin) > 0) break;
        pcName = line;
        *(pcName.rbegin()) = '\0';

//	thisone = lslgeneric::readVRML(pcName.c_str());


//	thisone1.width = thisone1.height = 0;
        cout<<"reading "<<pcName<<endl;
        reader.read<pcl::PointXYZ>(pcName,thisone);
        //cout<<"filtering density..."<<thisone1.points.size()<<endl;
        //filterDensity(thisone1,thisone);
        cout<<" --> "<<thisone.points.size()<<endl;
        lslgeneric::writeToVRML("/home/tsv/ndt_tmp/last.wrl",thisone);


        Tloc = Tprev.inverse();
        Tprev = Eigen::Translation<double,3>(xd,yd,zd)*Eigen::AngleAxis<double>(q4,Eigen::Vector3d(q1,q2,q3));
        Tloc = Tloc*Tprev;

        transformPointCloudInPlace(Tloc,thisone);
        cout<<"old pose is t:"<<Tprev.translation().transpose()<<" r "<<Tprev.rotation().eulerAngles(0,1,2).transpose()<<endl;
        cout<<"local pose is t:"<<Tloc.translation().transpose()<<" r "<<Tloc.rotation().eulerAngles(0,1,2).transpose()<<endl;

        if(!first)
        {
            //register
            Treg.setIdentity();
            matcher.match(thisone,prev,Treg);
            cout<<"registration pose is t:"<<Treg.translation().transpose()<<" r "<<Treg.rotation().eulerAngles(0,1,2).transpose()<<endl;
            Tglob = Tglob*Tloc*Treg;
            cout<<"new global pose t:"<<Tglob.translation().transpose()<<" r "<<Tglob.rotation().eulerAngles(0,1,2).transpose()<<endl;
        }
        else
        {
            Tglob = Tloc;
            first = false;
        }

        //set prev to thisone
        prev = lslgeneric::readVRML(pcName.c_str());
        thisone = prev;
        //transform thisone
        transformPointCloudInPlace(Tglob,thisone);

        cout<<"read vrml file at "<<pcName<<endl;
        cloudFinal += thisone;
    }

    fclose(fin);
    lslgeneric::writeToVRML("/home/tsv/ndt_tmp/final.wrl",cloudFinal);
    lslgeneric::writeToVRML("/home/tsv/ndt_tmp/last.wrl",thisone);

#if 0
    //completely unrelated, needed it for other debug TSV
    pcl::PointCloud<pcl::PointXYZ> cl;
    pcl::PointXYZ pt;
    for(double x = 0; x<2; x+=0.1)
    {
        for(double y = 0; y<2; y+=0.1)
        {
            for(double z = 0; z<2; z+=0.1)
            {
                pt.x=x;
                pt.y=y;
                pt.z=z;
                cl.points.push_back(pt);
            }
        }
    }
    lslgeneric::writeToVRML("/home/tsv/ndt_tmp/final2.wrl",cl);
    LazzyGrid pr(0.33);
    NDTMap ndt2( &pr );
    ndt2.loadPointCloud( cl );
    ndt2.computeNDTCells();
    ndt2.writeToVRML("/home/tsv/ndt_tmp/example.wrl");

#endif


}

