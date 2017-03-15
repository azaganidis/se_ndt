#include<NDTMapBuilder.hh>
#include <OctTree.hh>
#include <NDTMatcherF2F.hh>
#include <PointCloudUtils.hh>

using namespace std;
using namespace lslgeneric;

int main(int argc, char**argv)
{

    if(argc != 7)
    {
        cout<<"Usage: "<<argv[0]<<" N_CLOUDS clouds_prefix (bKinect|bLaser) (doHistogram=0|1) (bD2D|bP2D|bICP) offset\n \
	    \t bKinect -> use settings more suitable for a kinect\n \
	    \t bLaser -> use settings more suitable for a kinect\n \
	    \t bD2D . p2D ICP -> choice of registration algorithm\n";
        return -1;
    }

    bool bKinect = (strncmp(argv[3],"bKinect",7) == 0);
    bool doHistogram = (atoi(argv[4]) == 1);
    bool bICP = (strncmp(argv[5],"bICP",4) == 0);
    bool bP2D = (strncmp(argv[5],"bP2D",4) == 0);
    bool bD2D = ((strncmp(argv[5],"bD2D",4) == 0) || !(bICP || bP2D));
    int offset = atoi(argv[6]);

    lslgeneric::OctTree histogramPrototype;
    std::vector<double> resolutions;

    double __res1[] = {0.5, 1, 2, 4};
    double __res2[] = {0.1,0.2,1,2,4};
    if(!bKinect)
    {
        resolutions = std::vector<double>(__res1, __res1+sizeof(__res1)/sizeof(double));
        lslgeneric::OctTree::BIG_CELL_SIZE = 2;
        lslgeneric::OctTree::SMALL_CELL_SIZE = 0.2;

    }
    else
    {
        resolutions = std::vector<double>(__res2, __res2+sizeof(__res2)/sizeof(double));
        lslgeneric::OctTree::BIG_CELL_SIZE = 0.5;
        lslgeneric::OctTree::SMALL_CELL_SIZE = 0.1;

    }
    lslgeneric::OctTree::parametersSet = true;
    lslgeneric::NDTMatcherF2F matcherF2F(false, false, false, resolutions);
    lslgeneric::NDTMatcher matcherP2F(resolutions);

    NDTMapBuilder mapper(doHistogram);
    if(bD2D)
    {
        mapper.setMatcherF2F(&matcherF2F);
        cout<<"setting to D2D matcher\n";
    }
    else if(bP2D)
    {
        mapper.setMatcherP2F(&matcherP2F);
        cout<<"setting to P2D matcher\n";
    }
    else
    {
        mapper.setICP();
        cout<<"setting to ICP matcher\n";
    }

    mapper.tr = histogramPrototype;

    int N_CLOUDS = atoi(argv[1]);
    char fname[600];

    double MAX_DIST = 26;

    if(bKinect) MAX_DIST = 5;

    for (int i=offset; i<N_CLOUDS; i++)
    {
        snprintf(fname,600,"%s%03d.wrl",argv[2],i);
        cout<<fname<<endl;
        pcl::PointCloud<pcl::PointXYZ> cl = lslgeneric::readVRML(fname);
        pcl::PointCloud<pcl::PointXYZ> filtered;
        for(int q=0; q<cl.points.size(); q++)
        {
            double dist = sqrt(pow(cl.points[q].x,2)+pow(cl.points[q].y,2)+pow(cl.points[q].z,2));
            if(dist<MAX_DIST)
            {
                if(!bKinect)
                {
                    filtered.points.push_back(cl.points[q]);
                }
                else
                {
                    pcl::PointXYZ pNew;
                    pNew.x = cl.points[q].z;
                    pNew.y = cl.points[q].x;
                    pNew.z = -cl.points[q].y;
                    filtered.points.push_back(pNew);
                }

            }
        }

        cout<<"adding cloud number "<<i<<endl;
        mapper.addScan(filtered);
    }

    snprintf(fname,600,"%s.g2o",argv[2]);
    mapper.saveG2OlogFile(fname);

    snprintf(fname,600,"%s.dat",argv[2]);
    mapper.saveDatlogFile(fname);

    snprintf(fname,600,"%s_COMPLETE.wrl",argv[2]);
    mapper.theMotherOfAllPointClouds(fname);

    return 0;
}
