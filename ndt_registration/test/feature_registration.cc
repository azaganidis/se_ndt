#include <limits>
#include <fstream>
#include <vector>
#include <Eigen/Core>
#include "pcl/point_types.h"
#include "pcl/point_cloud.h"
#include "pcl/io/pcd_io.h"
#include "pcl/kdtree/kdtree_flann.h"
#include "pcl/filters/passthrough.h"
#include "pcl/filters/voxel_grid.h"
#include "pcl/features/normal_3d.h"
#include "pcl/features/fpfh.h"
#include "pcl/registration/ia_ransac.h"
#include "pcl/registration/icp.h"
#include <pcl/surface/mls.h>

#include <NDTMatcherF2F.hh>
#include <PointCloudUtils.hh>

class FeatureCloud
{
public:
    // A bit of shorthand
    typedef pcl::PointCloud<pcl::PointXYZ> PointCloud;
    typedef pcl::PointCloud<pcl::Normal> SurfaceNormals;
    typedef pcl::PointCloud<pcl::FPFHSignature33> LocalFeatures;
    //typedef pcl::KdTreeFLANN<pcl::PointXYZ> SearchMethod;


    FeatureCloud () :
        //search_method_xyz_ (new SearchMethod),
        normal_radius_ (0.05),
        feature_radius_ (0.05)
    {
    }

    ~FeatureCloud () {}

    // Process the given cloud
    void
    setInputCloud (PointCloud::Ptr xyz)
    {
        xyz_ = xyz;
        processInput ();
    }

    // Load and process the cloud in the given PCD file
    void
    loadInputCloud (const std::string &pcd_file)
    {
        xyz_ = PointCloud::Ptr (new PointCloud);
        pcl::io::loadPCDFile (pcd_file, *xyz_);
        processInput ();
    }

    // Get a pointer to the cloud 3D points
    PointCloud::Ptr
    getPointCloud () const
    {
        return (xyz_);
    }

    // Get a pointer to the cloud of 3D surface normals
    SurfaceNormals::Ptr
    getSurfaceNormals () const
    {
        return (normals_);
    }

    // Get a pointer to the cloud of feature descriptors
    LocalFeatures::Ptr
    getLocalFeatures () const
    {
        return (features_);
    }

protected:
    // Compute the surface normals and local features
    void
    processInput ()
    {
        computeSurfaceNormals ();
        computeLocalFeatures ();
    }

    // Compute the surface normals
    void
    computeSurfaceNormals ()
    {
        normals_ = SurfaceNormals::Ptr (new SurfaceNormals);

        pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> norm_est;
        norm_est.setInputCloud (xyz_);
        norm_est.setSearchMethod (search_method_xyz_);
        norm_est.setRadiusSearch (normal_radius_);
        norm_est.compute (*normals_);
    }

    // Compute the local feature descriptors
    void
    computeLocalFeatures ()
    {
        features_ = LocalFeatures::Ptr (new LocalFeatures);

        pcl::FPFHEstimation<pcl::PointXYZ, pcl::Normal, pcl::FPFHSignature33> fpfh_est;
        fpfh_est.setInputCloud (xyz_);
        fpfh_est.setInputNormals (normals_);
        fpfh_est.setSearchMethod (search_method_xyz_);
        fpfh_est.setRadiusSearch (feature_radius_);
        fpfh_est.compute (*features_);
    }

private:
    // Point cloud data
    PointCloud::Ptr xyz_;
    SurfaceNormals::Ptr normals_;
    LocalFeatures::Ptr features_;
    //SearchMethod::Ptr search_method_xyz_;
    pcl::Feature<pcl::PointXYZ, pcl::Normal>::KdTreePtr search_method_xyz_;

    // Parameters
    float normal_radius_;
    float feature_radius_;
};

class TemplateRegistration
{
public:

    // A struct for storing alignment results

    typedef Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> Result;

    TemplateRegistration () :
        min_sample_distance_ (0.05),
        max_correspondence_distance_ (0.01*0.01),
        nr_iterations_ (500)
    {
        // Intialize the parameters in the Sample Consensus Intial Alignment (SAC-IA) algorithm
        sac_ia_.setMinSampleDistance (min_sample_distance_);
        sac_ia_.setMaxCorrespondenceDistance (max_correspondence_distance_);
        sac_ia_.setMaximumIterations (nr_iterations_);
        sac_ia_.setNumberOfSamples (10);
    }

    ~TemplateRegistration () {}

    // Set the given cloud as the target to which the templates will be aligned
    void
    setTargetCloud (FeatureCloud &target_cloud)
    {
        target_ = target_cloud;
        sac_ia_.setInputTarget (target_cloud.getPointCloud ());
        sac_ia_.setTargetFeatures (target_cloud.getLocalFeatures ());
    }

    // Align the moving cloud to the target specified by setTargetCloud ()
    void
    align (FeatureCloud &moving_cloud, TemplateRegistration::Result &result)
    {
        sac_ia_.setInputCloud (moving_cloud.getPointCloud ());
        sac_ia_.setSourceFeatures (moving_cloud.getLocalFeatures ());

        pcl::PointCloud<pcl::PointXYZ> registration_output;
        sac_ia_.align (registration_output);

        result = sac_ia_.getFinalTransformation ().cast<double>();
    }

private:
    // A list of template clouds and the target to which they will be aligned
    FeatureCloud target_;

    // The Sample Consensus Initial Alignment (SAC-IA) registration routine and its parameters
    pcl::SampleConsensusInitialAlignment<pcl::PointXYZ, pcl::PointXYZ, pcl::FPFHSignature33> sac_ia_;
    float min_sample_distance_;
    float max_correspondence_distance_;
    float nr_iterations_;
};

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
    //*cloud_in = moving;
    //*cloud_out= fixed;


    pcl::IterativeClosestPoint<pcl::PointXYZ, pcl::PointXYZ> icp;

    icp.setMaximumIterations(10000);
    std::cout<<"max itr are "<<icp.getMaximumIterations()<<std::endl;
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

// Align two point clouds based on the features
int
main (int argc, char **argv)
{
    if (argc < 3)
    {
        printf ("No targets given!\n");
        return (-1);
    }

    struct timeval tv_start, tv_end1, tv_end2;
    TemplateRegistration::Result ToutFPFH, ToutICP, ToutNDT, Tout;

    // Load the target cloud PCD file
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloudM (new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr cloudF (new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ> Final, Final0, Final1;

    *cloudM = lslgeneric::readVRML(argv[1]);
    *cloudF = lslgeneric::readVRML(argv[2]);

    double __res[] = {0.2, 0.4, 1, 2};
    std::vector<double> resolutions (__res, __res+sizeof(__res)/sizeof(double));
    lslgeneric::NDTMatcherF2F matcherF2F(false, false, false, resolutions);
    bool ret = matcherF2F.match(*cloudF,*cloudM,ToutNDT);
    Final1 = lslgeneric::transformPointCloud(ToutNDT,*cloudM);

    /* char f[50];
     snprintf(f,49,"/home/tsv/ndt_tmp/c2_offset.wrl");
     FILE *fo = fopen(f,"w");
     fprintf(fo,"#VRML V2.0 utf8\n");
     lslgeneric::writeToVRML(fo,*cloudM,Eigen::Vector3d(1,0,0));
     lslgeneric::writeToVRML(fo,Final1,Eigen::Vector3d(0,1,1));
     lslgeneric::writeToVRML(fo,*cloudF,Eigen::Vector3d(1,1,1));
     fclose(fo);
     return (0);
     */
    //start timing
    gettimeofday(&tv_start,NULL);

    pcl::VoxelGrid<pcl::PointXYZ> gr1,gr2;
    gr1.setLeafSize(0.05,0.05,0.05);
    gr2.setLeafSize(0.05,0.05,0.05);

    gr1.setInputCloud(cloudM);
    gr2.setInputCloud(cloudF);

    gr1.filter(*cloudM);
    gr2.filter(*cloudF);

    cloudM->height = 1;
    cloudM->width = cloudM->points.size();
    cloudF->height = 1;
    cloudF->width = cloudF->points.size();
    cloudM->is_dense = false;
    cloudF->is_dense = false;

    // Assign to the target FeatureCloud
    FeatureCloud target_cloud, moving_cloud;
    target_cloud.setInputCloud (cloudF);
    moving_cloud.setInputCloud (cloudM);

    TemplateRegistration templateReg;
    templateReg.setTargetCloud (target_cloud);

    // Find the best template alignment
    templateReg.align(moving_cloud,ToutFPFH);
    //stop timing1
    gettimeofday(&tv_end1,NULL);

    std::cout<<"ToutFPFH: "<<ToutFPFH.translation().transpose()<<"\n"<<ToutFPFH.rotation()<<std::endl;
    Final0 = lslgeneric::transformPointCloud(ToutFPFH,*cloudM);

    bool converged = matchICP(*cloudF,Final0,ToutICP);
    Final = lslgeneric::transformPointCloud(ToutICP,Final0);
    //stop timing2
    gettimeofday(&tv_end2,NULL);


    Tout = ToutFPFH*(ToutNDT.inverse());
    std::cout<<"FPFH\n";
    std::cout<<"E translation "<<Tout.translation().transpose()
             <<" (norm) "<<Tout.translation().norm()<<std::endl;
    std::cout<<"E rotation "<<Tout.rotation().eulerAngles(0,1,2).transpose()
             <<" (norm) "<<Tout.rotation().eulerAngles(0,1,2).norm()<<std::endl;
    std::cout<<" TIME: "<<
             (tv_end1.tv_sec-tv_start.tv_sec)*1000.+(tv_end1.tv_usec-tv_start.tv_usec)/1000.<<std::endl;

    Tout = ToutICP*ToutFPFH*(ToutNDT.inverse());
    std::cout<<"FPFH+ICP\n";
    std::cout<<"E translation "<<Tout.translation().transpose()
             <<" (norm) "<<Tout.translation().norm()<<std::endl;
    std::cout<<"E rotation "<<Tout.rotation().eulerAngles(0,1,2).transpose()
             <<" (norm) "<<Tout.rotation().eulerAngles(0,1,2).norm()<<std::endl;
    std::cout<<" TIME: "<<
             (tv_end2.tv_sec-tv_start.tv_sec)*1000.+(tv_end2.tv_usec-tv_start.tv_usec)/1000.<<std::endl;

    char fname[50];
    snprintf(fname,49,"/home/tsv/ndt_tmp/c2_offset.wrl");
    FILE *fout = fopen(fname,"w");
    fprintf(fout,"#VRML V2.0 utf8\n");
    lslgeneric::writeToVRML(fout,*cloudM,Eigen::Vector3d(1,0,0));
    lslgeneric::writeToVRML(fout,Final0,Eigen::Vector3d(0,0,1));
    lslgeneric::writeToVRML(fout,Final1,Eigen::Vector3d(0,1,1));
    lslgeneric::writeToVRML(fout,Final,Eigen::Vector3d(0,1,0));
    lslgeneric::writeToVRML(fout,*cloudF,Eigen::Vector3d(1,1,1));
    fclose(fout);



    return (0);
}
