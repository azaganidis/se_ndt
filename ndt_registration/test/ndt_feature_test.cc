// It works if the --p2p option is given!??
#include <iostream>
#include <boost/program_options.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/random.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>

#include <ndt_registration/ndt_matcher_d2d.h>
#include <ndt_registration/ndt_matcher_d2d_feature.h>
#include <ndt_map/ndt_map.h>
#include <ndt_map/cell_vector.h>
#include <ndt_map/lazy_grid.h>
#include <pointcloud_vrml/pointcloud_utils.h>

using namespace std;
using namespace lslgeneric;
namespace po = boost::program_options;

int main(int argc, char** argv)
{
    cout << "--------------------------------------------------" << endl;
    cout << "Small test program of CellVector + F2F matcher " << endl;
    cout << "--------------------------------------------------" << endl;

    po::options_description desc("Allowed options");
    Eigen::Matrix<double,6,1> pose_increment_v;
    string static_file_name, moving_file_name;
    int nb_clusters, nb_points;
    double std_dev, min, max;
    desc.add_options()
    ("help", "produce help message")
    ("x", po::value<double>(&pose_increment_v(0))->default_value(1.), "x pos gt offset")
    ("y", po::value<double>(&pose_increment_v(1))->default_value(1.), "y pos gt offset")
    ("z", po::value<double>(&pose_increment_v(2))->default_value(1.), "z pos gt offset")
    ("X", po::value<double>(&pose_increment_v(3))->default_value(0.1), "x axis rot gt offset")
    ("Y", po::value<double>(&pose_increment_v(4))->default_value(0.1), "y axis rot gt offset")
    ("Z", po::value<double>(&pose_increment_v(5))->default_value(0.1), "z axis rot gt offset")
    ("nb_clusters", po::value<int>(&nb_clusters)->default_value(20), "number of clusters")
    ("nb_points", po::value<int>(&nb_points)->default_value(10), "number of points per clusters")
    ("std_dev", po::value<double>(&std_dev)->default_value(0.1), "standard deviation of the points drawn from a normal distribution (here independent on the axes")
    ("min", po::value<double>(&min)->default_value(-10), "minimum center point")
    ("max", po::value<double>(&max)->default_value(10), "maximum center point")
    ("lazzy", "use lazzygrid")
    ("p2preg", "calculate NDTMatchF2F using two pointclouds")
    ("irregular_grid", "use irregular grid in the p2p registration")
    ("nfeaturecorr", "feature_f2f should NOT be evaluated (NDTMatcherFeatureF2F)")
    ("singleres", "use single resolution in the 'p2p' matching")
    ("nf2f", "if f2f should NOT be evaluated")
    ("usegt", "if the ground truth should be used as initial estimate")
    ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
        cout << desc << "\n";
        return 1;
    }
    bool use_lazzy = vm.count("lazzy");
    bool use_p2preg = vm.count("p2preg");
    bool use_irregular_grid = vm.count("irregular_grid");
    bool use_featurecorr = !vm.count("nfeaturecorr");
    bool use_singleres = vm.count("singleres");
    bool use_f2f = !vm.count("nf2f");
    bool usegt = vm.count("usegt");
    pcl::PointCloud<pcl::PointXYZ> static_pc, moving_pc, tmp_pc;
    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> gt_transform;

    // Generate the static point cloud + indices...
    boost::mt19937 rng;
    boost::uniform_real<> ud(min,max);
    boost::normal_distribution<> nd(0.0, std_dev);
    boost::variate_generator<boost::mt19937&, boost::normal_distribution<> > var_nor(rng, nd);
    boost::variate_generator<boost::mt19937&, boost::uniform_real<> > var_uni(rng, ud);
    std::vector<std::pair<int,int> > corresp;

    size_t index = 0;
    std::vector<std::vector<size_t> > all_indices;
    for (int i = 0; i < nb_clusters; i++)
    {
        std::vector<size_t> indices;
        double c_x = var_uni();
        double c_y = var_uni();
        double c_z = var_uni();
        for (int j = 0; j < nb_points; j++)
        {
            static_pc.push_back(pcl::PointXYZ(c_x+var_nor(), c_y+var_nor(), c_z+var_nor()));
            indices.push_back(index);
            index++;
        }
        all_indices.push_back(indices);
        corresp.push_back(std::pair<int,int>(i,nb_clusters-1-i)); // nb_clusters-1-i -> To check the da functions is working.
    }
    std::vector<std::vector<size_t> > all_indices_moving(all_indices.size());  // Reverse the correspondances (To check the da functions)
    std::reverse_copy(all_indices.begin(), all_indices.end(), all_indices_moving.begin());

    for (int i = 0; i < 1/*nb_clusters*/; i++)
    {
        tmp_pc.push_back(pcl::PointXYZ(i,i,i));
    }

    // Specify some offset...
    gt_transform = Eigen::Translation<double,3>(pose_increment_v(0),pose_increment_v(1),pose_increment_v(2))*
                   Eigen::AngleAxis<double>(pose_increment_v(3),Eigen::Vector3d::UnitX()) *
                   Eigen::AngleAxis<double>(pose_increment_v(4),Eigen::Vector3d::UnitY()) *
                   Eigen::AngleAxis<double>(pose_increment_v(5),Eigen::Vector3d::UnitZ()) ;


    std::vector<double> resolutions;
    if (use_singleres)
        resolutions.push_back(1.);
    NDTMatcherD2D<pcl::PointXYZ,pcl::PointXYZ> matcher(use_irregular_grid,!use_singleres,resolutions);
    Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> T_f2f,T_p2p, T_feat;
    T_f2f.setIdentity();
    T_p2p.setIdentity();
    T_feat.setIdentity();
    if (usegt)
    {
        T_f2f = gt_transform;
        T_p2p = gt_transform;
        T_feat = gt_transform;
    }

    moving_pc = transformPointCloud(gt_transform, static_pc);

    if (use_p2preg)
    {
        matcher.match(moving_pc, static_pc, T_p2p);
    }
    {
        double current_resolution = 1;

        SpatialIndex<pcl::PointXYZ>* index = NULL;
        if (use_lazzy)
        {
            index = new LazyGrid<pcl::PointXYZ>(current_resolution);
        }
        else
        {
            index = new CellVector<pcl::PointXYZ>();
        }

        NDTMap<pcl::PointXYZ> ndt(index);
        if (!use_lazzy)
            ndt.loadPointCloud( static_pc, all_indices );
        else
            ndt.loadPointCloud (static_pc );
        ndt.computeNDTCells();

        NDTMap<pcl::PointXYZ> mov(index);
        if (!use_lazzy)
            mov.loadPointCloud( moving_pc, all_indices_moving );
        else
            mov.loadPointCloud( moving_pc );
        mov.computeNDTCells();

        if (use_f2f)
        {
            matcher.match( mov, ndt, T_f2f );
        }
        if (use_featurecorr)
        {
            NDTMatcherFeatureD2D<pcl::PointXYZ,pcl::PointXYZ> matcher_feat(corresp);
            matcher_feat.match( mov, ndt, T_feat );
        }
        delete index;
    }

    std::cout<<"GT translation "<<gt_transform.translation().transpose()
             <<" (norm) "<<gt_transform.translation().norm()<<std::endl;
    std::cout<<"GT rotation "<<gt_transform.rotation().eulerAngles(0,1,2).transpose()
             <<" (norm) "<<gt_transform.rotation().eulerAngles(0,1,2).norm()<<std::endl;

    if (use_f2f)
    {
        std::cout<<"f2f translation "<<T_f2f.translation().transpose()
                 <<" (norm) "<<T_f2f.translation().norm()<<std::endl;
        std::cout<<"f2f rotation "<<T_f2f.rotation().eulerAngles(0,1,2).transpose()
                 <<" (norm) "<<T_f2f.rotation().eulerAngles(0,1,2).norm()<<std::endl;
    }

    if (use_featurecorr)
    {
        std::cout<<"feat translation "<<T_feat.translation().transpose()
                 <<" (norm) "<<T_feat.translation().norm()<<std::endl;
        std::cout<<"feat rotation "<<T_feat.rotation().eulerAngles(0,1,2).transpose()
                 <<" (norm) "<<T_feat.rotation().eulerAngles(0,1,2).norm()<<std::endl;
    }

    if (use_p2preg)
    {
        std::cout<<"p2preg translation "<<T_p2p.translation().transpose()
                 <<" (norm) "<<T_p2p.translation().norm()<<std::endl;
        std::cout<<"p2preg rotation "<<T_p2p.rotation().eulerAngles(0,1,2).transpose()
                 <<" (norm) "<<T_p2p.rotation().eulerAngles(0,1,2).norm()<<std::endl;
    }
    cout << "done." << endl;
}
