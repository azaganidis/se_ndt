
#ifndef OTHER_HPP
#define OTHER_HPP
std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> getSegmentsFast(pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudIn, int NumInputs);
perception_oru::NDTMap*** allocateMap(std::vector<float> &resolutions, std::vector<float> &size,int nIn);
void destroyMap(perception_oru::NDTMap ***map,unsigned int nRes,unsigned int nIn);
void loadMap(perception_oru::NDTMap **map,
        std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> input_clouds,
        float sensor_range=100);
Eigen::Matrix<double,6,6> getHes(Eigen::Matrix<double,6,6> Hessian,Eigen::Matrix<double,6,1> score_gradient);
Eigen::Matrix<double,7,6> getJacobian(Eigen::VectorXd v);
void print_vals();
std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> getSegmentsFast(
        pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudIn, int NumInputs);
#endif
