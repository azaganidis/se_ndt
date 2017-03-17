namespace lslgeneric
{

template< typename PointT>
pcl::PointCloud<PointT> transformPointCloud(Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> &Tr, const pcl::PointCloud<PointT> &pc)
{
    Eigen::Transform<float,3,Eigen::Affine,Eigen::ColMajor> T = Tr.cast<float>();
    pcl::PointCloud<PointT> cloud;
    for(unsigned int pit=0; pit<pc.points.size(); ++pit)
    {
        PointT thisPoint = pc.points[pit];
        Eigen::Map<Eigen::Vector3f> pt((float*)&thisPoint,3);
        pt = T*pt;
        cloud.points.push_back(thisPoint);
    }
    cloud.width = pc.width;
    cloud.height = pc.height;
    return cloud;
}

template< typename PointT>
void transformPointCloudInPlace(Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> &Tr, pcl::PointCloud<PointT> &pc)
{
    Eigen::Transform<float,3,Eigen::Affine,Eigen::ColMajor> T = Tr.cast<float>();
    for(unsigned int pit=0; pit<pc.points.size(); ++pit)
    {
        Eigen::Map<Eigen::Vector3f> pt((float*)&pc.points[pit],3);
		//std::cout<<-atan2(pc.points[pit].y,pc.points[pit].x)<<std::endl;
        pt = T*pt;
    }
}
template< typename PointT>
void transformPointCloudsInPlace(Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> &Tr, pcl::PointCloud<PointT> *pc,int NumInputs)
{
    Eigen::Transform<float,3,Eigen::Affine,Eigen::ColMajor> T = Tr.cast<float>();
	for(unsigned int j=0;j<NumInputs;j++)
		for(unsigned int pit=0; pit<pc[j].points.size(); ++pit)
		{
			Eigen::Map<Eigen::Vector3f> pt((float*)&pc[j].points[pit],3);
			//std::cout<<-atan2(pc.points[pit].y,pc.points[pit].x)<<std::endl;
			pt = T*pt;
		}
}
template< typename PointT>
double geomDist(PointT p1, PointT p2)
{
    Eigen::Vector3d v;
    v << p1.x-p2.x, p1.y-p2.y, p1.z-p2.z;
    return v.norm();
}

}
