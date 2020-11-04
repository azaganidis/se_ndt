#include <ndt_map/ndt_map.h>
#include <pcl/point_cloud.h>

void loadMap(perception_oru::NDTMap **map,std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> input_clouds,float sensor_range)
{
    #pragma omp parallel num_threads(N_THREADS)
    {
        #pragma omp for
        for(size_t i=0;i<input_clouds.size();i++)
        {
            map[i]->loadPointCloud(*input_clouds[i],sensor_range);
            map[i]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE);
        }
    }
}
