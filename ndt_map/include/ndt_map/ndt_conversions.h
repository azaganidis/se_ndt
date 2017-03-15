#ifndef NDT_CONVERSIONS_HH
#define NDT_CONVERSIONS_HH

#include <ros/ros.h>
#include <vector>
#include <Eigen/Eigen>
#include <ndt_map/ndt_map.h>
#include <ndt_map/ndt_cell.h>
#include <ndt_map/NDTMapMsg.h>
#include <ndt_map/NDTCellMsg.h>
#include <nav_msgs/OccupancyGrid.h>
#include <string>


namespace lslgeneric{
  /** 
   *
   * \brief Message building fucntion
   * \details Converts an object of type NDTMap into NDTMapMsg
   * message. Message contain all the data strored in the object. 
   * @param[in] map Pointer to NDTMap object.
   * @param[out] msg formated message
   * @param[in] frame_name name of the coordination frame for the transformed map
   *   
   */
  bool toMessage(NDTMap *map, ndt_map::NDTMapMsg &msg,std::string frame_name){
    std::vector<lslgeneric::NDTCell*> map_vector=map->getAllInitializedCells();
    msg.header.stamp=ros::Time::now();
    msg.header.frame_id=frame_name;//is it in *map?    
    if(!map->getGridSizeInMeters(msg.x_size,msg.y_size,msg.z_size)){
	  ROS_ERROR("NO GRID SIZE");
	  return false;
    }
    if(!map->getCentroid(msg.x_cen,msg.y_cen,msg.z_cen)){
	  ROS_ERROR("NO GRID CENTER");
	  return false;
    }
    if(!map->getCellSizeInMeters(msg.x_cell_size,msg.y_cell_size,msg.z_cell_size)){
	  ROS_ERROR("NO CELL SIZE");
	  return false;
    }
    for (int cell_idx=0;cell_idx<map_vector.size();cell_idx++){
	  if(map_vector[cell_idx]->hasGaussian_){ //we only send a cell with gaussian
        ndt_map::NDTCellMsg cell;
        Eigen::Vector3d means=map_vector[cell_idx]->getMean();
        cell.mean_x=means(0);
        cell.mean_y=means(1);
        cell.mean_z=means(2);
        cell.occupancy=map_vector[cell_idx]->getOccupancyRescaled();
        Eigen::Matrix3d cov=map_vector[cell_idx]->getCov();
        for(int i=0;i<3;i++){
		  for(int j=0;j<3;j++){
            cell.cov_matrix.push_back(cov(i,j));
		  }
        }
        cell.N=map_vector[cell_idx]->getN();
        msg.cells.push_back(cell);
	  }
	  delete map_vector[cell_idx];
    }
    return true;
  }
  /** 
   *
   * \brief from message to NDTMap object
   * \details Converts ndt map message into a NDTMap object
   * @param[in,out] idx Pointer to lazy grid of the new NDTMap
   * @param[out] map Pointer to NDTMap object
   * @param[in] msg message to be converted
   * @param[out] frame_name name of the coordination frame of the map
   *   
   */
  bool fromMessage(LazyGrid* &idx, NDTMap* &map, ndt_map::NDTMapMsg msg, std::string &frame_name){
    if(!(msg.x_cell_size==msg.y_cell_size&&msg.y_cell_size==msg.z_cell_size)){ //we assume that voxels are cubes
	  ROS_ERROR("SOMETHING HAS GONE VERY WRONG YOUR VOXELL IS NOT A CUBE"); 
	  return false;
    }
    idx=new LazyGrid(msg.x_cell_size);
    map = new NDTMap(idx,msg.x_cen,msg.y_cen,msg.z_cen,msg.x_size,msg.y_size,msg.z_size);
    frame_name=msg.header.frame_id;
    int gaussians=0;
    for(int itr=0;itr<msg.cells.size();itr++){
	  Eigen::Vector3d mean;
	  Eigen::Matrix3d cov;
	  mean<<msg.cells[itr].mean_x,msg.cells[itr].mean_y,msg.cells[itr].mean_z;
	  int m_itr=0;
	  for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
		  cov(i,j)=msg.cells[itr].cov_matrix[m_itr];
		  m_itr++;
        }
	  }
	  map->addDistributionToCell(cov,mean,msg.cells[itr].N);
    }
    return true;
  }
  /**
   *
   * \brief builds ocuupancy grid message
   * \details Builds 2D occupancy grid map based on 2D NDTMap
   * @param[in] ndt_map 2D ndt map to conversion
   * @param[out] occ_grid 2D cost map
   * @param[in] resolution desired resolution of occupancy map
   * @param[in] name of cooridnation frame for the map (same as the NDT map has)
   * 
   */
  bool toOccupancyGrid(NDTMap *ndt_map, nav_msgs::OccupancyGrid &occ_grid, double resolution,std::string frame_id){//works only for 2D case
    double size_x, size_y, size_z;
    int size_x_cell_count, size_y_cell_count;
    double cen_x, cen_y, cen_z;
    double orig_x, orig_y;
    ndt_map->getGridSizeInMeters(size_x,size_y,size_z);
    ndt_map->getCentroid(cen_x,cen_y,cen_z);
    orig_x=cen_x-size_x/2.0;
    orig_y=cen_y-size_y/2.0;
    size_x_cell_count=int(size_x/resolution);
    size_y_cell_count=int(size_y/resolution);
    occ_grid.info.width=size_x_cell_count;
    occ_grid.info.height=size_y_cell_count;
    occ_grid.info.resolution=resolution;
    occ_grid.header.stamp=ros::Time::now();
    occ_grid.header.frame_id=frame_id;
    for(double px=orig_x+resolution/2.0;px<orig_x+size_x;px+=resolution){
      for(double py=orig_y+resolution/2.0;py<orig_y+size_x;py+=resolution){
        pcl::PointXYZ pt(px,py,0);
        lslgeneric::NDTCell *cell;
        if(!ndt_map->getCellAtPoint(pt, cell)){
          occ_grid.data.push_back(-1);
        }
        else if(cell == NULL){
          occ_grid.data.push_back(-1);
        }
        else{
          Eigen::Vector3d vec (pt.x,pt.y,pt.z);
          vec = vec-cell->getMean();                  
          double likelihood = vec.dot(cell-> getInverseCov()*vec);
          char s_likelihood;
          if(cell->getOccupancy()!=0.0){
            if(cell->getOccupancy()>0.0){
            if(std::isnan(likelihood)) s_likelihood = -1;
            likelihood = exp(-likelihood/2.0) + 0.1;
            likelihood = (0.5+0.5*likelihood);
            s_likelihood=char(likelihood*100.0);
            if(likelihood >1.0) s_likelihood =100;
            occ_grid.data.push_back(s_likelihood);
            }
            else{
              occ_grid.data.push_back(0);
            }
          }
          else{
             occ_grid.data.push_back(-1);
          }
        }
      }
    }    
    return true;
  } 
}
#endif
