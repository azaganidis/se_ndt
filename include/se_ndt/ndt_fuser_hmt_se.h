#ifndef NDT_FUSER_HMT_SE_HH
#define NDT_FUSER_HMT_SE_HH
#define NO_NDT_VIZ 1
#ifndef NO_NDT_VIZ
#include <ndt_visualisation/ndt_viz.h>
#endif
#include <ndt_map/ndt_map.h>
#include <ndt_map/ndt_map_hmt.h>
#include <ndt_registration/ndt_matcher_d2d_2d.h>
#include "se_ndt/ndt_matcher_d2d_se.h"
#include <ndt_map/pointcloud_utils.h>
#include <se_ndt/se_ndt.hpp>

#include <Eigen/Eigen>
#include <pcl/point_cloud.h>
#include <sys/time.h>

//#define BASELINE

namespace perception_oru {
/**
  * \brief This class fuses new point clouds into a common ndt map reference, keeping tack of the 
  * camera postion.
  * \author Jari, Todor
  */
class NDTFuserHMT_SE : NDTMatch_SE{
    public:
		unsigned int NumInputs;
		Eigen::Affine3d Tnow, Tlast_fuse, Todom; ///< current pose
		perception_oru::NDTMap ***map;		 ///< da map
		perception_oru::NDTMap ***mapLocal;		 ///< da map
		bool checkConsistency;			 ///perform a check for consistency against initial estimate
		vector<float> resolutions;
		initializer_list<float> ignore,size;
		initializer_list<int> resolutions_order,tails;
		double max_translation_norm, max_rotation_norm;
		double sensor_range;
		float removeProbability;
		bool fuseIncomplete, beHMT,canUpdate;
		std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud_c;
		int ctr;
		std::string prefix;
		std::string hmt_map_dir;
#ifndef NO_NDT_VIZ
	NDTViz *viewer;
#endif
	FILE *fAddTimes, *fRegTimes;

	NDTFuserHMT_SE(Eigen::Affine3d a,initializer_list<float> b,initializer_list<int> c,initializer_list<float> d,initializer_list<int> e,initializer_list<float> ig,float removeP,int max_iter);
	    
	~NDTFuserHMT_SE()
	{
		for(unsigned int i=0;i<resolutions.size();i++)
		{
			for(auto j=0;j<NumInputs;j++)
			{
				delete map[i][j];
				delete mapLocal[i][j];
			}
			delete[] map[i];
			delete[] mapLocal[i];
		}
	    delete[] map;
#ifndef NO_NDT_VIZ
	    delete viewer;
#endif
	    if(fAddTimes!=NULL) fclose(fAddTimes);
	    if(fRegTimes!=NULL) fclose(fRegTimes);
	}

	double getDoubleTime()
	{
	    struct timeval time;
	    gettimeofday(&time,NULL);
	    return time.tv_sec + time.tv_usec * 1e-6;
	}
	void setSensorPose(Eigen::Affine3d spose){
	    sensor_pose = spose;
	}
	
	bool saveMap() {
	    if(map == NULL) return false;
		bool result=true;
	    if(beHMT) 
		{
			for(unsigned int j=0;j<resolutions.size();j++)
			{
			for(unsigned int i=0;i<NumInputs;i++)
			{
				perception_oru::NDTMapHMT *map_hmt = dynamic_cast<perception_oru::NDTMapHMT*> (map[j][i]);
				if(map_hmt==NULL) return false;
				result&=(map_hmt->writeTo()==0);
			}
			}
			return result;
	    }
		else 
		{
			for(unsigned int j=0;j<resolutions.size();j++)
			{
				for(unsigned int i=0;i<NumInputs;i++)
				{
					char fname[1000];
					snprintf(fname,999,"%s/%s_map%d_%d.jff",hmt_map_dir.c_str(),prefix.c_str(),i,j);
					result&= (map[j][i]->writeToJFF(fname)==0);
				}
			}
			return result;
	    }
	}

	/**
	 *
	 *
	 */
	Eigen::Affine3d update(Eigen::Affine3d Tmotion, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,initializer_list<vector<double> > attributes);
	Eigen::Affine3d match(Eigen::Affine3d Tmotion, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud,initializer_list<vector<double> > attributes);
	bool updateMap();
	Eigen::Matrix<double,6,6> getPoseCovariance(Eigen::Affine3d T);
	perception_oru::NDTMatcherD2D_SE matcher;
    private:

	double translation_fuse_delta, rotation_fuse_delta;
	bool visualize,firstRun;

	Eigen::Affine3d sensor_pose;
	perception_oru::NDTMatcherD2D_2D matcher2D;
	Eigen::Vector3d localMapSize;

    public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW

};
}
#endif
