#ifndef RVIZ_NDT_H
#define RVIZ_NDT_H
#include <ndt_map/ndt_map.h>
#include <ros/ros.h>
#include <visualization_msgs/Marker.h>
std::array<float,3> value_to_rgb(float i, float max)
{
    float H=i/max*360;
    float C=1;
    float X=1-abs( int(H/60)%2 - 1);
    std::array<float,3> r;
    if(H<60){ r[0]=C; r[1]=X; r[2]=0;}
    else if(H<120){ r[0]=X; r[1]=C; r[2]=0;}
    else if(H<180){ r[0]=0; r[1]=C; r[2]=X;}
    else if(H<240){ r[0]=0; r[1]=X; r[2]=C;}
    else if(H<300){ r[0]=X; r[1]=0; r[2]=C;}
    else if(H<360){ r[0]=C; r[1]=0; r[2]=X;}
    return r;
}
class ndt_rviz{
    std::vector<ros::Publisher> marker_pub;
    int NumSem;
    int *max_num=NULL;
    public:
    ros::Duration dur;
    ndt_rviz(ros::NodeHandle &n, int n_res, std::string prefix = "vm_res"):dur()
    {
        for(int i=0;i<n_res;i++)
            marker_pub.push_back(n.advertise<visualization_msgs::Marker>(prefix+std::to_string(i),10000));
    }
    void show_cell(const Eigen::Matrix3d &m_cov, const Eigen::Vector3d &m_mean,float occupancy,ros::Time& cl_time,int iRes, int iSem, int ID){
        const double d=m_cov.determinant();
        if (d==0 || d!=d) // Note: "d!=d" is a great test for invalid numbers, don't remove!
            return;
        Eigen::EigenSolver<Eigen::Matrix3d> es(m_cov);
        Eigen::Matrix3d m_eigVec;
        Eigen::Vector3d m_eigVal; 
        m_eigVal = es.pseudoEigenvalueMatrix().cwiseSqrt().diagonal();
        m_eigVec = es.pseudoEigenvectors();
        //if(!(m_eigVal(0,0) != 0.0 && m_eigVal(1,1) != 0.0 && m_eigVal(2,2) != 0.0))
        //    return;
        std::vector<int> idx(3);
        iota(idx.begin(),idx.end(),0);
        std::sort(idx.begin(), idx.end(), [&m_eigVal](int i1,int i2){return m_eigVal(i1)<m_eigVal(i2);});
        Eigen::Vector3d v1, v2=m_eigVec.row(idx[0]);
        v1<<1,0,0;
        Eigen::Quaterniond qR = Eigen::Quaterniond::FromTwoVectors(v1,v2);
        visualization_msgs::Marker marker;
        marker.header.frame_id = "/world";
        marker.header.stamp = cl_time;
        marker.ns = "sem"+std::to_string(iSem);
        marker.id=ID;
        marker.type=visualization_msgs::Marker::SPHERE;
        if(max_num[iSem]>ID)
            marker.action = visualization_msgs::Marker::MODIFY;
        else
        {
            max_num[iSem]++;
            marker.action = visualization_msgs::Marker::ADD;
        }
        marker.pose.position.x = m_mean(0);
        marker.pose.position.y = m_mean(1);
        marker.pose.position.z = m_mean(2);
        marker.pose.orientation.x = qR.x();
        marker.pose.orientation.y = qR.y();
        marker.pose.orientation.z = qR.z();
        marker.pose.orientation.w = qR.w();
        marker.scale.x = m_eigVal(idx[0]);
        marker.scale.y = m_eigVal(idx[1])>0.1?m_eigVal(idx[1]):(m_eigVal(idx[0])/10);
        marker.scale.z = m_eigVal(idx[2])>0.1?m_eigVal(idx[2]):(m_eigVal(idx[0])/10);
        marker.scale.x*=3;
        marker.scale.y*=3;
        marker.scale.z*=3;
        std::array<float,3> rgb_val=value_to_rgb(iSem,NumSem );
        marker.color.r = rgb_val[0];
        marker.color.g = rgb_val[1];
        marker.color.b = rgb_val[2];
        //std::cerr<<rgb_val[0]<<" "<<rgb_val[1]<<" "<<rgb_val[1]<<"\n";
        marker.color.a = sqrt(occupancy/150);
        //marker.color.a = 1;
        //if(marker.color.a<0.5)
        //    marker.color.a=0.5;
        marker.lifetime = dur;
        marker_pub[iRes].publish(marker);
    }
    void clearMarkers(int n_res){
        visualization_msgs::Marker marker;
        marker.header.frame_id = "/world";
        marker.header.stamp = ros::Time();
        marker.action = visualization_msgs::Marker::DELETEALL;
        for(int i=0;i<n_res;i++)
            marker_pub[i].publish(marker);
    }
	void plotNDTs(perception_oru::NDTMap ***maps,int n_res, int n_sem, ros::Time cl_time){
        clearMarkers(n_res);
        NumSem=n_sem;
        for(int iRes=0;iRes<n_res;iRes++)
            for(int iSem=0;iSem<n_sem;iSem++)
            {
				std::vector<perception_oru::NDTCell*> tempMap=(maps[iRes][iSem]->getAllCells());
				for(unsigned int i=0;i<tempMap.size();i++)
                {
					Eigen::Vector3d m = tempMap[i]->getMean();
					float occupancy = tempMap[i]->getOccupancy();
                    Eigen::Matrix3d cov = tempMap[i]->getCov();
                    show_cell(cov, m, occupancy,cl_time,iRes, iSem, i);
                    delete tempMap[i];
				}
            }
    }
	void plotNDTs(std::vector<std::vector<perception_oru::NDTCell*> > &tempMap)
    {
        //clearMarkers(0);
        ros::Time tn=ros::Time::now();
        NumSem=tempMap.size();
        if(!max_num)
            max_num= (int*)calloc(NumSem,sizeof(int));
        int dbg_max=0;
        for(unsigned int j=0;j<tempMap.size();j++)
        {
            if(dbg_max<max_num[j])dbg_max=max_num[j];
            for(unsigned int i=0;i<tempMap[j].size();i++)
            {
                Eigen::Vector3d m = tempMap[j][i]->getMean();
                float occupancy = tempMap[j][i]->getOccupancy();
                if(occupancy<32)
                    continue;
                Eigen::Matrix3d cov = tempMap[j][i]->getCov();
                show_cell(cov, m, occupancy,tn,0, j, i);
                delete tempMap[j][i];
            }
            tempMap[j].clear();
        }
        std::cout<<"****************\n"<<dbg_max<<std::endl;
    }
};
#endif
