#include "pose_optimizer.h"
#include <iostream>
#include <fstream>
#include <Eigen/Core>
using namespace std;
int main(int argc, char** argv)
{
    PoseOptimizer pose_graph;
    ifstream pose_in("LC1.txt");
    Eigen::Affine3d T;
    T.setIdentity();
    pose_graph.addPose(T,0);
    int lpcount=0;
    int last_loop=0;
    while(pose_in.good()){
        int fId=0, tId=0;
        double pose[7]{};
        pose_in>>fId>>tId;
        if(fId==tId)
            break;
        if(fId>atoi(argv[4]))
            break;
        if(tId-fId==1)
        {
            for(int i=0;i<7;i++)
                pose_in>>pose[i];
            Eigen::Map<Eigen::Matrix<double,3,1> > poseT(pose);
            Eigen::Map<Eigen::Quaterniond > poseQ(pose+3);
            Eigen::Affine3d Td;
            Td.translation()= poseT;
            Td.linear()=poseQ.toRotationMatrix();
            T=T*Td;
            pose_graph.addPose(T,tId);
            pose_graph.addConstraint(Td,fId, tId, atoi(argv[3])*Eigen::Matrix<double, 7,7>::Identity());
        }
        else{
            float hist_score, match_score;
            pose_in>>hist_score>>match_score;
            for(int i=0;i<7;i++)
                pose_in>>pose[i];
            Eigen::Map<Eigen::Matrix<double,3,1> > poseT(pose);
            Eigen::Map<Eigen::Quaterniond > poseQ(pose+3);
            Eigen::Affine3d Td;
            Td.translation()= poseT;
            Td.linear()=poseQ.toRotationMatrix();

            Eigen::Affine3d Tf=pose_graph.getT(fId);
            Eigen::Affine3d Tt=pose_graph.getT(tId);
            Eigen::Affine3d Tc;
            Tc.setIdentity();
            Tc= Tf.inverse()*Td*Tt;//not terrible
            //Tc.translation()=Td.translation()-Tf.translation()+Tt.translation();
            Tc=Td;


            if(hist_score<atof(argv[1]) && match_score < atof(argv[2]) && tId>last_loop+atoi(argv[5]))
            {
                pose_graph.addConstraint(Tc,fId, tId, Eigen::Matrix<double, 7,7>::Identity());
                last_loop=tId;
                lpcount++;
            }
        }
    }
    std::cout<<"N loop closures: "<<lpcount<<std::endl;
    return 0;
}

