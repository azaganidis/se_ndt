#include <vector>
#include <ndt_map/ndt_map.h>
#include <pcl/point_cloud.h>
#include <other.hpp>
using namespace perception_oru;
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
Eigen::Matrix<double,6,6> getHes(Eigen::Matrix<double,6,6> Hessian,Eigen::Matrix<double,6,1> score_gradient)
{
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double,6,6> > Sol (Hessian);
        Eigen::Matrix<double,6,1> evals = Sol.eigenvalues().real();
        double minCoeff = evals.minCoeff();
        double maxCoeff = evals.maxCoeff();
        if(minCoeff < 0)  //|| evals.minCoeff()) // < 10e-5*evals.maxCoeff()) 
		{
			Eigen::Matrix<double,6,6> evecs = Sol.eigenvectors().real();
			double regularizer = score_gradient.norm();
			regularizer = regularizer + minCoeff > 0 ? regularizer : 0.001*maxCoeff - minCoeff;
			//double regularizer = 0.001*maxCoeff - minCoeff;
			Eigen::Matrix<double,6,1> reg;
			//ugly
			reg<<regularizer,regularizer,regularizer,regularizer,regularizer,regularizer;
			evals += reg;
			Eigen::Matrix<double,6,6> Lam;
			Lam = evals.asDiagonal();
			Hessian = evecs*Lam*(evecs.transpose());
		}
		return Hessian;
}

Eigen::Matrix<double,7,6> getJacobian(Eigen::VectorXd v)
{
	Eigen::Matrix<double,7,6> m;
	double x=v(0)/2;
	double y=v(1)/2;
	double z=v(2)/2;
	double cx=cos(x);
	double cy=cos(y);
	double cz=cos(z);
	double sx=sin(x);
	double sy=sin(y);
	double sz=sin(z);
	m(0,3)=0;
	m(1,3)=0;
	m(2,3)=0;
	m(3,3)=-sx*sy*cz+cx*cy*sz;
	m(4,3)=-sx*cy*sz-cx*sy*cz;
	m(5,3)=+cx*cy*cz+sx*sy*sz;
	m(6,3)=-sx*cy*cz+cx*sy*sz;

	m(0,4)=0;
	m(1,4)=0;
	m(2,4)=0;
	m(3,4)=+cx*cy*cz-sx*sy*sz;
	m(4,4)=-cx*sy*sz-sx*cy*cz;
	m(5,4)=-sx*sy*cz-cx*cy*sz;
	m(6,4)=-cx*sy*cz+sx*cy*sz;

	m(0,5)=0;
	m(1,5)=0;
	m(2,5)=0;
	m(3,5)=-cx*sy*sz+sx*cy*cz;
	m(4,5)=+cx*cy*cz+sx*sy*sz;
	m(5,5)=-sx*cy*sz-cx*sy*cz;
	m(6,5)=-cx*cy*sz+sx*sy*cz;

	m=m/2;

	m(0,0)=1;//1
	m(1,0)=0;
	m(2,0)=0;
	m(3,0)=0;
	m(4,0)=0;
	m(5,0)=0;
	m(6,0)=0;

	m(0,1)=0;
	m(1,1)=1;//1	
	m(2,1)=0;
	m(3,1)=0;
	m(4,1)=0;
	m(5,1)=0;
	m(6,1)=0;

	m(0,2)=0;
	m(1,2)=0;
	m(2,2)=1;//1
	m(3,2)=0;
	m(4,2)=0;
	m(5,2)=0;
	m(6,2)=0;

	return m;
}

std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> getSegmentsFast(pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudIn, int NumInputs)
{
	int cloudSize = laserCloudIn->points.size();
    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud;
	for(int i=0;i<NumInputs;i++)
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT(new pcl::PointCloud<pcl::PointXYZ>);
		laserCloud.push_back(cloudT);
	}
	for(int i=0;i<cloudSize;i++)
	{
            auto pnt=laserCloudIn->points[i];
            int atr = round(pnt.intensity);
            if(NumInputs==1)
                atr=0;
            pcl::PointXYZ point;
            point.x=pnt.x;
            point.y=pnt.y;
            point.z=pnt.z;
            laserCloud[atr]->points.push_back(point);
    }
	return laserCloud;
}
NDTMap*** allocateMap(std::vector<float> &resolutions, std::vector<float> &size,int nIn,bool whatisthat)
{
    NDTMap ***map=new NDTMap **[resolutions.size()];
    for(unsigned int i=0;i<resolutions.size();i++)
    {
        #pragma omp parallel num_threads(N_THREADS)
        {
            #pragma omp for
            for(unsigned int j=0;j<nIn;j++)
            {
                LazyGrid *grid = new LazyGrid(resolutions.at(i));
                grid->semantic_index=j;
                map[i][j]=new NDTMap(grid,whatisthat);
                map[i][j]->guessSize(0,0,0,size[i],size[i],size[i]);
                map[i][j]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE);
            }
        }
    }
    return map;
}
void destroyMap(NDTMap ***map,unsigned int nRes,unsigned int nIn){
    for(unsigned int i=0;i<nRes;i++)
    {
        for(unsigned int j=0;j<nIn;j++)
            delete map[i][j];
        delete[] map[i];
    }
    delete[] map;
}
/*
void NDTMatch_SE::print_vals()
{
    std::ofstream fout("similarityLog.txt");
    int n_poses=hists.size();
    for(int i=0;i<n_poses;i++)
        for(int j=i+1;j<n_poses;j++)
        {
            fout<<j-i<<" ";
            fout<<pcl::geometry::distance(poses->at(j), poses->at(i))<<" ";
            fout<<hists[i]->getSimilarity(*hists[j])<<" ";
            //fout<<(gfeat[poses->at(j).label]-gfeat[poses->at(i).label]).norm()/1024<<" ";
            fout<<std::endl;
        }
    fout.close();
};
*/
