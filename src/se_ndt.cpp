#include <se_ndt/se_ndt.hpp>
#include <numeric>
using namespace std;
typedef Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> ET;
size_t* sort_pointcloud(vector<double> &in,float disregard)
{
	vector<size_t> idx(in.size());
	size_t *idx2 = new size_t[in.size()];
	iota(idx.begin(),idx.end(),0);
	sort(idx.begin(), idx.end(),[&in](size_t i1,size_t i2){return in[i1]<in[i2];});
	std::vector<size_t>::iterator up= std::upper_bound (idx.begin(), idx.end(), disregard,[&in](double v,size_t i1){return in[i1]>v;}); //                   
	size_t mid = (idx.end()-up)/2;
	  std::swap_ranges(up, up+mid, idx.begin());
	  int j=0;
	  for(auto i:idx)
		  idx2[i]=j++;
	return idx2;
}
inline size_t index_selector(size_t **I,int p,int num,std::vector<int> Tails,size_t number_points)
{
	size_t minI=0,minA=I[0][p];
	size_t maxI=0,maxA=I[0][p];
	for(auto i=0;i<num;i++)
	{
		if(*(Tails.begin()+i)&1)if(I[i][p]<=minA) { minA=I[i][p];minI=i; }
		if(*(Tails.begin()+i)&2)if(I[i][p]>=maxA) { maxA=I[i][p];maxI=i; }
	}
	return number_points-maxA<minA?maxI+num:minI;
}
inline bool checkInLimits(size_t **in,int p,int num,int cu,int cl)
{
	for(int i=0;i<num;i++)
		if(in[i][p]>cu||in[i][p]<cl)
			return true;
	return false;
}
vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> getSegments(pcl::PointCloud<pcl::PointXYZ>::Ptr laserCloudIn,initializer_list<vector<double> >& attributes_,initializer_list<int > distribution_tails_,initializer_list<float> disregard_, float rejectPerc)
{
	vector<vector<double> > attributes(attributes_);
	vector<float> disregard(disregard_);
	vector<int> distribution_tails(distribution_tails_);
	int cloudSize = laserCloudIn->points.size();
	int num_attributes=attributes.size();
	size_t *sorted_indexes[num_attributes];
	for(auto i=0;i< num_attributes;i++)
		sorted_indexes[i]=sort_pointcloud(attributes.at(i),disregard.at(i));

	int cut_off_l=(1-rejectPerc)/2*cloudSize;
	int cut_off_u=(1+rejectPerc)/2*cloudSize;
	////////
	//Create look-up
	size_t look_up[2*num_attributes];
	size_t num_tails=0;
	for(int i=0;i<num_attributes;i++)
	{
		if(distribution_tails[i]&1)
		{
			look_up[2*i]= num_tails;
			num_tails++;
		}else look_up[2*i]=-1;
		if(distribution_tails[i]&2)
		{
			look_up[2*i+1]=num_tails;
			num_tails++;
		}else look_up[2*i+1]=-1;
	}
	///////
	vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud;
	for(int i=0;i<num_tails;i++)
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT(new pcl::PointCloud<pcl::PointXYZ>);
		laserCloud.push_back(cloudT);
	}
	for(int i=0;i<cloudSize;i++)
	{
		if(!checkInLimits(sorted_indexes,i,num_attributes,cut_off_u,cut_off_l)) continue;
		pcl::PointXYZ point;
		point.x=laserCloudIn->points[i].x;
		point.y=laserCloudIn->points[i].y;
		point.z=laserCloudIn->points[i].z;
		//point.intensity=laserCloudIn->points[i].intensity;
		int index=look_up[index_selector(sorted_indexes, i,num_attributes,distribution_tails,cloudSize)];
		if(index!=-1)
			laserCloud[index]->points.push_back(point);
	}
	for(auto i=0;i<num_attributes;i++)
		delete[] sorted_indexes[i];


	return laserCloud;
}
size_t count_tails(vector<int>& distribution_tails)
{
	size_t number_tails1=count(distribution_tails.begin(),distribution_tails.end(),1);
	size_t number_tails2=count(distribution_tails.begin(),distribution_tails.end(),2);
	size_t number_tails3=count(distribution_tails.begin(),distribution_tails.end(),3);
	return number_tails1+number_tails2+2*number_tails3;
}

lslgeneric::NDTMap **initMap(initializer_list<int> distribution_tails_,initializer_list<float> resolutions_, initializer_list<float>size_)
{
	vector<int> distribution_tails(distribution_tails_);
	vector<float> size(size_),resolutions(resolutions_);
	size_t number_tails=count_tails(distribution_tails);
	if(number_tails!=resolutions.size()&&resolutions.size()!=1)
		cerr<<"Number of resolutions different than number of segments. Taking resolution index modulus."<<endl;
	if(size.size()!=3)
	{
		float max_size=*max_element(size.begin(),size.end());
		size.resize(3);
		size[0]=max_size;
		size[1]=max_size;
		size[2]=max_size;
		cerr<<"Wrong size parameter. Using size x=y=z="<<max_size<<endl;
	}
	lslgeneric::NDTMap **map = new lslgeneric::NDTMap * [number_tails];
	for(size_t i=0;i<number_tails;i++)
	{
		lslgeneric::LazyGrid *grid = new lslgeneric::LazyGrid(resolutions[i%resolutions.size()]);
		lslgeneric::NDTMap *mapTMP=new lslgeneric::NDTMap(grid);
		map[i]=mapTMP;
		map[i]->guessSize(0,0,0,size[0],size[1],size[2]);
		map[i]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE);
	}
	return map;
}
void loadMap(lslgeneric::NDTMap **map,std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> input_clouds,size_t number_tails,float sensor_range)
{
	for(size_t i=0;i<number_tails;i++)
	{
		map[i]->loadPointCloud(*input_clouds[i],sensor_range);
		map[i]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE);
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
NDTMatch_SE::NDTMatch_SE(initializer_list<float> b,initializer_list<int> c,initializer_list<float> d,initializer_list<int> e,initializer_list<float> ig,float removeP,int max_iter):resolutions(b),resolutions_order(c),size(d),tails(e),ignore(ig),removeProbability(removeP)
{
	vector<int> tails_t(tails);
	NumInputs=count_tails(tails_t);
	firstRun=true;

		matcher.NumInputs=NumInputs;
		matcher.ITR_MAX =max_iter;
		matcher.step_control=true;

	map=new lslgeneric::NDTMap ** [resolutions.size()];
	mapLocal=new lslgeneric::NDTMap ** [resolutions.size()];
	for(auto i=0;i<resolutions.size();i++)
	{
		map[i]=initMap(tails,{resolutions.at(i)},size);
		mapLocal[i]=initMap(tails,{resolutions.at(i)},size);
	}
}
Eigen::Affine3d NDTMatch_SE::match(Eigen::Affine3d Tinit, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2,initializer_list<vector<double> > attributes)
{
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud1=getSegments(cloud1,attributes,tails,ignore,removeProbability);
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud2=getSegments(cloud2,attributes,tails,ignore,removeProbability);
	for(int i=0;i<resolutions.size();i++)
	{
		loadMap(map[i],laserCloud1,NumInputs);
		loadMap(mapLocal[i],laserCloud2,NumInputs);
	}
	for(auto i:resolutions_order)
	{
		matcher.current_resolution=resolutions.at(i);
		matcher.match(map[i],mapLocal[i],Tinit,true);
	}
	//std::cout<<getHes(matcher.HessianF,matcher.score_gradientF).inverse()<<std::endl;
	return Tinit;
}
Eigen::Affine3d NDTMatch_SE::match(Eigen::Affine3d Tinit, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, initializer_list<vector<double> > attributes)
{
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud=getSegments(cloud,attributes,tails,ignore,removeProbability);
	for(int i=0;i<resolutions.size();i++)
		loadMap(mapLocal[i],laserCloud,NumInputs);
	if(!firstRun)
	{
		for(auto i:resolutions_order)
		{
			matcher.current_resolution=resolutions.at(i);
			matcher.match(map[i],mapLocal[i],Tinit,true);
		}
		//std::cout<<getHes(matcher.HessianF,matcher.score_gradientF).inverse()<<std::endl;
	}
	else firstRun=false;
	lslgeneric::NDTMap ***mapT;
	mapT=map;
	map=mapLocal;
	mapLocal=mapT;
	return Tinit;
}
Eigen::Affine3d NDTMatch_SE::match(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, initializer_list<vector<double> > attributes)
{
	Eigen::Affine3d T;
	T.setIdentity();
	return this->match(T,cloud,attributes);
}
