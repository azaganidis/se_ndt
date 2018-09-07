#include <se_ndt/se_ndt.hpp>
#include <numeric>
using namespace std;
typedef Eigen::Transform<double,3,Eigen::Affine,Eigen::ColMajor> ET;

template <>
pcl::PointCloud<pcl::PointXYZI>::Ptr getCloud<pcl::PointXYZI>(string filename,char IFS, bool skip)
{
	ifstream infile(filename); // for example
	string line = "";
	pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloud(new pcl::PointCloud<pcl::PointXYZI>);
	if(skip)getline(infile, line);
	while (getline(infile, line)){
		stringstream strstr(line);
		string word = "";
		pcl::PointXYZI point;
		getline(strstr,word, IFS);
		point.x=stof(word);
		getline(strstr,word, IFS);
		point.y=stof(word);
		getline(strstr,word, IFS);
		point.z=stof(word);
		getline(strstr,word);
		point.intensity=stof(word);
		(*laserCloud).points.push_back(point);
	}
    return laserCloud;
}
template <>
pcl::PointCloud<pcl::PointXYZ>::Ptr getCloud<pcl::PointXYZ>(string filename,char IFS, bool skip)
{
	ifstream infile(filename); // for example
	string line = "";
	pcl::PointCloud<pcl::PointXYZ>::Ptr laserCloud(new pcl::PointCloud<pcl::PointXYZ>);
	if(skip)getline(infile, line);
	while (getline(infile, line)){
		stringstream strstr(line);
		string word = "";
		pcl::PointXYZ point;
		getline(strstr,word, IFS);
		point.x=stof(word);
		getline(strstr,word, IFS);
		point.y=stof(word);
		getline(strstr,word, IFS);
		point.z=stof(word);
		(*laserCloud).points.push_back(point);
	}
    return laserCloud;
}
Eigen::Affine3d readTransform(istream &infile)
{
	string line = "";
	Eigen::Affine3d T;
	T.setIdentity();
	for(int i=0;i<4;i++)
	{
		for(int j=0;j<4;j++)
		{
			string word="";
			char c=' ';
			while(c!='.'&&c!='-'&&(c<'0'||c>'9'))
			{
				c=infile.get();
			}
			while(c=='e'||c=='.'||c=='-'||(c>='0'&&c<='9'))
			{
				word+=c;
				c=infile.get();
			}
			T(i,j)=stof(word);
		}
	}
    return T;
}
size_t* sort_pointcloud(vector<double> &in,float disregard)
{
	vector<size_t> idx(in.size());
	size_t *idx2 = new size_t[in.size()];
	std::iota(idx.begin(),idx.end(),0);
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
		if(*(Tails.begin()+i)<4)
		{
			if(*(Tails.begin()+i)&1)if(I[i][p]<=minA) { minA=I[i][p];minI=i; }
			if(*(Tails.begin()+i)&2)if(I[i][p]>=maxA) { maxA=I[i][p];maxI=i; }
		}
	}
	return number_points-maxA<minA?maxI+num:minI;
}
inline bool checkInLimits(size_t **in,int p,int num,int cu,int cl)
{
	for(int i=0;i<num;i++)
		if(in[i]!=NULL)
			if(in[i][p]>cu||in[i][p]<cl)
				return true;
	return false;
}

vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> NDTMatch_SE::getSegments(pcl::PointCloud<pcl::PointXYZ>::Ptr laserCloudIn,initializer_list<vector<double> >& attributes_,initializer_list<int > distribution_tails_,initializer_list<float> disregard_, float rejectPerc)
{
	vector<vector<double> > attributes(attributes_);
	vector<float> disregard(disregard_);
	vector<int> distribution_tails(distribution_tails_);
	int cloudSize = laserCloudIn->points.size();
	int num_attributes=attributes.size();
	size_t *sorted_indexes[num_attributes];
	for(auto i=0;i< num_attributes;i++)
	{
		if(distribution_tails[i]<4)
			sorted_indexes[i]=sort_pointcloud(attributes.at(i),disregard.at(i));
		else sorted_indexes[i]=NULL;
	}

	int cut_off_l=(1-rejectPerc)/2*cloudSize;
	int cut_off_u=(1+rejectPerc)/2*cloudSize;
	////////
	//Create look-up
	int look_up[2*num_attributes];
	size_t num_tails=0;
	for(int i=0;i<num_attributes;i++)
	{
		if(distribution_tails[i]<4)
		{
			if(distribution_tails[i]==1||distribution_tails[i]==3)
			{
				look_up[i]= num_tails;
				num_tails++;
			}else look_up[i]=-1;
			if(distribution_tails[i]==2)
			{
				look_up[i+num_attributes]=num_tails;
				num_tails++;
			}else look_up[i+num_attributes]=-1;
		}
		else if(distribution_tails[i]==(int )'e'||distribution_tails[i]=='>'||distribution_tails[i]=='<'||distribution_tails[i]=='*'||distribution_tails[i]=='=')
		{
			look_up[i]=num_tails;
			num_tails++;
			look_up[i+num_attributes]=-1;
		}
		else look_up[i]=-1;
	}
	///////
	vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud;
	for(int i=0;i<num_tails+semantic_labels.size();i++)
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT(new pcl::PointCloud<pcl::PointXYZ>);
		laserCloud.push_back(cloudT);
	}
	for(int i=0;i<cloudSize;i++)
	{
		int index=-1;
		for(auto j=0;j< num_attributes;j++)
		{
			if(distribution_tails.at(j)==(int )'e'||distribution_tails.at(j)==(int )'=')
			{
				if(attributes.at(j).at(i)==disregard.at(j))
				{
					index=look_up[j];
					continue;
				}
			}
			else if(distribution_tails.at(j)==(int )'>')
			{
				if(attributes.at(j).at(i)>=disregard.at(j))
				{
					index=look_up[j];
					continue;
				}
			}
			else if(distribution_tails.at(j)==(int )'<')
			{
				if(attributes.at(j).at(i)<=disregard.at(j))
				{
					index=look_up[j];
					continue;
				}
			}
			else if (distribution_tails.at(j)==117&&attributes.at(j).at(i)!=disregard.at(j))
			{
				std::vector<int>::iterator my_iter;
				int sem_val=attributes.at(j).at(i);
				my_iter=find(semantic_labels.begin(),semantic_labels.end(),sem_val);
				if(semantic_labels.end()==my_iter)
				{
					pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT(new pcl::PointCloud<pcl::PointXYZ>);
					laserCloud.push_back(cloudT);
					semantic_labels.push_back(sem_val);
					my_iter=semantic_labels.end()-1;
				}
				int nummm=distance(semantic_labels.begin(),my_iter);
				index=num_tails+distance(semantic_labels.begin(),my_iter);
				if(index>=NumInputs)cerr<<"too many labels"<<endl;
				continue;
			}
			else if(distribution_tails.at(j)==(int )'*'&&index==-1)
				index=look_up[j];
		}


		if(index==-1&&!checkInLimits(sorted_indexes,i,num_attributes,cut_off_u,cut_off_l)) continue;
		pcl::PointXYZ point;
		point.x=laserCloudIn->points[i].x;
		point.y=laserCloudIn->points[i].y;
		point.z=laserCloudIn->points[i].z;
		//point.intensity=laserCloudIn->points[i].intensity;
		if(index==-1)index=look_up[index_selector(sorted_indexes, i,num_attributes,distribution_tails,cloudSize)];
		if(index!=-1)
			laserCloud[index]->points.push_back(point);
	}
	for(auto i=0;i<num_attributes;i++)
		delete[] sorted_indexes[i];


	return laserCloud;
}
size_t count_tails(vector<int>& distribution_tails)
{
	size_t number_tails0=count(distribution_tails.begin(),distribution_tails.end(),0);
	size_t number_tails3=count(distribution_tails.begin(),distribution_tails.end(),3);
	return distribution_tails.size()+number_tails3-number_tails0;
}

lslgeneric::NDTMap **initMap(int number_tails,initializer_list<float> resolutions_, initializer_list<float>size_)
{
	vector<float> size(size_),resolutions(resolutions_);
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
//#include <chrono>
//#include <iomanip>
void loadMap(lslgeneric::NDTMap **map,std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> input_clouds,float sensor_range)
{
//    cerr<<std::setprecision(3);
//    double one=0, two=0;
	for(size_t i=0;i<input_clouds.size();i++)
	{
//        auto start = std::chrono::high_resolution_clock::now();
		map[i]->loadPointCloud(*input_clouds[i],sensor_range);
//        std::chrono::duration<double> elapsed1 = std::chrono::high_resolution_clock::now()- start;
//        start = std::chrono::high_resolution_clock::now();
		map[i]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE);
//        std::chrono::duration<double> elapsed2 = std::chrono::high_resolution_clock::now()- start;
//        one+=elapsed1.count();
//        two+=elapsed2.count();
	}
//    cerr<<std::fixed<<std::setprecision(3)<<one<< " "<<two<<endl;
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
	NumInputs=count_tails(tails_t)+70*std::count(tails_t.begin(),tails_t.end(),'u');
	firstRun=true;

		matcher.NumInputs=NumInputs;
		matcher.ITR_MAX =max_iter;
		matcher.step_control=true;
		matcher.n_neighbours=2;

	map=new lslgeneric::NDTMap ** [resolutions.size()];
	mapLocal=new lslgeneric::NDTMap ** [resolutions.size()];
	for(unsigned int i=0;i<resolutions.size();i++)
	{
		map[i]=initMap(NumInputs,{resolutions.at(i)},size);
		mapLocal[i]=initMap(NumInputs,{resolutions.at(i)},size);
	}
}
#ifdef VISUALIZE
void NDTMatch_SE::visualize()
{
    thread t1(&NDTMatch_SE::visualize_thread, this);
    t1.detach();
}
void NDTMatch_SE::visualize_thread()
{
    if(viewer==NULL)
    {
        viewer = new NDTViz(true);
        viewer->win3D->start_main_loop_own_thread();
    }
    float occupancy=128;
    viewer->plotNDTSAccordingToOccupancy(occupancy, map,std::vector<int>{0}, std::vector<int>{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14});
    while(viewer->win3D->isOpen())
    {
         usleep(10000);
        if (viewer->win3D->keyHit())
        {
            int key = viewer->win3D->getPushedKey();
            switch (key){
                case '+':occupancy++;break;
                case '-':occupancy--;break;
                case '*':occupancy*=2;break;
                case '/':occupancy/=2;break;
                case 'q':delete viewer;break;
            }
            viewer->plotNDTSAccordingToOccupancy(occupancy, map,std::vector<int>{0}, std::vector<int>{0,1,2,3,4,5,6,7,8,9,10,11,12,13,14});
            viewer->win3D->keyHitReset();
        }
    }
}
#endif
Eigen::Affine3d NDTMatch_SE::match(Eigen::Affine3d Tinit, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2,initializer_list<vector<double> > attributes1,initializer_list<vector<double> > attributes2)
{
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud1=getSegments(cloud1,attributes1,tails,ignore,removeProbability);
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud2=getSegments(cloud2,attributes2,tails,ignore,removeProbability);
#define USE_OMP
#define n_threads 12
#pragma omp parallel num_threads(n_threads)
{
    #pragma omp for
	for(unsigned int i=0;i<resolutions.size();i++)
	{
		loadMap(map[i],laserCloud1);
		loadMap(mapLocal[i],laserCloud2);
	}
}
	for(auto i:resolutions_order)
	{
		matcher.current_resolution=resolutions.at(i);
		matcher.match(map[i],mapLocal[i],Tinit,true);
	}
	//std::cout<<getHes(matcher.HessianF,matcher.score_gradientF).inverse()<<std::endl;
	return Tinit;
}
vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> NDTMatch_SE::getSegmentsFast(pcl::PointCloud<pcl::PointXYZ>::Ptr laserCloudIn,vector<double> &attributes)
{
	int cloudSize = laserCloudIn->points.size();
	int num_attributes=ignore.size();

	vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud;
	for(int i=0;i<num_attributes;i++)
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT(new pcl::PointCloud<pcl::PointXYZ>);
		laserCloud.push_back(cloudT);
	}
	for(int i=0;i<cloudSize;i++)
	{
            int atr = int(attributes.at(i));
            pcl::PointXYZ point;
            point.x=laserCloudIn->points[i].x;
            point.y=laserCloudIn->points[i].y;
            point.z=laserCloudIn->points[i].z;
            laserCloud[atr]->points.push_back(point);
    }
	return laserCloud;
}
vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> NDTMatch_SE::getSegmentsFast(pcl::PointCloud<pcl::PointXYZI>::Ptr laserCloudIn)
{
	int cloudSize = laserCloudIn->points.size();
	int num_attributes=ignore.size();

	vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud;
	for(int i=0;i<num_attributes;i++)
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloudT(new pcl::PointCloud<pcl::PointXYZ>);
		laserCloud.push_back(cloudT);
	}
	for(int i=0;i<cloudSize;i++)
	{
            auto pnt=laserCloudIn->points[i];
            int atr = round(pnt.intensity);
            pcl::PointXYZ point;
            point.x=pnt.x;
            point.y=pnt.y;
            point.z=pnt.z;
            laserCloud[atr]->points.push_back(point);
    }
	return laserCloud;
}
Eigen::Affine3d NDTMatch_SE::matchFast(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2,vector<double> &attributes1,vector<double> &attributes2)
{
    Eigen::Affine3d Tinit;
    Tinit.setIdentity();
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud1=getSegmentsFast(cloud1,attributes1);
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud2=getSegmentsFast(cloud2,attributes2);
#define USE_OMP
#define n_threads 12
#pragma omp parallel num_threads(n_threads)
{
    #pragma omp for
	for(unsigned int i=0;i<resolutions.size();i++)
	{
		loadMap(map[i],laserCloud1);
		loadMap(mapLocal[i],laserCloud2);
	}
}
	for(auto i:resolutions_order)
	{
		matcher.current_resolution=resolutions.at(i);
		matcher.match(map[i],mapLocal[i],Tinit,true);
	}
	//std::cout<<getHes(matcher.HessianF,matcher.score_gradientF).inverse()<<std::endl;
	return Tinit;
}

Eigen::Affine3d NDTMatch_SE::matchFaster(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, vector<double> &attributes)
{
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud=getSegmentsFast(cloud,attributes);
    Eigen::Affine3d Tinit;
    Tinit.setIdentity();
#pragma omp parallel num_threads(8)
{
    #pragma omp for
	for(unsigned int i=0;i<resolutions.size();i++)
		loadMap(mapLocal[i],laserCloud);
}
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
Eigen::Affine3d NDTMatch_SE::matchFaster(Eigen::Affine3d Tinit, pcl::PointCloud<pcl::PointXYZI>::Ptr cloud)
{
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud=getSegmentsFast(cloud);
	for(unsigned int i=0;i<resolutions.size();i++)
		loadMap(mapLocal[i],laserCloud);
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
Eigen::Affine3d NDTMatch_SE::matchFaster(pcl::PointCloud<pcl::PointXYZI>::Ptr cloud)
{
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud=getSegmentsFast(cloud);
    Eigen::Affine3d Tinit;
    Tinit.setIdentity();

	for(unsigned int i=0;i<resolutions.size();i++)
		loadMap(mapLocal[i],laserCloud);
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
Eigen::Affine3d NDTMatch_SE::matchFaster_OM(Eigen::Affine3d Tinit, pcl::PointCloud<pcl::PointXYZI>::Ptr cloud)
{
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud=getSegmentsFast(cloud);
	for(unsigned int i=0;i<resolutions.size();i++)
		loadMap(mapLocal[i],laserCloud);
	if(!firstRun)
	{
		for(auto i:resolutions_order)
		{
			matcher.current_resolution=resolutions.at(i);
			matcher.match(map[i],mapLocal[i],Tinit,true);
		}
	}
	else firstRun=false;

    vector<thread> tc;
    for(size_t i=0;i<laserCloud.size();i++)
        tc.push_back( std::thread([i,&Tinit, laserCloud](){
                    lslgeneric::transformPointCloudInPlace(Tinit, *laserCloud[i]);
                }));
    for(auto& t:tc)t.join();
    for(size_t i=0;i<laserCloud.size();i++)
    {
		for(auto rez:resolutions_order)
		{
            map[rez][i]->addPointCloud(Tinit.translation(), *laserCloud[i],0.06, 2.5);
            map[rez][i]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE,1e5,255,Tinit.translation(),0.1);
        }
    }
	return Tinit;
}
struct MatchPathSeparator
{
    bool operator()( char ch ) const
    {
        return ch == '/';
    }
};
std::string basename(std::string const& pathname)
{
	return std::string(std::find_if(pathname.rbegin(),pathname.rend(),MatchPathSeparator() ).base(),pathname.end());
}
Eigen::Affine3d NDTMatch_SE::match(Eigen::Affine3d Tinit, std::string cloudF1, std::string cloudF2,initializer_list<vector<double> > attributes1,initializer_list<vector<double> > attributes2)//Eigen::Affine3d Tinit, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2)
{
	for(unsigned int i=0;i<resolutions.size();i++)
	{
		for(unsigned int j=0;j<NumInputs;j++)
		{
			std::string fname=cloudF1;
			std::string fname_jff=precomputed_ndt_folder+
				basename(fname+"."+std::to_string(resolutions[i])+"."+std::to_string(j)+".jff");
			std::ifstream file((fname_jff).c_str());
			if(file.good()&&useSaved)
			{
				/*
				if(map!=NULL)
				{
				map[i][j]->unsetFirstLoad();
				delete map[i][j];
				}
				lslgeneric::LazyGrid *grid = new lslgeneric::LazyGrid(resolutions[i]);
				map[i][j]=new lslgeneric::NDTMap(grid);
				*/
				map[i][j]->loadFromJFF((fname_jff).c_str());
			}
			else
			{
				std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >
					laserClouds=getSegments(getCloud<pcl::PointXYZ>(fname,IFS,skip),attributes1,tails,ignore,removeProbability);
				for(j=j;j<NumInputs;j++)//intentional the same variable, load the remaining
				{
					map[i][j]->loadPointCloud(*laserClouds[j],sensor_range);
					map[i][j]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE);
					fname_jff=precomputed_ndt_folder+
						basename(fname+"."+std::to_string(resolutions[i])+"."+std::to_string(j)+".jff");
					map[i][j]->writeToJFF(fname_jff.c_str());
				}
			}
		}
		for(unsigned int j=0;j<NumInputs;j++)
		{
			std::string fname=cloudF2;
			std::string fname_jff=precomputed_ndt_folder+
				basename(fname+"."+std::to_string(resolutions[i])+"."+std::to_string(j)+".jff");
			std::ifstream file((fname_jff).c_str());
			if(file.good()&&useSaved)
			{
				/*
				if(mapLocal!=NULL)
				{
				mapLocal[i][j]->unsetFirstLoad();
				delete mapLocal[i][j];
				}
				lslgeneric::LazyGrid *grid = new lslgeneric::LazyGrid(resolutions[i]);
				mapLocal[i][j]=new lslgeneric::NDTMap(grid);
				*/
				mapLocal[i][j]->loadFromJFF((fname_jff).c_str());
			}
			else
			{
				std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >
					laserClouds=getSegments(getCloud<pcl::PointXYZ>(fname,IFS,skip),attributes2,tails,ignore,removeProbability);
				for(j=j;j<NumInputs;j++)//intentional the same variable, load the remaining
				{
					mapLocal[i][j]->loadPointCloud(*laserClouds[j],sensor_range);
					mapLocal[i][j]->computeNDTCells(CELL_UPDATE_MODE_SAMPLE_VARIANCE);
					fname_jff=precomputed_ndt_folder+
						basename(fname+"."+std::to_string(resolutions[i])+"."+std::to_string(j)+".jff");
					mapLocal[i][j]->writeToJFF(fname_jff.c_str());
				}
			}
		}
	}
	for(auto i:resolutions_order)
	{
		matcher.current_resolution=resolutions.at(i);
		matcher.match(map[i],mapLocal[i],Tinit,true);
	}
	//std::cout<<getHes(matcher.HessianF,matcher.score_gradientF).inverse()<<std::endl;
	return Tinit;
}
Eigen::Matrix<double, 6,6> NDTMatch_SE::getPoseCovariance(Eigen::Affine3d T)
{
	Eigen::MatrixXd Covariance(6,6);
	Eigen::Matrix<double,6,6> Covariance_;
	matcher.covariance(map,mapLocal,T,resolutions,Covariance);
	Covariance_=Covariance;
	return Covariance_;
}
Eigen::Affine3d NDTMatch_SE::match(Eigen::Affine3d Tinit, pcl::PointCloud<pcl::PointXYZ>::Ptr cloud, initializer_list<vector<double> > attributes)
{
	std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr >laserCloud=getSegments(cloud,attributes,tails,ignore,removeProbability);
	for(unsigned int i=0;i<resolutions.size();i++)
		loadMap(mapLocal[i],laserCloud);
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
Eigen::Affine3d NDTMatch_SE::match(pcl::PointCloud<pcl::PointXYZ>::Ptr cloud1,pcl::PointCloud<pcl::PointXYZ>::Ptr cloud2, initializer_list<vector<double> > attributes1, initializer_list<vector<double> > attributes2)
{
	Eigen::Affine3d T;
	T.setIdentity();
	return this->match(T,cloud1,cloud2,attributes1,attributes2);
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
