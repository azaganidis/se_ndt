#include "pcl/point_cloud.h"
#include <boost/program_options.hpp>
#include "pcl/io/pcd_io.h"

namespace po = boost::program_options;
using namespace std;

void parseLogFile(std::string filename, pcl::PointCloud<pcl::PointXYZ> &pcloud) {
    FILE *fin = fopen (filename.c_str(), "r");
    if(fin == NULL) return;
    char *line = NULL;
    size_t len;
    getline(&line,&len,fin);
    //ignore the first line of gibberish
    free(line);
    double t,in;
    int ids,idp;
    pcl::PointXYZ pt;

    pcloud.points.clear();
    while(fscanf(fin,"%lf,%f,%f,%f,%lf,%d,%d",&t,&pt.x,&pt.y,&pt.z,&in,&ids,&idp) > 0) {
	pcloud.points.push_back(pt);
    }
    pcloud.width=1;
    pcloud.height=pcloud.points.size();
    //pcloud.dense = false;

    fclose(fin);
}

int main(int argc, char **argv) {
    std::string base_name, point_type;
    double resolution;
    int n_scans;
    bool doFuser;

    po::options_description desc("Allowed options");
    desc.add_options()
	("help", "produce help message")
	("base-name", po::value<string>(&base_name), "pre-string of the files to be loaded, relative to current working directory")
	("number", po::value<int>(&n_scans)->default_value(1), "number of scans to process")
	("fuser", po::value<bool>(&doFuser)->default_value(false), "run ndt fuser instead")
	("resolution", po::value<double>(&resolution)->default_value(1.), "resolution of the map for ndt fuser")
	;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);

    if (!vm.count("base-name") || !vm.count("number"))
    {
	cout << "Missing arguments.\n";
	cout << desc << "\n";
	return 1;
    }
    if (vm.count("help"))
    {
	cout << desc << "\n";
	return 1;
    }

    pcl::PointCloud<pcl::PointXYZ> cloud;

    for(int i=0; i<n_scans; ++i) {
	char buf[1000];
	snprintf(buf,999,"%s%d.csv",base_name.c_str(),i);
	std::cerr<<"parsing file "<<buf<<std::endl;
	parseLogFile(buf,cloud);
	snprintf(buf,999,"%s%d.pcd",base_name.c_str(),i);
	pcl::io::savePCDFileBinary (buf, cloud);

    }
}

