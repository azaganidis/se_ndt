Dependencies: Boost, Eigen, PCL, Ceres solver

Build with:
```bash
git clone --recurse-submodules https://github.com/azaganidis/se_ndt
mkdir se_ndt/build
cd se_ndt/build
cmake -DCMAKE_BUILD_TYPE=Release -DWITH_GL=True ..
make -j4
```

Example use on [KITTI](http://www.cvlibs.net/download.php?file=data_odometry_velodyne.zip) using the labels from [RangeNet++](http://www.ipb.uni-bonn.de/wp-content/papercite-data/pdf/milioto2019iros.pdf) [download](http://www.ipb.uni-bonn.de/html/projects/bonnetal/lidar/semantic/predictions/darknet53-knn.tar.gz)
```bash
./kitti_slam -p /mnt/external/Datasets/kitti/sequences/00/velodyne/* -l /mnt/external/Datasets/kitti/darknet53-knn/darknet53-knn/sequences/00/predictions/* -v
```

Relevant papers:
```
@inproceedings{milioto2019iros,
  author    = {A. Milioto and I. Vizzo and J. Behley and C. Stachniss},
  title     = {{RangeNet++: Fast and Accurate LiDAR Semantic Segmentation}},
  booktitle = {IEEE/RSJ Intl.~Conf.~on Intelligent Robots and Systems (IROS)},
  year      = 2019,
  codeurl   = {https://github.com/PRBonn/lidar-bonnetal},
  videourl  = {https://youtu.be/wuokg7MFZyU},
}
```
