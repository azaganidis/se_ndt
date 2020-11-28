# Semantic assisted Normal Distributions Transform 
Point cloud registration and loop closure detection.

## How to run the code.
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
To visualize the map, press any key in the NDTVizGlut window. Exit with `q`.
The mapping results are written in `pose_graph_out.txt`.

-----
## Relevant papers:
For registration and loop closure detection:
```
@inproceedings{Zaganidis2019,
  title={{Semantic Assisted Loop Closure in SLAM using NDT Histograms}},
  author={Zaganidis, Anestis and Zerntev, Alexandros and Duckett, Tom and Cielniak, Grzegorz},
  year={2019},
  publisher={IEEE},
  booktitle = {2019 IEEE/RSJ Int. Conf. on Intelligent Robots and Syst.}
}
@article{Zaganidis2018,
  author={A. {Zaganidis} and L. {Sun} and T. {Duckett} and G. {Cielniak}},
  journal={IEEE Robotics and Automation Letters},
  title={{Integrating Deep Semantic Segmentation Into 3-D Point Cloud Registration}},
  year={2018},
  volume={3},
  number={4},
  pages={2942-2949},
  doi={10.1109/LRA.2018.2848308},
  ISSN={2377-3766},
  month={Oct},
}
@inproceedings{Zaganidis2017,
  title={{Semantic-assisted 3D Normal Distributions Transform for scan registration in environments with limited structure}},
  author={Zaganidis, Anestis and Magnusson, Martin and Duckett, Tom and Cielniak, Grzegorz},
  year={2017},
  publisher={IEEE},
  booktitle = {2017 IEEE/RSJ Int. Conf. on Intelligent Robots and Syst.}
}
```
For the classifier:
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
