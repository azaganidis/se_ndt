FROM ros:melodic-perception
RUN apt-get update && apt-get install -y \
    libeigen3-dev \
    libpcl-dev \
    ros-melodic-tf-conversions \
    ros-melodic-tf2-ros \
    ros-melodic-roscpp \
    ros-melodic-pcl-ros \
    screen \
    vim

RUN /bin/bash -c 'source $HOME/.bashrc;source /ros_entrypoint.sh;mkdir -p /root/catkin_ws/src;cd /root/catkin_ws/src;catkin_init_workspace'
RUN /bin/bash -c 'source $HOME/.bashrc;source /ros_entrypoint.sh; cd /root/catkin_ws/src; git clone https://github.com/tradr-project/tensorflow_ros_cpp;git clone https://github.com/ethz-asl/tensorflow_catkin;git clone https://github.com/catkin/catkin_simple;cd ..;catkin_make'
#up to here installation of ros and tensorflow
COPY se_ndt /root/catkin_ws/src/se_ndt/
RUN /bin/bash -c 'source $HOME/.bashrc;source /ros_entrypoint.sh; cd /root/catkin_ws/;catkin_make'
RUN echo source ~/catkin_ws/devel/setup.bash >> /root/.bashrc
WORKDIR /root/
CMD ["bash"]
