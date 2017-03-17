# Catkin User Guide: http://www.ros.org/doc/groovy/api/catkin/html/user_guide/user_guide.html
# Catkin CMake Standard: http://www.ros.org/doc/groovy/api/catkin/html/user_guide/standards.html
cmake_minimum_required(VERSION 2.8.3)
#set(CMAKE_BUILD_TYPE RelWithDebug) 
set(CMAKE_BUILD_TYPE Release) 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O0")
project(ndt_map)
# Load catkin and all dependencies required for this package

find_package(cmake_modules REQUIRED )
find_package(Eigen3 REQUIRED )
find_package(catkin REQUIRED COMPONENTS cv_bridge roscpp pcl_ros pcl_conversions 
  std_msgs
  nav_msgs
   message_generation  )

add_message_files(
  FILES
  NDTMapMsg.msg   
  NDTCellMsg.msg
)

generate_messages(
  DEPENDENCIES geometry_msgs std_msgs
)

catkin_package(
    CATKIN_DEPENDS cv_bridge message_runtime roscpp pcl_ros pcl_conversions 
    INCLUDE_DIRS include
    LIBRARIES ${PROJECT_NAME}
)

set (${PROJECT_NAME}_LIB_SRCS
	src/ndt_cell.cpp
	src/ndt_map.cpp
	src/cell_vector.cpp
	src/lazy_grid.cpp
	src/ndt_map_hmt.cpp
	src/ndt_histogram.cpp
#	src/oc_tree.cpp
#	src/adaptive_oc_tree.cpp
    )
add_executable(histTest test/test_ndt_histogram.cc)
add_library(${PROJECT_NAME} ${${PROJECT_NAME}_LIB_SRCS})
#add_executable(jffSaveTest test/jffLazyGridTest_saving.cc)
#add_executable(jffLoadTest test/jffLazyGridTest_loading.cc)
add_executable(test_map_topic test/test_map_topic.cpp)
add_executable(test_map_topic_1 test/test_map_topic_1.cpp)
add_executable(ndt_builder test/ndtMapBuilder.cc)
add_executable(test_occ_map_topic test/test_occ_map_topic.cpp)

find_package(OpenCV REQUIRED)
include_directories(include ${catkin_INCLUDE_DIRS})
include_directories(include ${EIGEN_INCLUDE_DIRS})
include_directories(include ${PCL_INCLUDE_DIRS})
add_definitions(${EIGEN_DEFINITIONS})
link_directories(${catkin_LIBRARY_DIRS})

#uncomment if you have defined messages

#uncomment if you have defined services
#add_service_files(
#  FILES
  # TODO: List your msg files here
#)

#make a node that listens to point clouds and makes ndts

#target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES})
target_link_libraries(histTest ${PROJECT_NAME} ${catkin_LIBRARIES})
#target_link_libraries(jffSaveTest ${catkin_LIBRARIES})
#target_link_libraries(jffLoadTest ${catkin_LIBRARIES})
target_link_libraries(test_map_topic ${PROJECT_NAME} ${catkin_LIBRARIES})
target_link_libraries(test_map_topic_1 ${PROJECT_NAME} ${catkin_LIBRARIES})
target_link_libraries(ndt_builder ${PROJECT_NAME} ${catkin_LIBRARIES})
target_link_libraries(test_occ_map_topic ${PROJECT_NAME} ${catkin_LIBRARIES})

add_dependencies(${PROJECT_NAME} ${PROJECT_NAME}_generate_messages_cpp)
target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES})

#add_executable(likelihood_test test/ndtLikelihoodTester.cc)
#target_link_libraries(likelihood_test ${PROJECT_NAME}  pointcloud_vrml)
##add_executable(ltest test/likelihoodSingleScan.cc)
#add_executable(batchTestHist test/batchTesterHistogramFPFH.cc)
#target_link_libraries(batchTestHist ${PROJECT_NAME} pointcloud_vrml)
## Generate added messages and services with any dependencies listed here

install(TARGETS ${PROJECT_NAME}
        ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
        LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
        RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
       )

install(DIRECTORY include/${PROJECT_NAME}/
	DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
       )

