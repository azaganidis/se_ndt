cmake_minimum_required(VERSION 2.8.3)
project(se_ndt)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
## Find catkin macros and libraries
## if COMPONENTS list like find_package(catkin REQUIRED COMPONENTS xyz)
## is used, also find other catkin packages
find_package(ndt_map)
if( NOT ndt_map_FOUND)
	set(ndt_map_SOURCE_DIR "${PROJECT_SOURCE_DIR}/perception_oru_minimal/ndt_map/")
	add_subdirectory("ndt_map")
	find_package(ndt_map REQUIRED)
endif()
find_package(ndt_registration)
if( NOT ndt_registration_FOUND)
	set(ndt_registration_SOURCE_DIR
		"${PROJECT_SOURCE_DIR}/perception_oru_minimal/ndt_registration/")
	add_subdirectory("ndt_registration")
	find_package(ndt_registration REQUIRED)
endif()
find_package(catkin REQUIRED COMPONENTS
  geometry_msgs
  nav_msgs
  roscpp
  rospy
  sensor_msgs
  std_msgs
  pcl_ros
  tf
  geodesy
  rosbag
)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
## System dependencies are found with CMake's conventions
# find_package(Boost REQUIRED COMPONENTS system)
#find_package(PCL REQUIRED COMPONENTS common io)
find_package(OpenCV REQUIRED)
find_package(Eigen3 REQUIRED)

include_directories(
  include
	${catkin_INCLUDE_DIRS}
	${Boost_INCLUDE_DIRS}
	${ndt_map_INCLUDE_DIRS}
	${ndt_registration_INCLUDE_DIRS}
#	${EIGEN3_INCLUDE_DIR}
	${PCL_INCLUDE_DIRS})
link_directories ( ${Boost_LIBRARY_DIRS} )
catkin_package(
  CATKIN_DEPENDS geometry_msgs nav_msgs roscpp rospy std_msgs rosbag
  pcl_ros
  INCLUDE_DIRS include
)
#include(pcl_find_sse.cmake)
#PCL_CHECK_FOR_SSE()
#find_package(Boost REQUIRED)

set (${PROJECT_NAME}_LIB_SRCS
	src/ndt_matcher_d2d_se.cpp
	src/ndt_matcher_p2d_se.cpp
	)

add_library(${PROJECT_NAME} ${${PROJECT_NAME}_LIB_SRCS})
target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES})

add_executable(ETH_d_ci src/ETH_d_ci.cpp)
target_link_libraries(ETH_d_ci ${PROJECT_NAME} ${ndt_map_LIBRARIES}
	${ndt_registration_LIBRARIES} ${catkin_LIBRARIES} ${PCL_LIBRARIES} ${OpenCV_LIBS} ${Boost_LIBRARIES} -lboost_serialization)

add_executable(smoothness src/smoothness.cpp)
target_link_libraries(smoothness ${PROJECT_NAME} ${ndt_map_LIBRARIES}
	${ndt_registration_LIBRARIES} ${PCL_LIBRARIES} ${OpenCV_LIBS} ${Boost_LIBRARIES} -lboost_serialization)

add_executable(polesRegister src/polesRegister.cpp)
target_link_libraries(polesRegister ${PROJECT_NAME} ${ndt_map_LIBRARIES}
	${ndt_registration_LIBRARIES} ${PCL_LIBRARIES} ${OpenCV_LIBS} ${Boost_LIBRARIES} -lboost_serialization)

add_executable(polesRegister_ndt src/polesRegister_ndt.cpp)
target_link_libraries(polesRegister_ndt ${PROJECT_NAME} ${ndt_map_LIBRARIES}
	${ndt_registration_LIBRARIES} ${PCL_LIBRARIES} ${OpenCV_LIBS} ${Boost_LIBRARIES} -lboost_serialization)

add_executable(polesRegister_ndt_r src/polesRegister_ndt_r.cpp)
target_link_libraries(polesRegister_ndt_r ${PROJECT_NAME} ${ndt_map_LIBRARIES}
	${ndt_registration_LIBRARIES} ${PCL_LIBRARIES} ${OpenCV_LIBS} ${Boost_LIBRARIES} -lboost_serialization)

add_executable(plot_path src/plot_path.cpp)
target_link_libraries(plot_path ${PROJECT_NAME} ${ndt_map_LIBRARIES}
	${ndt_registration_LIBRARIES} ${PCL_LIBRARIES} ${OpenCV_LIBS} ${Boost_LIBRARIES} -lboost_serialization)

add_executable(playCloud src/playCloud.cpp)
target_link_libraries(playCloud ${PROJECT_NAME} ${ndt_map_LIBRARIES}
	${catkin_LIBRARIES} ${ndt_registration_LIBRARIES} ${PCL_LIBRARIES} ${OpenCV_LIBS} ${Boost_LIBRARIES} -lboost_serialization)

#add_executable(multimap_NDT  src/multimap_NDT.cpp)
#target_link_libraries(multimap_NDT ${PROJECT_NAME} ${ndt_map_LIBRARIES}
#	${ndt_registration_LIBRARIES} ${catkin_LIBRARIES} ${PCL_LIBRARIES} ${OpenCV_LIBS} ${Boost_LIBRARIES} -lboost_serialization)

#add_executable(multimap_NDT_fuser  src/multimap_NDT_fuser.cpp)
#message(FATAL_ERROR "BOOOST ${Boost_LIBRARIES} ${Boost_LIBRARY_DIRS}")

#add_executable(ETH_d_ci src/ETH_d_ci.cpp)
#target_link_libraries(ETH_d_ci  ${Boost_LIBRARIES} ${PCL_LIBRARIES}
#	${catkin_LIBRARIES} -lboost_serialization)

#add_executable(evaluate src/evaluate.cpp)
#target_link_libraries(evaluate ${catkin_LIBRARIES} ${Boost_LIBRARIES} ${PCL_LIBRARIES} -lboost_serialization)

#add_executable(NDT_OM_ORU src/NDT-OM_ORU.cpp)
#target_link_libraries(NDT_OM_ORU ${catkin_LIBRARIES} ${PCL_LIBRARIES} ${OpenCV_LIBS} ${Boost_LIBRARIES} -lboost_serialization)

#add_executable(NDT_ORU src/NDT_ORU.cpp)
#target_link_libraries(NDT_ORU ${catkin_LIBRARIES} ${PCL_LIBRARIES} ${OpenCV_LIBS} ${Boost_LIBRARIES} -lboost_serialization)

#target_link_libraries(multimap_NDT_fuser ${catkin_LIBRARIES} ${PCL_LIBRARIES} ${OpenCV_LIBS} ${Boost_LIBRARIES} -lboost_serialization)

#add_executable(multimap_NDT_P2D  src/multimap_NDT_P2D.cpp)
#target_link_libraries(multimap_NDT_P2D ${catkin_LIBRARIES} ${PCL_LIBRARIES} ${OpenCV_LIBS} ${Boost_LIBRARIES} -lboost_serialization)

#add_executable(cornerness src/cornerness.cpp)
#target_link_libraries(cornerness ${catkin_LIBRARIES} ${PCL_LIBRARIES} ${OpenCV_LIBS} ${Boost_LIBRARIES} -lboost_serialization)

## Uncomment this if the package has a setup.py. This macro ensures
## Uncomment this if the package has a setup.py. This macro ensures
## modules and global scripts declared therein get installed
## See http://ros.org/doc/api/catkin/html/user_guide/setup_dot_py.html
# catkin_python_setup()

################################################
## Declare ROS messages, services and actions ##
################################################

## To declare and build messages, services or actions from within this
## package, follow these steps:
## * Let MSG_DEP_SET be the set of packages whose message types you use in
##   your messages/services/actions (e.g. std_msgs, actionlib_msgs, ...).
## * In the file package.xml:
##   * add a build_depend tag for "message_generation"
##   * add a build_depend and a run_depend tag for each package in MSG_DEP_SET
##   * If MSG_DEP_SET isn't empty the following dependency has been pulled in
##     but can be declared for certainty nonetheless:
##     * add a run_depend tag for "message_runtime"
## * In this file (CMakeLists.txt):
##   * add "message_generation" and every package in MSG_DEP_SET to
##     find_package(catkin REQUIRED COMPONENTS ...)
##   * add "message_runtime" and every package in MSG_DEP_SET to
##     catkin_package(CATKIN_DEPENDS ...)
##   * uncomment the add_*_files sections below as needed
##     and list every .msg/.srv/.action file to be processed
##   * uncomment the generate_messages entry below
##   * add every package in MSG_DEP_SET to generate_messages(DEPENDENCIES ...)

## Generate messages in the 'msg' folder
# add_message_files(
#   FILES
#   Message1.msg
#   Message2.msg
# )

## Generate services in the 'srv' folder
# add_service_files(
#   FILES
#   Service1.srv
#   Service2.srv
# )

## Generate actions in the 'action' folder
# add_action_files(
#   FILES
#   Action1.action
#   Action2.action
# )

## Generate added messages and services with any dependencies listed here
# generate_messages(
#   DEPENDENCIES
#   geometry_msgs#   nav_msgs#   sensor_msgs#   std_msgs
# )

################################################
## Declare ROS dynamic reconfigure parameters ##
################################################

## To declare and build dynamic reconfigure parameters within this
## package, follow these steps:
## * In the file package.xml:
##   * add a build_depend and a run_depend tag for "dynamic_reconfigure"
## * In this file (CMakeLists.txt):
##   * add "dynamic_reconfigure" to
##     find_package(catkin REQUIRED COMPONENTS ...)
##   * uncomment the "generate_dynamic_reconfigure_options" section below
##     and list every .cfg file to be processed

## Generate dynamic reconfigure parameters in the 'cfg' folder
# generate_dynamic_reconfigure_options(
#   cfg/DynReconf1.cfg
#   cfg/DynReconf2.cfg
# )

###################################
## catkin specific configuration ##
###################################
## The catkin_package macro generates cmake config files for your package
## Declare things to be passed to dependent projects
## INCLUDE_DIRS: uncomment this if you package contains header files
## LIBRARIES: libraries you create in this project that dependent projects also need
## CATKIN_DEPENDS: catkin_packages dependent projects also need
## DEPENDS: system dependencies of this project that dependent projects also need

## Declare a C++ library
# add_library(anestis_loam_mod
#   src/${PROJECT_NAME}/anestis_loam_mod.cpp
# )

## Add cmake target dependencies of the library
## as an example, code may need to be generated before libraries
## either from message generation or dynamic reconfigure
# add_dependencies(anestis_loam_mod ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

## Declare a C++ executable
# add_executable(anestis_loam_mod_node src/anestis_loam_mod_node.cpp)

## Add cmake target dependencies of the executable
## same as for the library above
# add_dependencies(anestis_loam_mod_node ${${PROJECT_NAME}_EXPORTED_TARGETS} ${catkin_EXPORTED_TARGETS})

## Specify libraries to link a library or executable target against
# target_link_libraries(anestis_loam_mod_node
#   ${catkin_LIBRARIES}
# )

#############
## Install ##
#############

# all install targets should use catkin DESTINATION variables
# See http://ros.org/doc/api/catkin/html/adv_user_guide/variables.html

## Mark executable scripts (Python etc.) for installation
## in contrast to setup.py, you can choose the destination
# install(PROGRAMS
#   scripts/my_python_script
#   DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
# )

## Mark executables and/or libraries for installation
# install(TARGETS anestis_loam_mod anestis_loam_mod_node
#   ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
#   RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
# )

## Mark cpp header files for installation
# install(DIRECTORY include/${PROJECT_NAME}/
#   DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
#   FILES_MATCHING PATTERN "*.h"
#   PATTERN ".svn" EXCLUDE
# )

## Mark other files for installation (e.g. launch and bag files, etc.)
# install(FILES
#   # myfile1
#   # myfile2
#   DESTINATION ${CATKIN_PACKAGE_SHARE_DESTINATION}
# )

#############
## Testing ##
#############

## Add gtest based cpp test target and link libraries
# catkin_add_gtest(${PROJECT_NAME}-test test/test_anestis_loam_mod.cpp)
# if(TARGET ${PROJECT_NAME}-test)
#   target_link_libraries(${PROJECT_NAME}-test ${PROJECT_NAME})
# endif()

## Add folders to be run by python nosetests
# catkin_add_nosetests(test)
