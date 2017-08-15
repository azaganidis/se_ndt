cmake_minimum_required(VERSION 2.8.3)
project(se_ndt)
find_package(OpenMP REQUIRED)
find_package(OpenCV REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#find_package(ndt_map)
#find_package(ndt_registration)
if( NOT ndt_map_FOUND OR NOT ndt_registration_FOUND)
endif()
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)
## System dependencies are found with CMake's conventions
find_package(Boost REQUIRED COMPONENTS program_options system)
#find_package(PCL REQUIRED COMPONENTS common io)

include_directories(
  include
	${Boost_INCLUDE_DIRS}
	${PROJECT_SOURCE_DIR}/oru_minimal/include
	${PCL_INCLUDE_DIRS}
	${ndt_map_INCLUDE_DIRS}
	${ndt_registration_INCLUDE_DIRS}
	${EIGEN3_INCLUDE_DIR}
	${PCL_INCLUDE_DIRS})
link_directories ( ${Boost_LIBRARY_DIRS} )
include(pcl_find_sse.cmake)
PCL_CHECK_FOR_SSE()
#find_package(Boost REQUIRED)

set (${PROJECT_NAME}_LIB_SRCS
	oru_minimal/src/ndt_map.cpp
	oru_minimal/src/lazy_grid.cpp
	oru_minimal/src/ndt_matcher_d2d.cpp
	oru_minimal/src/ndt_matcher_d2d_2d.cpp
	oru_minimal/src/ndt_cell.cpp
	oru_minimal/src/cell_vector.cpp
	src/ndt_matcher_d2d_se.cpp
	src/ndt_matcher_p2d_se.cpp
	src/se_ndt.cpp
	)

add_library(${PROJECT_NAME} ${${PROJECT_NAME}_LIB_SRCS})
target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES})

add_library(ndt_fuser src/ndt_fuser_hmt_se.cpp)
target_link_libraries(ndt_fuser ${PROJECT_NAME})

add_executable(polesRegister_fuser src/polesRegister_fuser.cpp)
target_link_libraries(polesRegister_fuser ndt_fuser ${Boost_LIBRARIES}
	${OpenCV_LIBS} ${PCL_LIBRARIES})

add_executable(smoothness src/smoothness.cpp)
target_link_libraries(smoothness ${PROJECT_NAME} ${ndt_map_LIBRARIES}
	${ndt_registration_LIBRARIES} ${PCL_LIBRARIES} ${OpenCV_LIBS} ${Boost_LIBRARIES} -lboost_serialization)

add_executable(polesRegister src/polesRegister.cpp)
target_link_libraries(polesRegister ${PROJECT_NAME} ${ndt_map_LIBRARIES}
	${ndt_registration_LIBRARIES} ${PCL_LIBRARIES} ${OpenCV_LIBS} ${Boost_LIBRARIES} -lboost_serialization)

add_executable(polesRegister_ndt src/polesRegister_ndt.cpp)
target_link_libraries(polesRegister_ndt ${PROJECT_NAME} ${ndt_map_LIBRARIES}
	${ndt_registration_LIBRARIES} ${PCL_LIBRARIES} ${OpenCV_LIBS} ${Boost_LIBRARIES} -lboost_serialization)

add_executable(concat_transform src/concat_transform.cpp)
target_link_libraries(concat_transform ${PROJECT_NAME} ${ndt_map_LIBRARIES}
	${ndt_registration_LIBRARIES} ${PCL_LIBRARIES} ${OpenCV_LIBS} ${Boost_LIBRARIES} -lboost_serialization)

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
