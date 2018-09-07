cmake_minimum_required(VERSION 2.8.3)
project(se_ndt)
if(NOT CMAKE_BUILD_TYPE)
      set(CMAKE_BUILD_TYPE Release)
  endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(OpenMP REQUIRED)
#find_package(OpenCV REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
find_package(Boost REQUIRED COMPONENTS program_options system)

find_package(roscpp REQUIRED)
find_package(rosconsole REQUIRED)
find_package(tf2_ros REQUIRED)
find_package(pcl_ros REQUIRED)
find_package(VTK REQUIRED)
find_package(GLUT REQUIRED)
find_package(OpenGL REQUIRED)
set(ROS_REQUIRED_LIBS
	${roscpp_LIBRARIES}
	${rosconsole_LIBRARIES}
	${tf2_ros_LIBRARIES}
	${pcl_ros_LIBRARIES}
    ${VTK_LIBRARIES}
	)

set(ROS_REQUIRED_DIRS
    ${roscpp_INCLUDE_DIRS}
    ${rosconsole_INCLUDE_DIRS}
    ${tf2_ros_INCLUDE_DIRS}
    ${pcl_ros_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIRS}
	)

include_directories(
  include
	${Boost_INCLUDE_DIRS}
	${PROJECT_SOURCE_DIR}/oru_minimal/include
	${PCL_INCLUDE_DIRS}
	${ndt_map_INCLUDE_DIRS}
	${ndt_registration_INCLUDE_DIRS}
	${EIGEN3_INCLUDE_DIR}
    ${ROS_REQUIRED_DIRS}
	${PCL_INCLUDE_DIRS})
link_directories ( ${Boost_LIBRARY_DIRS} ${G2O_LIBRARY_DIRS})
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
add_library(NDTVizGlut oru_minimal/src/NDTVizGlut.cc)
target_link_libraries(NDTVizGlut ${GLUT_LIBRARY} ${OPENGL_LIBRARIES})

add_library(${PROJECT_NAME} ${${PROJECT_NAME}_LIB_SRCS})
target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES})

add_library(ndt_fuser src/ndt_fuser_hmt_se.cpp)
target_link_libraries(ndt_fuser ${PROJECT_NAME})

add_executable(semanticRegister src/semanticRegister.cpp)
target_link_libraries(semanticRegister ${PROJECT_NAME} ${ndt_map_LIBRARIES}
    ndt_fuser ${ndt_registration_LIBRARIES} ${PCL_COMMON_LIBRARIES}
    ${PCL_FILTERS_LIBRARIES}
    NDTVizGlut
	${Boost_LIBRARIES} -lboost_serialization)

add_executable(ros_register src/ros_register.cpp)
target_link_libraries(ros_register ${PROJECT_NAME} ${ndt_map_LIBRARIES}
    ndt_fuser ${ndt_registration_LIBRARIES} ${PCL_COMMON_LIBRARIES}
    ${PCL_FILTERS_LIBRARIES}
    NDTVizGlut
    ${ROS_REQUIRED_LIBS}
	${Boost_LIBRARIES} -lboost_serialization)

