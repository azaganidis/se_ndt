cmake_minimum_required(VERSION 2.8.3)
project(se_ndt)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

##set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS "-Wno-dev -Wno-ignored-attributes -Wno-deprecated-declarations -Wno-narrowing")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

if(${CMAKE_VERSION} VERSION_GREATER "3.14.0")
    add_compile_definitions(N_THREADS=4)
else()
    add_definitions(-DN_THREADS=4)
endif()
find_package(Ceres REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++14")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
find_package(Boost REQUIRED COMPONENTS program_options system serialization
    filesystem iostreams)
link_directories ( ${Boost_LIBRARY_DIRS} ${CERES_LIBRARY_DIRS})
include(pcl_find_sse.cmake)
PCL_CHECK_FOR_SSE()

set (${PROJECT_NAME}_LIB_SRCS
	oru_minimal/src/ndt_map.cpp
	oru_minimal/src/lazy_grid.cpp
	oru_minimal/src/ndt_matcher_d2d.cpp
	oru_minimal/src/ndt_cell.cpp
	oru_minimal/src/cell_vector.cpp
	src/ndt_matcher_d2d_se.cpp
)

add_library(pose_optimizer src/pose_graph_3d.cc)
target_link_libraries(pose_optimizer ${EIGEN3_LIBRARIES} ${CERES_LIBRARIES})

add_library(ndt_libs ${${PROJECT_NAME}_LIB_SRCS})
target_link_libraries(ndt_libs ${Boost_LIBRARIES} ${PCL_COMMON_LIBRARIES} ${PCL_LIBRARIES})

add_library(other_fin src/other.cpp)
target_link_libraries(other_fin ndt_libs ${PCL_COMMON_LIBRARIES} ${PCL_LIBRARIES} )

add_library(histogram src/ndt_histogram.cpp)
target_link_libraries(histogram ndt_libs ${PCL_COMMON_LIBRARIES} ${PCL_LIBRARIES} )

add_library(${PROJECT_NAME} src/se_ndt.cpp)
target_link_libraries(${PROJECT_NAME} histogram other_fin pose_optimizer ndt_libs ${PCL_COMMON_LIBRARIES} ${PCL_LIBRARIES} )

option(WITH_GL "Enable/Disable OpenGL visualization of distributions" OFF)
if(WITH_GL)
    set(OpenGL_GL_PREFERENCE "GLVND")
    add_definitions(-DGL_VISUALIZE)
    find_package(GLUT REQUIRED)
    find_package(OpenGL REQUIRED)
    add_definitions(-DVISUALIZE)
    add_library(NDTVizGlut oru_minimal/src/NDTVizGlut.cc)
    target_link_libraries(NDTVizGlut ${GLUT_LIBRARY} ${OPENGL_LIBRARIES})
    list(APPEND LINK_LIBS "NDTVizGlut")
endif()

add_executable(opo src/offline_graph_opt.cpp)
target_link_libraries(opo pose_optimizer)

add_executable(kitti_slam src/kitti_slam.cpp)
target_link_libraries(kitti_slam ${PROJECT_NAME} other_fin
    ${PCL_FILTERS_LIBRARIES}
    ${LINK_LIBS}
	)

include_directories(
	include
	${Boost_INCLUDE_DIRS}
	${PROJECT_SOURCE_DIR}/oru_minimal/include
	${EIGEN3_INCLUDE_DIR}
	${PCL_INCLUDE_DIRS})
