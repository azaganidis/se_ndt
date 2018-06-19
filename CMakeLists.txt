cmake_minimum_required(VERSION 2.8.3)
project(se_ndt)
include(GNUInstallDirs)
if(NOT CMAKE_BUILD_TYPE)
      set(CMAKE_BUILD_TYPE Release)
  endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(PCL REQUIRED COMPONENTS kdtree)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -std=c++11")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
find_package(Boost REQUIRED COMPONENTS program_options system)

include_directories(
  include
	${Boost_INCLUDE_DIRS}
	${PROJECT_SOURCE_DIR}/oru_minimal/include
	${EIGEN3_INCLUDE_DIR}
	${PCL_INCLUDE_DIRS})
link_directories ( ${Boost_LIBRARY_DIRS})
include(pcl_find_sse.cmake)
PCL_CHECK_FOR_SSE()

set(NDT_FILES oru_minimal/src/ndt_map.cpp
    oru_minimal/src/cell_vector.cpp
    oru_minimal/src/ndt_cell.cpp
    oru_minimal/src/ndt_matcher_d2d.cpp
    oru_minimal/src/lazy_grid.cpp)
add_library(ndt_map SHARED ${NDT_FILES})
# target_include_directories(ndt_map PUBLIC
#    "oru_minimal/include")
target_link_libraries(ndt_map ${PCL_KDTREE_LIBRARIES} )# ${PCL_LIBRARIES} )
install(DIRECTORY oru_minimal/include/ndt_map DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})
install(DIRECTORY oru_minimal/include/ndt_registration DESTINATION ${CMAKE_INSTALL_INCLUDEDIR})

add_library(${PROJECT_NAME} SHARED src/se_ndt.cpp src/ndt_matcher_d2d_se.cpp)
target_link_libraries(${PROJECT_NAME} ndt_map)

#target_include_directories(${PROJECT_NAME} PUBLIC
#    "${CMAKE_INSTALL_INCLUDEDIR}/include/se_ndt")
set_target_properties(${PROJECT_NAME} PROPERTIES PUBLIC_HEADER
    "include/se_ndt/ndt_matcher_d2d_se.h;include/se_ndt/se_ndt.hpp")
INSTALL(TARGETS ${PROJECT_NAME} ndt_map EXPORT ${PROJECT_NAME} 
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_ARCHIVEDIR}
    PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/${PROJECT_NAME})
list(APPEND se_ndt_LIBRARIES lib${PROJECT_NAME}.so)
list(APPEND se_ndt_LIBRARIES libndt_map.so)
configure_file(se_ndt-config.cmake se_ndtConfig.cmake @ONLY)
install(FILES ${CMAKE_BINARY_DIR}/se_ndtConfig.cmake DESTINATION
    "${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME}")
install(EXPORT ${PROJECT_NAME} DESTINATION
   ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME} 
   )
#install(FILES se_ndt-config.cmake DESTINATION
#   ${CMAKE_INSTALL_LIBDIR}/${PROJECT_NAME})

export(PACKAGE ${PROJECT_NAME})
