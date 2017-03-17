# Catkin User Guide: http://www.ros.org/doc/groovy/api/catkin/html/user_guide/user_guide.html
# Catkin CMake Standard: http://www.ros.org/doc/groovy/api/catkin/html/user_guide/standards.html
cmake_minimum_required(VERSION 2.8.3)
set(CMAKE_BUILD_TYPE Release) 
project(ndt_registration)
# Load catkin and all dependencies required for this package
# TODO: remove all from COMPONENTS that are not catkin packages.
find_package(catkin REQUIRED COMPONENTS ndt_map tf tf_conversions pcl_ros)
catkin_package(
    CATKIN_DEPENDS ndt_map tf tf_conversions pcl_ros
    INCLUDE_DIRS include
    LIBRARIES ${PROJECT_NAME} 
)

find_package(Boost REQUIRED COMPONENTS program_options system)
include_directories(include ${Boost_INCLUDE_DIR} ${catkin_INCLUDE_DIRS})

#uncomment if you have defined messages
#add_message_files(
#  FILES
  # TODO: List your msg files here
#)
#uncomment if you have defined services
#add_service_files(
#  FILES
  # TODO: List your msg files here
#)

#add_executable(ndt_registration_test test/ndt_registration_test.cc)
#rosbuild_add_openmp_flags(ndt_registration_test)

find_package(OpenMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -g -O0")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -g -O0")

set (${PROJECT_NAME}_LIB_SRCS
	src/ndt_matcher_d2d.cpp
	src/ndt_matcher_d2d_2d.cpp	
	src/ndt_matcher_d2d_feature.cpp  
	src/ndt_matcher_p2d.cpp
    )
add_library(${PROJECT_NAME} ${${PROJECT_NAME}_LIB_SRCS})
target_link_libraries(${PROJECT_NAME} ${catkin_LIBRARIES})

add_executable(simple_test test/simple.cc)
target_link_libraries(simple_test ${PROJECT_NAME} ${Boost_LIBRARIES} ${catkin_LIBRARIES})

add_executable(d2d_deathmatch test/deathmatch.cc)
target_link_libraries(d2d_deathmatch ${PROJECT_NAME} ${Boost_LIBRARIES} ${catkin_LIBRARIES})

add_executable(p2d_deathmatch test/p2d_deathmatch.cc)
target_link_libraries(p2d_deathmatch ${PROJECT_NAME} ${Boost_LIBRARIES} ${catkin_LIBRARIES})

#add_executable(convert_logs test/convert.cc)
#target_link_libraries(convert_logs ${Boost_LIBRARIES} ${catkin_LIBRARIES})

#add_executable(from_jff test/from_jff.cc)
#target_link_libraries(from_jff ${Boost_LIBRARIES} ${catkin_LIBRARIES})
#rosbuild_add_openmp_flags(from_jff)

#add_executable(ndt_feature_test test/ndt_feature_test.cc)
#target_link_libraries(ndt_feature_test ${Boost_LIBRARIES} ${catkin_LIBRARIES})

#add_executable(feature_registration test/feature_registration.cc)
#target_link_libraries(feature_registration ndt_map pointcloud_vrml )

#add_executable(ndt_reconstruct test/reconstruct.cc)
#target_link_libraries(ndt_reconstruct ${PROJECT_NAME} ndt_map pointcloud_vrml )

#add_executable(ndt_matcher_node test/ndt_matcher_node.cc)
#target_link_libraries(ndt_matcher_node ${PROJECT_NAME} ndt_map pointcloud_vrml )

#add_executable(batchTestRangerAccuracy test/batchTestRangerAccuracy.cc)
#target_link_libraries(batchTestRangerAccuracy ${PROJECT_NAME} ndt_map  pointcloud_vrml)

#add_executable(narf_keypoint_extraction test/narf_extraction.cc)
#target_link_libraries(narf_keypoint_extraction ${PROJECT_NAME} ndt_map pointcloud_vrml )

#add_executable(narf_registration test/narf_registration.cc)
#target_link_libraries(narf_registration ${PROJECT_NAME} ndt_map pointcloud_vrml )

#add_executable(ndt_histogram_registration test/ndt_histogram_registration.cc)
#target_link_libraries(ndt_histogram_registration ${PROJECT_NAME} ndt_map pointcloud_vrml )

#add_executable(test_covariance test/test_covariance.cc)
#target_link_libraries(test_covariance ${PROJECT_NAME} ndt_map pointcloud_vrml )

#target_link_libraries(test_ndt_mapper ${PROJECT_NAME} ndt_map pointcloud_vrml )

add_executable(check_reg test/check_reg.cc)
target_link_libraries(check_reg
  ${PROJECT_NAME}
  #${Boost_LIBRARIES}
  ${catkin_LIBRARIES}
  )

## Generate added messages and services with any dependencies listed here
#generate_messages(
    #TODO DEPENDENCIES geometry_msgs std_msgs
#)

# catkin_package parameters: http://ros.org/doc/groovy/api/catkin/html/dev_guide/generated_cmake_api.html#catkin-package
# TODO: fill in what other packages will need to use this package
install(TARGETS simple_test ${PROJECT_NAME} 
        ARCHIVE DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
        LIBRARY DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION}
        RUNTIME DESTINATION ${CATKIN_PACKAGE_BIN_DESTINATION}
       )

install(DIRECTORY include/${PROJECT_NAME}/
	DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION}
       )
