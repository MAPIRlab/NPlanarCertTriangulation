project(NPlanarViewsTrian-Examples CXX)


# Find Eigen library
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/")
find_package(Eigen3 3.3.3 REQUIRED)
if(EIGEN3_FOUND)
message(STATUS "Found Eigen3 library (version ${EIGEN3_VERSION})")
message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}\n")
else()
message(STATUS "Eigen library not found!")
endif()


## For the iterative certifier


# Add the Eigen include directories
include_directories(${EIGEN3_INCLUDE_DIR})



add_library(exp_gen ../utils/generatePointCloud.h ../utils/generatePointCloud.cpp)




# Basic example
add_executable(example_base ${CMAKE_CURRENT_SOURCE_DIR}/example_base.cpp)
target_link_libraries(example_base 
                              NPlanarTrian
                              exp_gen          
                              )
                              

