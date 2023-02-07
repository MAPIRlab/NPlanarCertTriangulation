project(NPlanarViewsTrian-Tests CXX)

# Find Eigen library
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../cmake/")
find_package(Eigen3 3.3.3 REQUIRED)
if(EIGEN3_FOUND)
message(STATUS "Found Eigen3 library (version ${EIGEN3_VERSION})")
message(STATUS "Eigen3 include directory: ${EIGEN3_INCLUDE_DIR}\n")
else()
message(STATUS "Eigen library not found!")
endif()

# Add the Eigen include directories
include_directories(${EIGEN3_INCLUDE_DIR})



add_library(exp_help_test experimentsHelper.h ../utils/generatePointCloud.h experimentsHelper.cpp ../utils/generatePointCloud.cpp)


## Generic test synthetic
add_executable(generic_test ${CMAKE_CURRENT_SOURCE_DIR}/generic_test.cpp)
target_link_libraries(generic_test 
                        NPlanarTrian 
                        exp_help_test   
                        )
        

