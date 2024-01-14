cmake_minimum_required(VERSION 3.0)
project(SIMULATOR)
set(CMAKE_CXX_STANDARD 11) 
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC)

#设置头文件目录
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

add_executable(simulate ${SRC})

message("EIGEN3_INCLUDE_DIR: ${EIGEN3_INCLUDE_DIR}")

target_include_directories(simulate PUBLIC ${EIGEN3_INCLUDE_DIR})


