cmake_minimum_required(VERSION 3.0)
project(FIT)
set(CMAKE_CXX_STANDARD 11) 
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

aux_source_directory(${PROJECT_SOURCE_DIR}/src SRC)

add_executable(fit ${SRC})

message("EIGEN3_INCLUDE_DIR: ${EIGEN3_INCLUDE_DIR}")

target_include_directories(fit PUBLIC ${EIGEN3_INCLUDE_DIR})


