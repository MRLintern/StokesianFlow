cmake_minimum_required(VERSION 3.18.3)
project(viscousRiver)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_REQUIRED True)

set(SRC_FILES src/main.cpp)

include_directories(src)

add_executable(main ${SRC_FILES})
target_link_libraries(main Eigen3::Eigen)
