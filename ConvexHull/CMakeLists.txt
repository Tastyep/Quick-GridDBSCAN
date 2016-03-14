cmake_minimum_required(VERSION 2.8)

project (ConvexHull)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc)
file(GLOB SOURCES "src/*.cpp")

add_library(ConvexHull STATIC ${SOURCES})