cmake_minimum_required(VERSION 2.8)

project (QGDBSAN)
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -Wall -g")
SET(CMAKE_EXE_LINKER_FLAGS  "-g")

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/inc)
file(GLOB SOURCES "src/*.cpp")

add_library(QGDBSCAN STATIC ${SOURCES})