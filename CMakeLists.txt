cmake_minimum_required(VERSION 3.29)

project(pidgeot_main)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# function includes
include(${CMAKE_CURRENT_LIST_DIR}/cmake/setTargetOutputDirectories.cmake)

include(${CMAKE_CURRENT_LIST_DIR}/3rdparty/dependencies.cmake)

add_subdirectory(src)
