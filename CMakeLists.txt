cmake_minimum_required(VERSION 3.20...3.28)

project(pidgeot)

# set(CMAKE_CXX_STANDARD 26)
# add_compile_options(-std=gnu++26)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=gnu++26")
set(CMAKE_CXX_STANDARD_REQUIRED On)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# function includes
include(${CMAKE_CURRENT_LIST_DIR}/cmake/setTargetOutputDirectories.cmake)

include(${CMAKE_CURRENT_LIST_DIR}/3rdparty/dependencies.cmake)

add_subdirectory(src)
