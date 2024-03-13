cmake_minimum_required(VERSION 3.28)

project(pidgeot LANGUAGES CXX)
set(CMAKE_CXX_EXTENSIONS OFF)
include(${CMAKE_CURRENT_LIST_DIR}/cmake/BuildSystemChecks.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/cmake/setTargetOutputDirectories.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/cmake/StandardProjectOptions.cmake)

include(${CMAKE_CURRENT_LIST_DIR}/3rdparty/Dependencies.cmake)

add_subdirectory(src)
