# should handle find_package or fetch content make availabe all necessary dependencies.
# see kiss icp for reference. but i want my own stuff.
find_package(Eigen3 REQUIRED NO_MODULE)
include(${CMAKE_CURRENT_LIST_DIR}/cli11/cli11.cmake)
