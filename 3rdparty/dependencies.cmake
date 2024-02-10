# should handle find_package or fetch content make availabe all necessary
# dependencies. see kiss icp for reference. but i want my own stuff.

# in the end it ended up pretty much the same as the kiss icp one. surprise.
function(find_dep PACKAGE_NAME TARGET_NAME FETCH_CMAKE_FILE)
  string(TOUPPER ${PACKAGE_NAME} PACKAGE_NAME_UP)
  set(USE_FROM_SYSTEM_OPTION "USE_SYSTEM_${PACKAGE_NAME_UP}")
  if(${${USE_FROM_SYSTEM_OPTION}})
    message(STATUS "Searching for ${PACKAGE_NAME} in system. ")
    find_package(${PACKAGE_NAME} QUIET)
  endif()
  if(NOT TARGET ${TARGET_NAME})
    message(STATUS "${PACKAGE_NAME}'s target: ${TARGET_NAME} not found.")
    message(STATUS "Fetching it using FetchContent (${FETCH_CMAKE_FILE})")
    include(${FETCH_CMAKE_FILE})
  endif()
endfunction()

option(USE_SYSTEM_EIGEN3 "Use system Eigen" OFF)
find_dep("Eigen3" "Eigen3::Eigen" "${CMAKE_CURRENT_LIST_DIR}/eigen/eigen.cmake")
option(USE_SYSTEM_CLI11 "Use system CLI11" OFF)
find_dep("CLI11" "CLI11::CLI11" "${CMAKE_CURRENT_LIST_DIR}/cli11/cli11.cmake")
find_dep("pb_utils" "pb_utils"
         "${CMAKE_CURRENT_LIST_DIR}/pb_utils/pb_utils.cmake")
