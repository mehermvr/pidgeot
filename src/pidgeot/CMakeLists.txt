add_library(
  pidgeot
  measurement.cpp
  linearSystem.cpp
  dogLegSolver.cpp
  steepestDescentSolver.cpp
  gaussNewtonSolver.cpp)
target_link_libraries(pidgeot PUBLIC Eigen3::Eigen CLI11::CLI11 pb_utils)
target_compile_features(pidgeot PUBLIC cxx_std_23)
target_include_directories(pidgeot PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src>)
set_target_output_directories(pidgeot)
