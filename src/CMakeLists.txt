project(pidgeot)

add_library(pidgeot measurement.cpp utils.cpp)
target_link_libraries(pidgeot PUBLIC Eigen3::Eigen CLI11::CLI11 pb_utils)
target_compile_features(pidgeot PUBLIC cxx_std_23)
set_target_output_directories(pidgeot)

add_executable(2d_rotation_sync main.cpp)
target_link_libraries(2d_rotation_sync PRIVATE pidgeot)
set_target_output_directories(2d_rotation_sync)
