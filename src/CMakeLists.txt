project(pidgeot)

add_library(${PROJECT_NAME} measurement.cpp utils.cpp)
target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen CLI11::CLI11 pb_utils)
target_compile_options(${PROJECT_NAME} PUBLIC -std=gnu++26)
set_target_output_directories(${PROJECT_NAME})

add_executable(2d_rotation_sync main.cpp)
target_link_libraries(2d_rotation_sync PRIVATE ${PROJECT_NAME})
set_target_output_directories(2d_rotation_sync)
