project(2d_rotation_sync)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen CLI11::CLI11)

set_target_output_directories(${PROJECT_NAME})
