project(testsrc)

find_package(Eigen3 REQUIRED NO_MODULE)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen)

set_target_output_directories(${PROJECT_NAME})
