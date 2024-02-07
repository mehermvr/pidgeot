project(happly)

add_executable(${PROJECT_NAME} test.cpp)
# target_link_libraries(${PROJECT_NAME} PUBLIC Eigen3::Eigen )

add_custom_target(
  link_data ALL
  COMMAND ${CMAKE_COMMAND} -E create_symlink ${PROJECT_SOURCE_DIR}/test_data.ply
          ${BIN_DIR}/test_data.ply
  COMMENT
    "Symlinking ${PROJECT_SOURCE_DIR}/test_data.ply to ${BIN_DIR}/test_data.ply"
)
set_target_properties(${PROJECT_NAME} PROPERTIES RUNTIME_OUTPUT_DIRECTORY
                                                 ${BIN_DIR})
