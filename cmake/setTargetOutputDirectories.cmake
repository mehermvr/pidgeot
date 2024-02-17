# mostly to use to set the output binary directory to bin inside the top level
# of the build tree instead of somewhere in the source
function(set_target_output_directories target_name)
  set_target_properties(
    ${target_name}
    PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${PROJECT_NAME}/lib"
               LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${PROJECT_NAME}/lib"
               RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/${PROJECT_NAME}/bin")
endfunction()
