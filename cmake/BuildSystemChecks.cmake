# taken from https://github.com/infiniflow/infinity/blob/main/CMakeLists.txt
# ensure ninja is being used
if(NOT
   CMAKE_GENERATOR
   STREQUAL
   "Ninja")
  message(
    # FATAL_ERROR
    AUTHOR_WARNING
      "This project requires the Ninja generator. Refers to https://cmake.org/cmake/help/latest/manual/cmake-cxxmodules.7.html#generator-support"
  )
endif()

# ensure ~~clang~~ gcc-14 is being used
# execute_process(COMMAND ${CMAKE_CXX_COMPILER} --version OUTPUT_VARIABLE clang_full_version_string)
execute_process(COMMAND ${CMAKE_CXX_COMPILER} --version OUTPUT_VARIABLE gcc_full_version_string)
#   REPLACE ".*clang version ([0-9]+\\.[0-9]+).*"
#           "\\1"
string(
  REGEX
  REPLACE ".*GCC\\) ([0-9]+\\.[0-9]+).*"
          "\\1"
          GCC_VERSION_STRING
          ${gcc_full_version_string})
if(GCC_VERSION_STRING VERSION_GREATER_EQUAL 14)
  # Print CMake version and project name
  message(STATUS "Building ${PROJECT_NAME} with CMake version: ${CMAKE_VERSION} On GCC-${GCC_VERSION_STRING}")

else()
  message(FATAL_ERROR "Please use gcc version 14.0 and above")
endif()
