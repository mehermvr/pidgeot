# taken from https://github.com/infiniflow/infinity/blob/main/CMakeLists.txt
# ensure ninja is being used
if(NOT
   CMAKE_GENERATOR
   STREQUAL
   "Ninja")
  message(
    FATAL_ERROR
      "This project requires the Ninja generator. Refers to https://cmake.org/cmake/help/latest/manual/cmake-cxxmodules.7.html#generator-support"
  )
endif()

# ensure clang is being used
execute_process(COMMAND ${CMAKE_CXX_COMPILER} --version OUTPUT_VARIABLE clang_full_version_string)
string(
  REGEX
  REPLACE ".*clang version ([0-9]+\\.[0-9]+).*"
          "\\1"
          CLANG_VERSION_STRING
          ${clang_full_version_string})
if(CLANG_VERSION_STRING VERSION_GREATER 16)
  # Print CMake version and project name
  message(STATUS "Building ${PROJECT_NAME} with CMake version: ${CMAKE_VERSION} On CLANG-${CLANG_VERSION_STRING}")

else()
  message(FATAL_ERROR "Please use clang version 17.0 and above")
endif()
