cmake_minimum_required(VERSION 3.20...3.25)

project(lsq_hw)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# function includes
include(cmake/setTargetOutputDirectories.cmake)

include(FetchContent)
FetchContent_Declare(
    cli11_proj
    QUIET
    GIT_REPOSITORY https://github.com/CLIUtils/CLI11.git
    GIT_TAG f4d0731
)

FetchContent_MakeAvailable(cli11_proj)

add_subdirectory(src)
