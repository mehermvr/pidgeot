cmake_minimum_required(VERSION 3.20...3.25)

project(testground)
set(CMAKE_CXX_STANDARD 20)

# custom vars for later
set(BIN_DIR
    "${CMAKE_CURRENT_BINARY_DIR}/bin"
    CACHE STRING "Bin Directory.")

add_subdirectory(src/happly)
