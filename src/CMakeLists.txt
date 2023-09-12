cmake_minimum_required(VERSION 3.10)

project (main VERSION 1.0)

set(CMAKE_TOOLCHAIN_FILE ${CMAKE_CURRENT_SOURCE_DIR}/vcpkg/scripts/buildsystems/vcpkg.cmake
  CACHE STRING "Vcpkg toolchain file")

find_package(SEAL)
find_package(OpenMP)

add_executable(main pdte.cpp server.cpp client.cpp utils.cpp tree_utils.cpp parser.cpp)
add_executable(cmp_bench cmp_bench.cpp server.cpp utils.cpp tree_utils.cpp parser.cpp)

target_link_libraries(main SEAL::seal OpenMP::OpenMP_CXX)
target_link_libraries(cmp_bench SEAL::seal OpenMP::OpenMP_CXX)

set(CMAKE_CXX_FLAGS "-g")
