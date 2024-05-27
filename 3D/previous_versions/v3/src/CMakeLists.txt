add_compile_options(-O3)

cmake_minimum_required (VERSION 3.1)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_STANDARD 11)
project(cmdline_LBM)
include_directories("${CMAKE_SOURCE_DIR}/LBM" "${CMAKE_SOURCE_DIR}/LBM/headers" "${CMAKE_SOURCE_DIR}/cmdline/rapidjson" "${CMAKE_SOURCE_DIR}/cmdline/fast-cpp-csv-parser" "${CMAKE_SOURCE_DIR}/LBM/include")
find_package (Threads)
add_executable(main.o ${CMAKE_SOURCE_DIR}/cmdline/main.cpp ${CMAKE_SOURCE_DIR}/LBM/include/utils.cpp ${CMAKE_SOURCE_DIR}/LBM/include/profile.cpp)
target_sources(main.o
PRIVATE
  ${CMAKE_SOURCE_DIR}/cmdline/csv.h
  ${CMAKE_SOURCE_DIR}/LBM/headers/utils.h
  ${CMAKE_SOURCE_DIR}/LBM/headers/profile.h
)
target_link_libraries(main.o ${CMAKE_THREAD_LIBS_INIT})
