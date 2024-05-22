cmake_minimum_required (VERSION 3.1)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_STANDARD 11)
project(cmdline_LBM)
include_directories("${CMAKE_SOURCE_DIR}/LBM" "${CMAKE_SOURCE_DIR}/LBM/headers" "${CMAKE_SOURCE_DIR}/cmdline/rapidjson" "${CMAKE_SOURCE_DIR}/cmdline/fast-cpp-csv-parser")
find_package (Threads)
add_executable(main.o ${CMAKE_SOURCE_DIR}/LBM/source/LBM.cpp ${CMAKE_SOURCE_DIR}/LBM/source/vector3.cpp ${CMAKE_SOURCE_DIR}/cmdline/main.cpp)
target_sources(main.o
PRIVATE
  ${CMAKE_SOURCE_DIR}/LBM/source/LBM.cpp
  ${CMAKE_SOURCE_DIR}/LBM/source/vector3.cpp
  ${CMAKE_SOURCE_DIR}/cmdline/rapidjson/document.h
  ${CMAKE_SOURCE_DIR}/cmdline/rapidjson/writer.h
  ${CMAKE_SOURCE_DIR}/cmdline/rapidjson/stringbuffer.h
  ${CMAKE_SOURCE_DIR}/cmdline/fast-cpp-csv-parser/csv.h
)
target_link_libraries(main.o ${CMAKE_THREAD_LIBS_INIT})
