add_compile_options(-march=native)

option(MNx "Define MNx macro")
option(MNy "Define MNy macro")
option(MNz "Define MNy macro")
option(MNt "Define MNt macro")

# Print to console
message("MNx: ${MNx}")
message("MNy: ${MNy}")
message("MNz: ${MNz}")
message("MNt: ${MNt}")

# Parameters have to be passed in input to the cmake command
add_compile_definitions(MNx=${MNx} MNy=${MNy} MNz=${MNz} MNt=${MNt})


cmake_minimum_required (VERSION 3.1)
set(CMAKE_CXX_COMPILER g++)
set(CMAKE_CXX_STANDARD 11)
project(cmdline_LBM)
include_directories("${CMAKE_SOURCE_DIR}/LBM" "${CMAKE_SOURCE_DIR}/LBM/headers" "${CMAKE_SOURCE_DIR}/cmdline/rapidjson" "${CMAKE_SOURCE_DIR}/cmdline/fast-cpp-csv-parser" "${CMAKE_SOURCE_DIR}/LBM/include")
find_package (Threads)
add_executable(main_${MNx}_${MNy}_${MNz}_${MNt}.o ${CMAKE_SOURCE_DIR}/cmdline/main.cpp ${CMAKE_SOURCE_DIR}/LBM/include/utils.cpp ${CMAKE_SOURCE_DIR}/LBM/include/profile.cpp)
target_sources(main_${MNx}_${MNy}_${MNz}_${MNt}.o
PRIVATE
  ${CMAKE_SOURCE_DIR}/cmdline/csv.h
  ${CMAKE_SOURCE_DIR}/LBM/headers/utils.h
  ${CMAKE_SOURCE_DIR}/LBM/headers/profile.h
)
target_compile_definitions(main_${MNx}_${MNy}_${MNz}_${MNt}.o PRIVATE OUTPUT)
target_link_libraries(main_${MNx}_${MNy}_${MNz}_${MNt}.o ${CMAKE_THREAD_LIBS_INIT})

#add_compile_definitions(PROFILE)
add_executable(main_${MNx}_${MNy}_${MNz}_${MNt}_profiler.o ${CMAKE_SOURCE_DIR}/cmdline/main.cpp ${CMAKE_SOURCE_DIR}/LBM/include/utils.cpp ${CMAKE_SOURCE_DIR}/LBM/include/profile.cpp)
target_compile_definitions(main_${MNx}_${MNy}_${MNz}_${MNt}_profiler.o PRIVATE PROFILE)
target_sources(main_${MNx}_${MNy}_${MNz}_${MNt}_profiler.o
PRIVATE
  ${CMAKE_SOURCE_DIR}/cmdline/csv.h
  ${CMAKE_SOURCE_DIR}/LBM/headers/utils.h
  ${CMAKE_SOURCE_DIR}/LBM/headers/profile.h
)
target_link_libraries(main_${MNx}_${MNy}_${MNz}_${MNt}_profiler.o ${CMAKE_THREAD_LIBS_INIT})

find_library(PAPI papi)
if (PAPI)
  add_executable(main_${MNx}_${MNy}_${MNz}_${MNt}_benchmark.o ${CMAKE_SOURCE_DIR}/cmdline/main.cpp ${CMAKE_SOURCE_DIR}/LBM/include/utils.cpp ${CMAKE_SOURCE_DIR}/LBM/include/profile.cpp)
  target_compile_definitions(main_${MNx}_${MNy}_${MNz}_${MNt}_benchmark.o PRIVATE BENCHMARK)
  target_link_libraries(main_${MNx}_${MNy}_${MNz}_${MNt}_benchmark.o PRIVATE ${PAPI})
  target_sources(main_${MNx}_${MNy}_${MNz}_${MNt}_benchmark.o
  PRIVATE
    ${CMAKE_SOURCE_DIR}/cmdline/csv.h
    ${CMAKE_SOURCE_DIR}/LBM/headers/utils.h
    ${CMAKE_SOURCE_DIR}/LBM/headers/profile.h
  )
  target_link_libraries(main_${MNx}_${MNy}_${MNz}_${MNt}_benchmark.o ${CMAKE_THREAD_LIBS_INIT})
endif()