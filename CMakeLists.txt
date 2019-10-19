cmake_minimum_required(VERSION 3.15)
project(healcov)

set(CMAKE_CXX_STANDARD 17)

set(HEALPIX_DIR path/to/healpix)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

add_library(healcov SHARED src/healcov.h src/healcov.cpp)

target_include_directories(healcov PUBLIC ${HEALPIX_DIR}/include)

target_link_libraries(healcov ${HEALPIX_DIR}/lib/libhealpix_cxx.a)