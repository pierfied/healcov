cmake_minimum_required(VERSION 3.15)
project(healcov)

set(CMAKE_CXX_STANDARD 17)

set(HEALPIX_DIR path/to/healpix)

add_library(healcov SHARED src/healcov.h src/healcov.cpp)

target_include_directories(healcov PUBLIC ${HEALPIX_DIR}/include)

target_link_libraries(healcov ${HEALPIX_DIR}/lib/libhealpix_cxx.a)