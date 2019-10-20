cmake_minimum_required(VERSION 3.13)
project(healcov)

set(CMAKE_CXX_STANDARD 17)

find_package(OpenMP)
if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif ()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

add_executable(healcov src/healcov.h src/healcov.cpp)

target_include_directories(healcov PUBLIC $ENV{HEALPIX_DIR}/include)

target_link_libraries(healcov $ENV{HEALPIX_DIR}/lib/libhealpix_cxx.a $ENV{HEALPIX_DIR}/lib/libc_utils.a
        $ENV{HEALPIX_DIR}/lib/libcxxsupport.a $ENV{HEALPIX_DIR}/lib/libpocketfft.a $ENV{HEALPIX_DIR}/lib/libsharp.a)