cmake_minimum_required(VERSION 3.16)

project(FourierHHG)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_PREFIX_PATH "~/usr/local")
include_directories(
    $ENV{HOME}/usr/local/include
)

find_package(OpenMP REQUIRED)
#find_package(MPI REQUIRED)
find_package(mrock REQUIRED)

if(CLUSTER_BUILD)
    include(${CMAKE_SOURCE_DIR}/cmake/ClusterCompilerFlags.cmake)
else()
    include(${CMAKE_SOURCE_DIR}/cmake/DefaultCompilerFlags.cmake)
endif()
include(${CMAKE_CURRENT_SOURCE_DIR}/../../PhdUtility/cmake/mrock-information.cmake)

add_subdirectory(sources/FHHG/Laser)

add_executable(fourier_hhg 
    sources/fourier_hhg.cpp
)
target_sources(fourier_hhg PRIVATE
    sources/FHHG/DiracSystem.cpp
    sources/FHHG/State.cpp
)

SET_COMPILER_FLAGS(fourier_hhg)

target_link_libraries(fourier_hhg PRIVATE Laser OpenMP::OpenMP_CXX ${mrock_LIBRARIES})
target_include_directories(fourier_hhg PRIVATE ${mrock_INCLUDE_DIRS})