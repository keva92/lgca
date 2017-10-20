set(CMAKE_CXX_COMPILER "/opt/rh/devtoolset-4/root/bin/g++")
set(CMAKE_C_COMPILER "/opt/rh/devtoolset-4/root/bin/gcc")
set(CMAKE_CUDA_COMPILER "/usr/local/cuda/bin/nvcc")

project(lgca CXX CUDA)

cmake_minimum_required(VERSION 2.8)

# Use C++ 11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CUDA_STANDARD 11)

# Find packages
set(CUDA_TOOLKIT_ROOT_DIR "/usr/local/cuda")
find_package(CUDA 8.0 REQUIRED)

find_package(OpenMP REQUIRED)
find_package(PNG REQUIRED)
find_package(Freetype REQUIRED)

set(VTK_DIR "$ENV{HOME}/VTK-8.0.0/build" CACHE PATH "VTK directory override" FORCE)
find_package(VTK 8.0 REQUIRED)
include(${VTK_USE_FILE})

# Pass options to GCC
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -O3 -march=native")

# Pass options to NVCC
set(CUDA_NVCC_FLAGS
  ${CUDA_NVCC_FLAGS} -O3 -gencode arch=compute_61,code=sm_61
)

# Specify include directories
include_directories(
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_SOURCE_DIR}/lib
  ${PNG_INCLUDE_DIRS}
  ${FREETYPE_INCLUDE_DIRS}
)

# Add application source files
file(GLOB LGCA_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cpp ${PROJECT_SOURCE_DIR}/src/*.cu)
file(GLOB LGCA_HEADERS  ${PROJECT_SOURCE_DIR}/src/*.h ${PROJECT_SOURCE_DIR}/src/*.cuh)
file(GLOB PNGWRITER_SOURCES ${PROJECT_SOURCE_DIR}/lib/pngwriter/*.cc)
file(GLOB PNGWRITER_HEADERS ${PROJECT_SOURCE_DIR}/lib/pngwriter/*.h)
file(GLOB TCLAP_HEADERS ${PROJECT_SOURCE_DIR}/lib/tclap/*.h)

# Set executable output path
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
  ${PROJECT_SOURCE_DIR}/bin
)

# Specify target and source files to compile it from
cuda_add_executable(
  ${PROJECT_NAME}
  ${LGCA_SOURCES} ${LGCA_HEADERS}
  ${PNGWRITER_SOURCES} ${PNGWRITER_HEADERS}
  ${TCLAP_HEADERS}
)

# Specify target and libraries to link it with
target_link_libraries(
  ${PROJECT_NAME}
  ${OpenMP_C_LIBRARIES} ${OpenMP_CXX_LIBRARIES}
  ${PNG_LIBRARIES}
  ${FREETYPE_LIBRARIES}
  ${VTK_LIBRARIES}
)
