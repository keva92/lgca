cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(lgca LANGUAGES CXX CUDA)

# Instruct CMake to run moc and uic automatically when needed
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)

# Use C++ 11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CUDA_STANDARD 11)

# Find packages
set(CUDA_TOOLKIT_ROOT_DIR "/usr/local/cuda")
find_package(CUDA 8.0 REQUIRED)

set(CMAKE_MODULE_PATH /opt/vtk/VTK-8.0.1/CMake)
find_package(TBB REQUIRED)
find_package(OpenMP REQUIRED)

set(VTK_DIR "/opt/vtk/VTK-8.0.1/build")
find_package(VTK 8.0 REQUIRED)
include(${VTK_USE_FILE})

find_package(Qt5Widgets REQUIRED)

# Pass options to GCC
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -march=native")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

# Pass options to NVCC
set(CMAKE_CUDA_FLAGS
 "${CMAKE_CUDA_FLAGS} -arch=sm_30"
)

# Specify include directories
include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_SOURCE_DIR}/lib
  ${CUDA_INCLUDE_DIRS}
  ${TBB_INCLUDE_DIRS}
)

# Compile and install CUDA PTX files
file(GLOB LIB_LGCA_CUDA_SOURCES ${PROJECT_SOURCE_DIR}/src/*.cu)
file(GLOB LIB_LGCA_CUDA_HEADERS ${PROJECT_SOURCE_DIR}/src/*.cuh)
add_library(CudaPTX OBJECT ${PROJECT_SOURCE_DIR}/src/cu_lattice.cu)
set_property(TARGET CudaPTX PROPERTY CUDA_PTX_COMPILATION ON)
install(TARGETS CudaPTX
   OBJECTS DESTINATION ${PROJECT_SOURCE_DIR}/ptx
)

# Add library source files
file(GLOB LIB_LGCA_SOURCES  ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB LIB_LGCA_HEADERS  ${PROJECT_SOURCE_DIR}/src/*.h)
file(GLOB LIB_LGCA_UIS      ${PROJECT_SOURCE_DIR}/src/*.ui)
file(GLOB TCLAP_HEADERS     ${PROJECT_SOURCE_DIR}/lib/tclap/*.h)

# Set executable output path
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
  ${PROJECT_SOURCE_DIR}/bin
)

#add_subdirectory(${PROJECT_SOURCE_DIR}/apps/pipe)
#add_subdirectory(${PROJECT_SOURCE_DIR}/apps/karman)
#add_subdirectory(${PROJECT_SOURCE_DIR}/apps/diffusion)
#add_subdirectory(${PROJECT_SOURCE_DIR}/apps/box)
add_subdirectory(${PROJECT_SOURCE_DIR}/apps/single)
#add_subdirectory(${PROJECT_SOURCE_DIR}/apps/periodic)



