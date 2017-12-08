cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(lgca LANGUAGES CXX CUDA)

# Use C++ 11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CUDA_STANDARD 11)

# Find packages
set(CUDA_TOOLKIT_ROOT_DIR "/usr/local/cuda")
find_package(CUDA 8.0 REQUIRED)

find_package(OpenMP REQUIRED)

set(CMAKE_MODULE_PATH /opt/vtk/VTK-8.0.1/CMake)
find_package(TBB REQUIRED)

find_package(PNG REQUIRED)
find_package(Freetype REQUIRED)

set(VTK_DIR "/opt/vtk/VTK-8.0.1/build")
find_package(VTK 8.0 REQUIRED)
include(${VTK_USE_FILE})

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
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_SOURCE_DIR}/lib
  ${CUDA_INCLUDE_DIRS}
  ${TBB_INCLUDE_DIRS}
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
add_executable(
  ${PROJECT_NAME}
  ${LGCA_SOURCES} ${LGCA_HEADERS}
  ${PNGWRITER_SOURCES} ${PNGWRITER_HEADERS}
  ${TCLAP_HEADERS}
)

# Specify target and libraries to link it with
target_link_libraries(
  ${PROJECT_NAME}
  ${CUDA_LIBRARIES}
  ${OpenMP_C_LIBRARIES} ${OpenMP_CXX_LIBRARIES}
  ${TBB_LIBRARIES}
  ${PNG_LIBRARIES}
  ${FREETYPE_LIBRARIES}
  ${VTK_LIBRARIES}
)
