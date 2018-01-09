cmake_minimum_required(VERSION 3.8 FATAL_ERROR)

project(lgca LANGUAGES CXX)

# Instruct CMake to run moc and uic automatically when needed
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)

# Use C++ 11
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CUDA_STANDARD 11)

# Find packages
set(CMAKE_MODULE_PATH "/opt/vtk/VTK-8.1.0/CMake")
find_package(TBB REQUIRED)
find_package(OpenMP REQUIRED)

set(VTK_DIR "/opt/vtk/VTK-8.1.0/build")
find_package(VTK 8.0 REQUIRED)
include(${VTK_USE_FILE})

find_package(Qt5Widgets REQUIRED)

# Pass options to GCC
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -march=native")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

# Specify include directories
include_directories(
  ${CMAKE_CURRENT_BINARY_DIR}
  ${PROJECT_SOURCE_DIR}/src
  ${PROJECT_SOURCE_DIR}/lib
  ${TBB_INCLUDE_DIRS}
)

# Add library source files
file(GLOB LIB_LGCA_SOURCES  ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB LIB_LGCA_HEADERS  ${PROJECT_SOURCE_DIR}/src/*.h)
file(GLOB TCLAP_HEADERS     ${PROJECT_SOURCE_DIR}/lib/tclap/*.h)

# Set executable output path
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY
  ${PROJECT_SOURCE_DIR}/bin
)

add_subdirectory(${PROJECT_SOURCE_DIR}/apps/pipe)
add_subdirectory(${PROJECT_SOURCE_DIR}/apps/karman)
add_subdirectory(${PROJECT_SOURCE_DIR}/apps/diffusion)
#add_subdirectory(${PROJECT_SOURCE_DIR}/apps/box)
add_subdirectory(${PROJECT_SOURCE_DIR}/apps/single)
#add_subdirectory(${PROJECT_SOURCE_DIR}/apps/periodic)



