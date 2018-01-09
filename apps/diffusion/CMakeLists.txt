# Add application source files
file(GLOB LGCA_DIFFUSION_SOURCES *.cpp)
file(GLOB LGCA_DIFFUSION_HEADERS *.h)
file(GLOB LGCA_DIFFUSION_UIS *.ui)

# Specify target and source files to compile from
add_executable(
  ${PROJECT_NAME}-diffusion
  ${LGCA_DIFFUSION_SOURCES} ${LGCA_DIFFUSION_HEADERS} ${LGCA_DIFFUSION_UIS}
  ${LIB_LGCA_SOURCES} ${LIB_LGCA_HEADERS}
  ${TCLAP_HEADERS}
)

# Specify target and libraries to link with
target_link_libraries(
  ${PROJECT_NAME}-diffusion
  ${CUDA_LIBRARIES}
  ${OpenMP_C_LIBRARIES} ${OpenMP_CXX_LIBRARIES}
  ${TBB_LIBRARIES}
  ${VTK_LIBRARIES}
  Qt5::Widgets
)
