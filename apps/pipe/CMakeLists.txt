# Add application source files
file(GLOB LGCA_PIPE_SOURCES *.cpp)
file(GLOB LGCA_PIPE_HEADERS *.h)
file(GLOB LGCA_PIPE_UIS *.ui)

# Specify target and source files to compile from
add_executable(
  ${PROJECT_NAME}-pipe
  ${LGCA_PIPE_SOURCES} ${LGCA_PIPE_HEADERS} ${LGCA_PIPE_UIS}
  ${LIB_LGCA_SOURCES} ${LIB_LGCA_HEADERS}
  ${TCLAP_HEADERS}
)

# Specify target and libraries to link with
target_link_libraries(
  ${PROJECT_NAME}-pipe
  ${CUDA_LIBRARIES}
  ${OpenMP_C_LIBRARIES} ${OpenMP_CXX_LIBRARIES}
  ${TBB_LIBRARIES}
  ${VTK_LIBRARIES}
  Qt5::Widgets
)
