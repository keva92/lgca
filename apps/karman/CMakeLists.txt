# Add application source files
file(GLOB LGCA_KARMAN_SOURCES *.cpp *.cu)
file(GLOB LGCA_KARMAN_HEADERS  *.h *.cuh)

# Specify target and source files to compile from
add_executable(
  ${PROJECT_NAME}-karman
  ${LGCA_KARMAN_SOURCES} ${LGCA_KARMAN_HEADERS}
  ${LIB_LGCA_SOURCES} ${LIB_LGCA_HEADERS}
  ${PNGWRITER_SOURCES} ${PNGWRITER_HEADERS}
  ${TCLAP_HEADERS}
)

# Specify target and libraries to link with
target_link_libraries(
  ${PROJECT_NAME}-karman
  ${CUDA_LIBRARIES}
  ${OpenMP_C_LIBRARIES} ${OpenMP_CXX_LIBRARIES}
  ${TBB_LIBRARIES}
  ${PNG_LIBRARIES}
  ${FREETYPE_LIBRARIES}
  ${VTK_LIBRARIES}
)
