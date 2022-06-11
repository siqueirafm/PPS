set (DCEL_INC ${CMAKE_CURRENT_SOURCE_DIR}/headers)

file (GLOB HPP_FILES ${DCEL_INC}/*.hpp)

add_library(dcel INTERFACE)

target_sources(dcel INTERFACE ${HPP_FILES})

target_include_directories(dcel SYSTEM INTERFACE ${DCEL_INC})

source_group(
  TREE "${DCEL_INC}"
  PREFIX "Header Files"
  FILES ${HPP_FILES})

install(FILES ${HPP_FILES} DESTINATION include/dcel)
