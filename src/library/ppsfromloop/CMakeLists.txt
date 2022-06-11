set (PPSFROMLOOPLIB_INC ${CMAKE_CURRENT_SOURCE_DIR}/headers)
set (PPSFROMLOOPLIB_SRC ${CMAKE_CURRENT_SOURCE_DIR}/sources)
set (PPSFROMLOOPLIB_DAT ${CMAKE_CURRENT_SOURCE_DIR}/dat)

file (GLOB HPP_FILES ${PPSFROMLOOPLIB_INC}/*.hpp)
file (GLOB CPP_FILES ${PPSFROMLOOPLIB_SRC}/*.cpp)
file (GLOB DAT_FILES ${PPSFROMLOOPLIB_DAT}/*.dat)

# Make a static library (omit STATIC if you want user setting to decide)
add_library(ppsfromloop
  STATIC
    ${HPP_FILES}
    ${CPP_FILES}
  )

set (DCELLIB_INC ${CMAKE_CURRENT_SOURCE_DIR}/../dcel/headers)
set (PPSLIB_INC ${CMAKE_CURRENT_SOURCE_DIR}/../pps/headers)
set (UTILSLIB_INC ${CMAKE_CURRENT_SOURCE_DIR}/../utils/headers)

# Library need this directories, and its users will need them too.
target_include_directories(ppsfromloop
  PUBLIC 
    $<INSTALL_INTERFACE:include/ppsfromloop>
    ${PPSFROMLOOPLIB_INC}
  PRIVATE
    ${DCELLIB_INC}
    ${PPSLIB_INC}
    ${UTILSLIB_INC}
)

# All users of this library will need at list C++ 17
target_compile_features(ppsfromloop PRIVATE cxx_std_17)

if (MSVC)
  if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(ppsfromloop PRIVATE /W4 /WX- /O2 /MD)
  else ()
    target_compile_options(ppsfromloop PRIVATE /W4 /WX- /MDd)
  endif ()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(ppsfromloop PRIVATE -Wall -Wextra -Werror -pedantic -O3)
  else ()
    target_compile_options(ppsfromloop PRIVATE -Wall -Wextra -Werror -pedantic)
  endif ()
endif ()

source_group(
  TREE "${PPSFROMLOOPLIB_INC}"
  PREFIX "Header Files"
  FILES ${HPP_FILES})

source_group(
  TREE "${PPSFROMLOOPLIB_SRC}"
  PREFIX "Source Files"
  FILES ${CPP_FILES})

install(TARGETS ppsfromloop DESTINATION lib)
install(FILES ${HPP_FILES} DESTINATION include/ppsfromloop)
install(FILES ${DAT_FILES} DESTINATION dat)
