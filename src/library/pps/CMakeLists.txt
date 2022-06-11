set (PPSLIB_INC ${CMAKE_CURRENT_SOURCE_DIR}/headers)
set (PPSLIB_SRC ${CMAKE_CURRENT_SOURCE_DIR}/sources)

file (GLOB HPP_FILES ${PPSLIB_INC}/*.hpp)
file (GLOB CPP_FILES ${PPSLIB_SRC}/*.cpp)

# Make a static library (omit STATIC if you want user setting to decide)
add_library(pps
  STATIC
    ${HPP_FILES}
    ${CPP_FILES}
  )

set (UTILSLIB_INC ${CMAKE_CURRENT_SOURCE_DIR}/../utils/headers)

# Library need this directories, and its users will need them too.
target_include_directories(pps
  PUBLIC 
    $<INSTALL_INTERFACE:include/pps>
    ${PPSLIB_INC}
  PRIVATE
    ${UTILSLIB_INC}
)

# All users of this library will need at list C++ 17
target_compile_features(pps PUBLIC cxx_std_17)

# Set compiler options for MS Visual Studio and everything else
if (MSVC)
  if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(pps PRIVATE /W4 /WX- /O2 /MD)
  else ()
    target_compile_options(pps PRIVATE /W4 /WX- /MDd)
  endif ()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  list (APPEND flags )
  if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(pps PRIVATE -Wall -Wextra -Werror -pedantic -O3)
  else ()
    target_compile_options(pps PRIVATE -Wall -Wextra -Werror -pedantic)
  endif ()
endif ()

source_group(
  TREE "${PPSLIB_INC}"
  PREFIX "Header Files"
  FILES ${HPP_FILES})

source_group(
  TREE "${PPSLIB_SRC}"
  PREFIX "Source Files"
  FILES ${CPP_FILES})

install(TARGETS pps DESTINATION lib)
install(FILES ${HPP_FILES} DESTINATION include/pps)
