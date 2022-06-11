set (PPSFROMPNTLIB_INC ${CMAKE_CURRENT_SOURCE_DIR}/headers)
set (PPSFROMPNTLIB_SRC ${CMAKE_CURRENT_SOURCE_DIR}/sources)

file (GLOB HPP_FILES ${PPSFROMPNTLIB_INC}/*.hpp)
file (GLOB CPP_FILES ${PPSFROMPNTLIB_SRC}/*.cpp)

# Make a static library (omit STATIC if you want user setting to decide)
add_library(ppsfrompnt
  STATIC
    ${HPP_FILES}
    ${CPP_FILES}
  )

set (DCELLIB_INC ${CMAKE_CURRENT_SOURCE_DIR}/../dcel/headers)
set (PPSLIB_INC ${CMAKE_CURRENT_SOURCE_DIR}/../pps/headers)
set (UTILSLIB_INC ${CMAKE_CURRENT_SOURCE_DIR}/../utils/headers)

# Library need this directories, and its users will need them too.
target_include_directories(ppsfrompnt
  PUBLIC 
    $<INSTALL_INTERFACE:include/ppsfrompnt>
    ${PPSFROMPNTLIB_INC}
  PRIVATE
    ${DCELLIB_INC}
    ${PPSLIB_INC}
    ${UTILSLIB_INC}
)

# All users of this library will need at list C++ 17
target_compile_features(ppsfrompnt PRIVATE cxx_std_17)

if (MSVC)
  if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(ppsfrompnt PRIVATE /W4 /WX- /O2 /MD)
  else ()
    target_compile_options(ppsfrompnt PRIVATE /W4 /WX- /MDd)
  endif ()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(ppsfrompnt PRIVATE -Wall -Wextra -Werror -pedantic -O3)
  else ()
    target_compile_options(ppsfrompnt PRIVATE -Wall -Wextra -Werror -pedantic)
  endif ()
endif ()

source_group(
  TREE "${PPSFROMPNTLIB_INC}"
  PREFIX "Header Files"
  FILES ${HPP_FILES})

source_group(
  TREE "${PPSFROMPNTLIB_SRC}"
  PREFIX "Source Files"
  FILES ${CPP_FILES})

install(TARGETS ppsfrompnt DESTINATION lib)
install(FILES ${HPP_FILES} DESTINATION include/ppsfrompnt)
