set (SAMPLERLOOPAPP_SRC ${CMAKE_CURRENT_SOURCE_DIR}/sources)

set (DCELLIB_INC ${CMAKE_CURRENT_SOURCE_DIR}/../../library/dcel/headers)
set (OFFLIB_INC ${CMAKE_CURRENT_SOURCE_DIR}/../../library/off/headers)
set (PPSLIB_INC ${CMAKE_CURRENT_SOURCE_DIR}/../../library/pps/headers)
set (PPSFROMLOOP_INC ${CMAKE_CURRENT_SOURCE_DIR}/../../library/ppsfromloop/headers)
set (UTILSLIB_INC ${CMAKE_CURRENT_SOURCE_DIR}/../../library/utils/headers)

file (GLOB CPP_FILES ${SAMPLERLOOPAPP_SRC}/*.cpp)

add_executable(sampler-loop "")

target_compile_features(sampler-loop PRIVATE cxx_std_17)

if (MSVC)
  if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(sampler-loop PRIVATE /W4 /WX- /O2 /MD)
  else ()
    target_compile_options(sampler-loop PRIVATE /W4 /WX- /MDd)
  endif ()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(sampler-loop PRIVATE -Wall -Wextra -Werror -pedantic -O3)
  else ()
    target_compile_options(sampler-loop PRIVATE -Wall -Wextra -Werror -pedantic)
  endif ()
endif ()

target_sources(sampler-loop PRIVATE ${CPP_FILES})

target_include_directories(
  sampler-loop
    PRIVATE
       ${DCELLIB_INC}
       ${OFFLIB_INC}
       ${PPSLIB_INC}
       ${PPSFROMLOOP_INC}
       ${UTILSLIB_INC}
)

target_link_libraries(
  sampler-loop
    PRIVATE
      ppsfromloop
      pps
      off
      dcel
      utils
)

install(TARGETS sampler-loop DESTINATION bin)
