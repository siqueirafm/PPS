set (SAMPLERPNTAPP_SRC ${CMAKE_CURRENT_SOURCE_DIR}/sources)

set (DCELLIB_INC ${CMAKE_CURRENT_SOURCE_DIR}/../../library/dcel/headers)
set (OFFLIB_INC ${CMAKE_CURRENT_SOURCE_DIR}/../../library/off/headers)
set (PPSLIB_INC ${CMAKE_CURRENT_SOURCE_DIR}/../../library/pps/headers)
set (PPSFROMPNT_INC ${CMAKE_CURRENT_SOURCE_DIR}/../../library/ppsfrompnt/headers)
set (UTILSLIB_INC ${CMAKE_CURRENT_SOURCE_DIR}/../../library/utils/headers)

file (GLOB CPP_FILES ${SAMPLERPNTAPP_SRC}/*.cpp)

add_executable(sampler-pnt "")

target_compile_features(sampler-pnt PRIVATE cxx_std_17)

if (MSVC)
  if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(sampler-pnt PRIVATE /W4 /WX- /O2 /MD)
  else ()
    target_compile_options(sampler-pnt PRIVATE /W4 /WX- /MDd)
  endif ()
elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang" OR CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  if (NOT CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_options(sampler-pnt PRIVATE -Wall -Wextra -Werror -pedantic -O3)
  else ()
    target_compile_options(sampler-pnt PRIVATE -Wall -Wextra -Werror -pedantic)
  endif ()
endif ()

target_sources(sampler-pnt PRIVATE ${CPP_FILES})

target_include_directories(
  sampler-pnt
    PRIVATE
       ${DCELLIB_INC}
       ${OFFLIB_INC}
       ${PPSLIB_INC}
       ${PPSFROMPNT_INC}
       ${UTILSLIB_INC}
)

target_link_libraries(
  sampler-pnt
    PRIVATE
      ppsfrompnt
      pps
      off
      dcel
      utils
)

install(TARGETS sampler-pnt DESTINATION bin)
