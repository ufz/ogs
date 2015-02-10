# This is Googletest 1.6 compressed to 2 source files

if (MSVC)
  add_definitions (/D_VARIADIC_MAX=10)
endif ()

if (HAVE_PTHREADS)
  add_definitions (-DGTEST_HAS_PTHREAD=1)
else ()
  add_definitions (-DGTEST_HAS_PTHREAD=0)
endif ()

include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/include)
add_library (GTest src/gtest-all.cc)

set_target_properties(GTest PROPERTIES FOLDER Testing)
