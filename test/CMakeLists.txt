# xtensor-interpolate: https://github.com/rjsberry/xtensor-interpolate
#
# Copyright (C) 2018, Richard Berry <rjsberry@protonmail.com>
#
# Distributed under the terms of BSD 2-Clause "simplified" license. (See
# accompanying file LICENSE, or copy at
# https://github.com/rjsberry/xtensor-interpolate/blob/master/LICENSE)
#

find_package(xtensor REQUIRED CONFIG)
set(XTENSOR_INCLUDE_DIR ${xtensor_INCLUDE_DIRS})

configure_file(GTest.cmake.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
                RESULT_VARIABLE result
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
execute_process(COMMAND ${CMAKE_COMMAND} --build .
                RESULT_VARIABLE result
                WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                 ${CMAKE_CURRENT_BINARY_DIR}/googletest-build)
set(GTEST_INCLUDE_DIRS "${gtest_SOURCE_DIR}/include")
set(GTEST_BOTH_LIBRARIES gtest_main gtest)

find_package(Threads)
find_package(xtensor REQUIRED)

include_directories(${XTENSOR_INTERPOLATE_INCLUDE_DIR})
include_directories(${XTENSOR_INCLUDE_DIR})
include_directories(${GTEST_INCLUDE_DIRS})

set(XTENSOR_INTERPOLATE_TESTS
    main.cpp
    test_functional.cpp
    test_univariate_spline.cpp
)

set(XTENSOR_INTERPOLATE_TARGET test_xtensor_interpolate)

add_executable(${XTENSOR_INTERPOLATE_TARGET} ${XTENSOR_INTERPOLATE_TESTS} ${XTENSOR_INTERPOLATE_HEADER})
add_dependencies(${XTENSOR_INTERPOLATE_TARGET} gtest_main)

target_link_libraries(${XTENSOR_INTERPOLATE_TARGET}
    ${GTEST_BOTH_LIBRARIES}
    ${FITPACK_LIBRARY}
    ${CMAKE_THREAD_LIBS_INIT}
)

add_custom_target(xtest COMMAND test_xtensor_interpolate)
