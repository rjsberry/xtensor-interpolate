# xtensor-fitpack: https://github.com/rjsberry/xtensor-fitpack
#
# Copyright (C) 2018, Richard Berry <rjsberry@protonmail.com>
#
# Distributed under the terms of BSD 2-Clause "simplified" license. (See
# accompanying file LICENSE, or copy at
# https://github.com/rjsberry/xtensor-fitpack/blob/master/LICENSE)
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

include_directories(${XTENSOR_FITPACK_INCLUDE_DIR})
include_directories(${XTENSOR_INCLUDE_DIR})
include_directories(${GTEST_INCLUDE_DIRS})

set(XTENSOR_FITPACK_TESTS
    main.cpp
    test_univariate_spline.cpp
)

set(XTENSOR_FITPACK_TARGET test_xtensor_fitpack)

add_executable(${XTENSOR_FITPACK_TARGET} ${XTENSOR_FITPACK_TESTS} ${XTENSOR_FITPACK_HEADER})
add_dependencies(${XTENSOR_FITPACK_TARGET} gtest_main)

target_link_libraries(${XTENSOR_FITPACK_TARGET}
    ${GTEST_BOTH_LIBRARIES}
    ${FITPACK_LIBRARY}
    ${CMAKE_THREAD_LIBS_INIT}
)

add_custom_target(xtest COMMAND test_xtensor_fitpack)
