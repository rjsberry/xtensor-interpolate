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

if(DOWNLOAD_GTEST OR GTEST_SRC_DIR)
    if(DOWNLOAD_GTEST)
        configure_file(downloadGTest.cmake.in googletest-download/CMakeLists.txt)
    else()
        configure_file(copyGTest.cmake.in googletest-download/CMakeLists.txt)
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
                    RESULT_VARIABLE result
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
    if(result)
        message(FATAL_ERROR "CMake step for googletest failed: ${result}")
    endif()
    execute_process(COMMAND ${CMAKE_COMMAND} --build .
                    RESULT_VARIABLE result
                    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download)
    if(result)
        message(FATAL_ERROR "Build step for googletest failed: ${result}")
    endif()

    add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
                     ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
                     EXCLUDE_FROM_ALL)

    set(GTEST_INCLUDE_DIRS "${gtest_SOURCE_DIR}/include")
    set(GTEST_BOTH_LIBRARIES gtest_main gtest)
else()
    find_package(GTest REQUIRED)
endif()

if(WIN32 AND (MINGW OR MSYS OR MSVC))
    set(GTEST_SINGLE_THREAD TRUE)
endif()

if(NOT GTEST_SINGLE_THREAD)
    find_package(Threads)
endif()

find_package(xtensor REQUIRED)

include_directories(${XTENSOR_INTERPOLATE_INCLUDE_DIR})
include_directories(${XTENSOR_INCLUDE_DIR})
include_directories(${GTEST_INCLUDE_DIRS})

set(XTENSOR_INTERPOLATE_TESTS
    main.cpp
    test_xinterpolate.cpp
   )

set(XTENSOR_INTERPOLATE_TARGET test_xtensor_interpolate)

add_executable(${XTENSOR_INTERPOLATE_TARGET}
    ${XTENSOR_INTERPOLATE_TESTS}
    ${XTENSOR_INTERPOLATE_HEADER}
    ${XTENSOR_HEADERS}
   )

if(DOWNLOAD_GTEST OR GTEST_SRC_DIR)
    add_dependencies(${XTENSOR_INTERPOLATE_TARGET} gtest_main)
endif()

if(GTEST_SINGLE_THREAD)
    target_link_libraries(${XTENSOR_INTERPOLATE_TARGET}
        ${FITPACK_LIBRARY}
        ${GTEST_BOTH_LIBRARIES}
       )
else()
    target_link_libraries(${XTENSOR_INTERPOLATE_TARGET}
        ${FITPACK_LIBRARY}
        ${GTEST_BOTH_LIBRARIES}
        ${CMAKE_THREAD_LIBS_INIT}
       )
endif()

add_custom_target(xtest COMMAND test_xtensor_interpolate DEPENDS ${XTENSOR_INTERPOLATE_TARGET})
