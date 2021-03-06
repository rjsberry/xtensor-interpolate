# xtensor-interpolate: https://github.com/rjsberry/xtensor-interpolate
#
# Copyright (C) 2018, Richard Berry <rjsberry@protonmail.com>
#
# Distributed under the terms of BSD 2-Clause "simplified" license. (See
# accompanying file LICENSE, or copy at
# https://github.com/rjsberry/xtensor-interpolate/blob/master/LICENSE)
#

cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 14)

project(xtensor-interpolate)
set(${PROJECT_NAME}_VERSION 0.1.0)

set(XTENSOR_INTERPOLATE_INCLUDE_DIR
    ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(XTENSOR_INTERPOLATE_HEADERS
    ${XTENSOR_INTERPOLATE_INCLUDE_DIR}/xtensor-interpolate/xinterpolate.hpp)

set(XTENSOR_FITPACK_INTERFACE_INCLUDE_DIR
    ${XTENSOR_INTERPOLATE_INCLUDE_DIR}/xfitpack)

set(XTENSOR_FITPACK_INTERFACE_HEADERS
    ${XTENSOR_FITPACK_INTERFACE_INCLUDE_DIR}/xfcdec.h
    ${XTENSOR_FITPACK_INTERFACE_INCLUDE_DIR}/xfcmangle.h)

find_package(xtensor REQUIRED)

# FITPACK interface

add_subdirectory(${XTENSOR_FITPACK_INTERFACE_INCLUDE_DIR})

# Tests

OPTION(BUILD_TESTS "xtensor-interpolate test suite" OFF)
OPTION(DOWNLOAD_GTEST "build gtest from downloaded sources" OFF)

if(DOWNLOAD_GTEST OR GTEST_SRC_DIR)
    set(BUILD_TESTS ON)
endif()

if(BUILD_TESTS)
    enable_testing(true)
    add_subdirectory(test)
endif()

# Installation

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

install(FILES ${XTENSOR_INTERPOLATE_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/xtensor-interpolate)

install(FILES ${XTENSOR_FITPACK_INTERFACE_HEADERS}
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/xfitpack)

set(XTENSOR_INTERPOLATE_CMAKECONFIG_INSTALL_DIR "${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}" CACHE
    STRING "install path for xtensor-interpolateConfig.cmake")

configure_package_config_file(${PROJECT_NAME}Config.cmake.in
                              "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake"
                              INSTALL_DESTINATION ${XTENSOR_INTERPOLATE_CMAKECONFIG_INSTALL_DIR})

write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
                                 VERSION ${${PROJECT_NAME}_VERSION}
                                 COMPATIBILITY AnyNewerVersion)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}ConfigVersion.cmake
        DESTINATION ${XTENSOR_INTERPOLATE_CMAKECONFIG_INSTALL_DIR})
