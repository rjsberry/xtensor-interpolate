# xtensor-fitpack: https://github.com/rjsberry/xtensor-fitpack
#
# Copyright (C) 2018, Richard Berry <rjsberry@protonmail.com>
#
# Distributed under the terms of BSD 2-Clause "simplified" license. (See
# accompanying file LICENSE, or copy at
# https://github.com/rjsberry/xtensor-fitpack/blob/master/LICENSE)
#

cmake_minimum_required(VERSION 3.1)
set(CMAKE_CXX_STANDARD 14)


project(xtensor-fitpack)

set(XTENSOR_FITPACK_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(XTENSOR_FITPACK_HEADER
    ${XTENSOR_FITPACK_INCLUDE_DIR}/xtensor-fitpack/xfitpack.hpp
)
set(XTENSOR_XFITPACK_INTERFACE_INCLUDE_DIR
    ${XTENSOR_FITPACK_INCLUDE_DIR}/xfitpack
)

find_package(xtensor REQUIRED)

add_subdirectory(${XTENSOR_XFITPACK_INTERFACE_INCLUDE_DIR})

enable_testing(true)
add_subdirectory(test)
