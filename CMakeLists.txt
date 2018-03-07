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

set(XTENSOR_FITPACK_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
set(XTENSOR_FITPACK_HEADER
    ${XTENSOR_FITPACK_INCLUDE_DIR}/xtensor-fitpack/xfitpack.hpp
    ${XTENSOR_FITPACK_INCLUDE_DIR}/xtensor-fitpack/xfitpack_utils.hpp
)
project(xtensor-fitpack)

find_package(xtensor REQUIRED)

enable_language(Fortran)

if(${CMAKE_Fortran_COMPILER} MATCHES "ifort.*")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -fltconsistency")
    set(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -fltconsistency")
endif()

set(fitpack_SRC
    third_party/fitpack/bispeu.f
    third_party/fitpack/bispev.f
    third_party/fitpack/clocur.f
    third_party/fitpack/cocosp.f
    third_party/fitpack/concon.f
    third_party/fitpack/concur.f
    third_party/fitpack/cualde.f
    third_party/fitpack/curev.f
    third_party/fitpack/curfit.f
    third_party/fitpack/dblint.f
    third_party/fitpack/evapol.f
    third_party/fitpack/fourco.f
    third_party/fitpack/fpader.f
    third_party/fitpack/fpadno.f
    third_party/fitpack/fpadpo.f
    third_party/fitpack/fpback.f
    third_party/fitpack/fpbacp.f
    third_party/fitpack/fpbfout.f
    third_party/fitpack/fpbisp.f
    third_party/fitpack/fpbspl.f
    third_party/fitpack/fpchec.f
    third_party/fitpack/fpched.f
    third_party/fitpack/fpchep.f
    third_party/fitpack/fpclos.f
    third_party/fitpack/fpcoco.f
    third_party/fitpack/fpcons.f
    third_party/fitpack/fpcosp.f
    third_party/fitpack/fpcsin.f
    third_party/fitpack/fpcurf.f
    third_party/fitpack/fpcuro.f
    third_party/fitpack/fpcyt1.f
    third_party/fitpack/fpcyt2.f
    third_party/fitpack/fpdeno.f
    third_party/fitpack/fpdisc.f
    third_party/fitpack/fpfrno.f
    third_party/fitpack/fpgivs.f
    third_party/fitpack/fpgrdi.f
    third_party/fitpack/fpgrpa.f
    third_party/fitpack/fpgrre.f
    third_party/fitpack/fpgrsp.f
    third_party/fitpack/fpinst.f
    third_party/fitpack/fpintb.f
    third_party/fitpack/fpknot.f
    third_party/fitpack/fpopdi.f
    third_party/fitpack/fpopsp.f
    third_party/fitpack/fporde.f
    third_party/fitpack/fppara.f
    third_party/fitpack/fppasu.f
    third_party/fitpack/fpperi.f
    third_party/fitpack/fppocu.f
    third_party/fitpack/fppogr.f
    third_party/fitpack/fppola.f
    third_party/fitpack/fprank.f
    third_party/fitpack/fprati.f
    third_party/fitpack/fpregr.f
    third_party/fitpack/fprota.f
    third_party/fitpack/fprppo.f
    third_party/fitpack/fprpsp.f
    third_party/fitpack/fpseno.f
    third_party/fitpack/fpspgr.f
    third_party/fitpack/fpsphe.f
    third_party/fitpack/fpsuev.f
    third_party/fitpack/fpsurf.f
    third_party/fitpack/fpsysy.f
    third_party/fitpack/fptrnp.f
    third_party/fitpack/fptrpe.f
    third_party/fitpack/insert.f
    third_party/fitpack/parcur.f
    third_party/fitpack/parder.f
    third_party/fitpack/pardeu.f
    third_party/fitpack/parsur.f
    third_party/fitpack/percur.f
    third_party/fitpack/pogrid.f
    third_party/fitpack/polar.f
    third_party/fitpack/profil.f
    third_party/fitpack/regrid.f
    third_party/fitpack/spalde.f
    third_party/fitpack/spgrid.f
    third_party/fitpack/sphere.f
    third_party/fitpack/splder.f
    third_party/fitpack/splev.f
    third_party/fitpack/splint.f
    third_party/fitpack/sproot.f
    third_party/fitpack/surev.f
    third_party/fitpack/surfit.f
)

add_library(fitpack STATIC ${fitpack_SRC})
set(FITPACK_LIBRARY fitpack)

enable_language(CXX)

include(FortranCInterface)
FortranCInterface_HEADER(FCMangle.h
    MACRO_NAMESPACE "FC_"
    SYMBOL_NAMESPACE ""
    SYMBOLS curfit
)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

enable_testing(true)
add_subdirectory(test/)
