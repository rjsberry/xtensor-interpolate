# xtensor-fitpack: https://github.com/rjsberry/xtensor-fitpack
#
# Copyright (C) 2018, Richard Berry <rjsberry@protonmail.com>
#
# Distributed under the terms of BSD 2-Clause "simplified" license. (See
# accompanying file LICENSE, or copy at
# https://github.com/rjsberry/xtensor-fitpack/blob/master/LICENSE)
#

enable_language(Fortran)

if(${CMAKE_Fortran_COMPILER} MATCHES "ifort.*")
    set(CMAKE_Fortran_FLAGS_RELEASE "${CMAKE_Fortran_FLAGS} -fltconsistency")
    set(CMAKE_Fortran_FLAGS_DEBUG   "${CMAKE_Fortran_FLAGS} -fltconsistency")
endif()

set(fitpack_SRC
    FITPACK/bispeu.f
    FITPACK/bispev.f
    FITPACK/clocur.f
    FITPACK/cocosp.f
    FITPACK/concon.f
    FITPACK/concur.f
    FITPACK/cualde.f
    FITPACK/curev.f
    FITPACK/curfit.f
    FITPACK/dblint.f
    FITPACK/evapol.f
    FITPACK/fourco.f
    FITPACK/fpader.f
    FITPACK/fpadno.f
    FITPACK/fpadpo.f
    FITPACK/fpback.f
    FITPACK/fpbacp.f
    FITPACK/fpbfout.f
    FITPACK/fpbisp.f
    FITPACK/fpbspl.f
    FITPACK/fpchec.f
    FITPACK/fpched.f
    FITPACK/fpchep.f
    FITPACK/fpclos.f
    FITPACK/fpcoco.f
    FITPACK/fpcons.f
    FITPACK/fpcosp.f
    FITPACK/fpcsin.f
    FITPACK/fpcurf.f
    FITPACK/fpcuro.f
    FITPACK/fpcyt1.f
    FITPACK/fpcyt2.f
    FITPACK/fpdeno.f
    FITPACK/fpdisc.f
    FITPACK/fpfrno.f
    FITPACK/fpgivs.f
    FITPACK/fpgrdi.f
    FITPACK/fpgrpa.f
    FITPACK/fpgrre.f
    FITPACK/fpgrsp.f
    FITPACK/fpinst.f
    FITPACK/fpintb.f
    FITPACK/fpknot.f
    FITPACK/fpopdi.f
    FITPACK/fpopsp.f
    FITPACK/fporde.f
    FITPACK/fppara.f
    FITPACK/fppasu.f
    FITPACK/fpperi.f
    FITPACK/fppocu.f
    FITPACK/fppogr.f
    FITPACK/fppola.f
    FITPACK/fprank.f
    FITPACK/fprati.f
    FITPACK/fpregr.f
    FITPACK/fprota.f
    FITPACK/fprppo.f
    FITPACK/fprpsp.f
    FITPACK/fpseno.f
    FITPACK/fpspgr.f
    FITPACK/fpsphe.f
    FITPACK/fpsuev.f
    FITPACK/fpsurf.f
    FITPACK/fpsysy.f
    FITPACK/fptrnp.f
    FITPACK/fptrpe.f
    FITPACK/insert.f
    FITPACK/parcur.f
    FITPACK/parder.f
    FITPACK/pardeu.f
    FITPACK/parsur.f
    FITPACK/percur.f
    FITPACK/pogrid.f
    FITPACK/polar.f
    FITPACK/profil.f
    FITPACK/regrid.f
    FITPACK/spalde.f
    FITPACK/spgrid.f
    FITPACK/sphere.f
    FITPACK/splder.f
    FITPACK/splev.f
    FITPACK/splint.f
    FITPACK/sproot.f
    FITPACK/surev.f
    FITPACK/surfit.f
)

add_library(fitpack STATIC ${fitpack_SRC})
set(FITPACK_LIBRARY fitpack PARENT_SCOPE)

enable_language(CXX)

set(FC_SYMBOLS
    curfit
    fpcurf
    splev
)

include(FortranCInterface)
FortranCInterface_HEADER(${XTENSOR_XFITPACK_INTERFACE_INCLUDE_DIR}/FCMangle.h
    MACRO_NAMESPACE "FC_"
    SYMBOL_NAMESPACE "_fc_"
    SYMBOLS ${FC_SYMBOLS}
)
