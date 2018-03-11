#ifndef INCLUDE_XTENSOR_FITPACK_INTERFACE_DECLARATIONS_H_
#define INCLUDE_XTENSOR_FITPACK_INTERFACE_DECLARATIONS_H_

// xtensor-interpolate: https://github.com/rjsberry/xtensor-interpolate
//
// Copyright (C) 2018, Richard Berry <rjsberry@protonmail.com>
//
// Distributed under the terms of BSD 2-Clause "simplified" license. (See
// accompanying file LICENSE, or copy at
// https://github.com/rjsberry/xtensor-interpolate/blob/master/LICENSE)
//

// Declarations of free C functions that interface Fortran code.

#include "FCMangle.h"

extern "C" {

void _fc_curfit(int* iopt, int* m, double* x, double* y, double* w,
                double* xb, double* xe, int* k, double* s, int* nest,
                int* n, double* t, double* c, double* fp, double* wrk,
                int* lwrk, int* iwrk, int* ier);

void _fc_fpcurf(int* iopt, double* x, double* y, double* w, int* m,
                double* xb, double* xe, int* k, double* s, int* nest,
                double* tol, int* maxit, int* k1, int* k2, int* n,
                double* t, double* c, double* fp, double* fpint,
                double* wrk, double* wrkn, double* wrknk2, double* wrkn2k2,
                double* wrkn3k2, int* nrdata, int* ier);

void _fc_splev(double* t, int* n, double* c, int* k, double* x, double* y,
               int* m, int* e, int* ier);

}  // extern "C"

#endif  // INCLUDE_XTENSOR_FITPACK_INTERFACE_DECLARATIONS_H_
