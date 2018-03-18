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

void _fc_curfit(const int* iopt, const int* m, const double* x, const double* y,
                const double* w, const double* xb, const double* xe, const int* k,
                const double* s, const int* nest, int* n, double* t, double* c,
                double* fp, double* wrk, const int* lwrk, int* iwrk, int* ier);

void _fc_fpcurf(const int* iopt, const double* x, const double* y, const double* w,
                const int* m, const double* xb, const double* xe, const int* k,
                const double* s, const int* nest, const double* tol,
                const int* maxit, const int* k1, const int* k2, int* n, double* t,
                double* c, double* fp, double* fpint, double* wrk, double* wrkn,
                double* wrknk2, double* wrkn2k2, double* wrkn3k2, int* nrdata,
                int* ier);

void _fc_spalde(const double* t, const int* n, const double* c, const int* k1,
                const double* x, double* d, int* ier);

double _fc_splint(const double* t, const int* n, const double* c, const int* k1,
                  const double* a, const double* b, double* wrk);

void _fc_splev(const double* t, const int* n, const double* c, const int* k,
               const double* x, double* y, const int* m, const int* e, int* ier);

}  // extern "C"

#endif  // INCLUDE_XTENSOR_FITPACK_INTERFACE_DECLARATIONS_H_
