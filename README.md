# xtensor-interpolate

[![Travis](https://travis-ci.org/rjsberry/xtensor-interpolate.svg?branch=master)](https://travis-ci.org/rjsberry/xtensor-interpolate)

`xtensor-interpolate` is an ***unofficial*** extension to the `xtensor` library offering spline interpolation via
internal bindings to Paul Dierckx's Fortran package, [FITPACK](http://www.netlib.org/dierckx/).

## Introduction

The FITPACK library is used in [SciPy](https://github.com/scipy/scipy) to build the back-end of `scipy.interpolate`;
`xtensor-interpolate` aims to provide a familiar API to anyone experienced with this Python package.

## Development

`xtensor-interpolate` is under heavy development as is **not** currently considered usable.

## License

This software is licensed under the BSD-2-Clause license. See the
[LICENSE](https://github.com/rjsberry/xtensor-interpolate/blob/master/LICENSE) file for details.
