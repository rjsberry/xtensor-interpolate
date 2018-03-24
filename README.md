# xtensor-interpolate

[![Travis](https://travis-ci.org/rjsberry/xtensor-interpolate.svg?branch=master)](https://travis-ci.org/rjsberry/xtensor-interpolate)
[![AppVeyor](https://ci.appveyor.com/api/projects/status/hp9uya98ca0oavlk?svg=true)](https://ci.appveyor.com/project/rjsberry/xtensor-interpolate)
[![Read the Docs](https://readthedocs.org/projects/xtensor-interpolate/badge/?version=latest)](http://xtensor-interpolate.readthedocs.io/en/latest/?badge=latest)

`xtensor-interpolate` is an ***unofficial*** extension to the `xtensor` library offering spline interpolation via
internal bindings to Paul Dierckx's Fortran package, [FITPACK](http://www.netlib.org/dierckx/).

## Introduction

The FITPACK library is used in [SciPy](https://github.com/scipy/scipy) to build the back-end of `scipy.interpolate`;
`xtensor-interpolate` aims to provide a familiar API to anyone experienced with this Python package.

## Installation

### From Source

Installation from source is streamlined with CMake.

```
mkdir build
cd build
cmake ..
make install
```

The `cmake` step of the process can be customized with a number of flags.

- `-DCMAKE_INSTALLATION_PREFIX=...` to change where you want
  `xtensor-interpolate` to install to.

- `-DBUILD_TESTS=ON` to build the unit test suite. This requires `googletest`.

  - You can target the unit tests directly with `make xtest`, or just build the
    executable with `make test_xtensor_interpolate`.

- `-DDOWNLOAD_GTEST=ON` to download and compile `googletest` as part of the
  build step. Note that `googletest` will not actually be installed to your
  system directories, and is just statically linked into the test executable.

  - Use `Dgtest_disable_pthreads=ON` to compile a single-threaded version of
    `googletest`. This is **required** when using Windows and MinGW.

## Usage

### Interpolation using the functional interface

```cpp
#include <cmath>
#include "xtensor/xtensor.hpp"
#include "xtensor-interpolate/xinterpolate.hpp"

xt::xtensor<double, 1> x =
    { -3.141, -2.443, -1.745, -1.047, -0.349, 0.349,  1.047,  1.745,  2.443,  3.141 };

xt::xtensor<double, 1> y =
    { 0., -0.643, -0.985, -0.866, -0.342,  0.342,  0.866,  0.985, 0.643,  0. };

xt::xtensor<double, 1> xs = xt::linspace<double>(-M_PI, M_PI, 100);

auto tck = xt::interpolate::splrep(x, y);

auto ys = xt::interpolate::splev(xs, tck);
```

Plotted, this data looks like:

![](assets/usage.png?raw=true)

## Development

`xtensor-interpolate` is under heavy development as is **not** currently considered usable.

## License

This software is licensed under the BSD-2-Clause license. See the
[LICENSE](https://github.com/rjsberry/xtensor-interpolate/blob/master/LICENSE) file for details.
