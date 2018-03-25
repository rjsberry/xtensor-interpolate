.. Copyright (C) 2018, Richard Berry

   Distributed under the terms of the BSD-2-Clause License.

   The full license is in the file LICENSE, distributed with this software.

.. raw:: html

   <style>
   .rst-content .section>img {
       width: 30px;
       margin-bottom: 0;
       margin-top: 0;
       margin-right: 15px;
       margin-left: 15px;
       float: left;
   }
   </style>

Installation
============

.. image:: cmake.svg

CMake
-----

.. code::

    mkdir build
    cd build
    cmake ..
    make install

This will install the ``xtensor-interpolate`` headers, compile and install the
FITPACK interface libraries (``.so`` on Linux, ``.dll`` on Windows), and place
the CMake project configuration file in the appropriate location so that other
projects can use CMake's ``find_package(xtensor-interpolate)``.
