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

As well as the xtensor headers, all installation methods place the ``cmake``
project configuration file in the appropriate location so that other projects
can use cmake's ``find_package`` to locate ``xtensor-interpolate``.

.. image:: cmake.svg

CMake
-----

On Unix platforms, from the source directory:

.. code::

    mkdir build
    cd build
    cmake ..
    make install
