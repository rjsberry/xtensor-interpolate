.. Copyright (C) 2018, Richard Berry

   Distributed under the terms of the BSD-2-Clause License.

   The full license is in the file LICENSE, distributed with this software.


The ``xtensor-interpolate`` API is closely modelled after Python's ``scipy.interpolate``.

API Reference
=============

Defined in ``xtensor-interpolate/xinterpolate.hpp``.

.. doxygenfunction:: xt::interpolate::splrep(const xexpression<E1>&, const xexpression<E2>&, const xexpression<E3>&, double, double, int, double)
    :project: xtensor-interpolate

.. doxygenfunction:: xt::interpolate::splrep(const xexpression<E1>&, const xexpression<E2>&, int, double)
    :project: xtensor-interpolate

.. doxygenfunction:: xt::interpolate::splev
    :project: xtensor-interpolate

.. doxygenfunction:: xt::interpolate::splint
    :project: xtensor-interpolate

.. doxygenfunction:: xt::interpolate::spalde
    :project: xtensor-interpolate

.. doxygenfunction:: xt::interpolate::splder
    :project: xtensor-interpolate
