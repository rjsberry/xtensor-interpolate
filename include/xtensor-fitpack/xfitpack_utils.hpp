#ifndef INCLUDE_XTENSOR_FITPACK_XFITPACK_UTILS_HPP_
#define INCLUDE_XTENSOR_FITPACK_XFITPACK_UTILS_HPP_

// xtensor-fitpack: https://github.com/rjsberry/xtensor-fitpack
//
// Copyright (C) 2018, Richard Berry <rjsberry@protonmail.com>
//
// Distributed under the terms of BSD 2-Clause "simplified" license. (See
// accompanying file LICENSE, or copy at
// https://github.com/rjsberry/xtensor-fitpack/blob/master/LICENSE)
//

#include <algorithm>
#include <memory>

namespace xt {

template <typename T>
auto allocate_homogeneous_array(std::size_t len, T filler) {
    auto a = std::make_unique<T[]>(len);
    std::fill(a.get(), a.get() + len, filler);

    return a;
}

template <typename T>
void trim_array_ptr(std::unique_ptr<T[]>& x,
                    std::size_t orig_len,
                    std::size_t new_len) {
    std::copy_n(x.get(), std::min(orig_len, new_len), x.get());
}

}  // namespace xt

#endif  // INCLUDE_XTENSOR_FITPACK_XFITPACK_UTILS_HPP_
