/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2018 Alberto Griggio <alberto.griggio@gmail.com>
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include "array2D.h"
#include "sleef.c"

namespace rtengine {

void guidedFilter(const array2D<float> &guide, const array2D<float> &src, array2D<float> &dst, int r, float epsilon, bool multithread, int subsampling=0);

inline void guidedFilterLog(float base, array2D<float> &chan, int r, float eps, bool multithread, int subsampling=0)
{
#ifdef _OPENMP
#    pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < chan.height(); ++y) {
        for (int x = 0; x < chan.width(); ++x) {
            chan[y][x] = xlin2log(max(chan[y][x], 0.f), base);
        }
    }

    guidedFilter(chan, chan, chan, r, eps, multithread, subsampling);

#ifdef _OPENMP
#    pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < chan.height(); ++y) {
        for (int x = 0; x < chan.width(); ++x) {
            chan[y][x] = xlog2lin(max(chan[y][x], 0.f), base);
        }
    }
}

} // namespace rtengine
