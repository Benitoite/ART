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

#include "boxblur.h"
#include "rescale.h"
#include "array2D.h"

namespace rtengine {

template <class T>
void guidedFilter(T **guide, T **src, T **dst, int W, int H, int r, float epsilon, bool multithread, int subsampling=4)
{
    const int ww = W / subsampling;
    const int hh = H / subsampling;

    array2D<float> guide1(ww, hh);
    array2D<float> src1(ww, hh);

    rescaleBilinear(guide, W, H, guide1, ww, hh, multithread);
    rescaleBilinear(src, W, H, src1, ww, hh, multithread);

    float r1 = float(r) / subsampling;

    array2D<float> meanI(ww, hh);
    boxblur(guide1, meanI, r1, r1, ww, hh);

    array2D<float> meanp(ww, hh);
    boxblur(src, meanp, r1, r1, ww, hh);

    enum Op { MUL, DIV, ADD, SUB };

    const auto apply =
        [=](Op op, T **a, T **b, T **res, int w, int h) -> void
        {
#ifdef _OPENMP
            #pragma omp parallel for if (multithread)
#endif
            for (int y = 0; y < h; ++y) {
                for (int x = 0; x < w; ++x) {
                    T r;
                    T aa = a[y][x];
                    T bb = b[y][x];
                    switch (op) {
                    case MUL:
                        r = aa * bb;
                        break;
                    case DIV:
                        r = aa / (bb + epsilon);
                        break;
                    case ADD:
                        r = aa + bb;
                        break;
                    case SUB:
                        r = aa - bb;
                        break;
                    default:
                        assert(false);
                        r = 0;
                        break;
                    }
                    res[y][x] = r;
                }
            }
        };

    array2D<float> corrI(ww, hh);
    apply(MUL, guide1, guide1, corrI, ww, hh);
    boxblur(corrI, corrI, r1, r1, ww, hh);

    array2D<float> corrIp(ww, hh);
    apply(MUL, guide1, src1, corrIp, ww, hh);
    boxblur(corrIp, corrIp, r1, r1, ww, hh);

    array2D<float> varI(ww, hh);
    apply(MUL, meanI, meanI, varI, ww, hh);
    apply(SUB, corrI, varI, varI, ww, hh);
    
    array2D<float> covIp(ww, hh);
    apply(MUL, meanI, meanp, covIp, ww, hh);
    apply(SUB, corrIp, covIp, covIp, ww, hh);

    T **a = covIp;
    apply(DIV, covIp, varI, a, ww, hh);

    T **b = meanp;
    apply(MUL, a, meanI, meanI, ww, hh);
    apply(SUB, meanp, meanI, b, ww, hh);

    boxblur(a, a, r1, r1, ww, hh);
    boxblur(b, b, r1, r1, ww, hh);

    array2D<float> aa(W, H);
    rescaleBilinear(a, ww, hh, aa, W, H, multithread);
    
    array2D<float> bb(W, H);
    rescaleBilinear(b, ww, hh, bb, W, H, multithread);

    apply(MUL, aa, guide, dst, W, H);
    apply(ADD, dst, bb, dst, W, H);
}

} // namespace rtengine
