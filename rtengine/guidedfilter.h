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

    enum Op { MUL, DIV, ADD, SUB, ADDMUL, SUBMUL };

    const auto apply =
        [=](Op op, int w, int h, T **res, T **a, T **b, T **c=nullptr) -> void
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
                    case ADDMUL:
                        r = aa * bb + c[y][x];
                        break;
                    case SUBMUL:
                        r = c[y][x] - (aa * bb);
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
    apply(MUL, ww, hh, corrI, guide1, guide1);
    boxblur(corrI, corrI, r1, r1, ww, hh);

    array2D<float> corrIp(ww, hh);
    apply(MUL, ww, hh, corrIp, guide1, src1);
    boxblur(corrIp, corrIp, r1, r1, ww, hh);

    array2D<float> varI(ww, hh);
    apply(SUBMUL, ww, hh, varI, meanI, meanI, corrI);
    
    array2D<float> covIp(ww, hh);
    apply(SUBMUL, ww, hh, covIp, meanI, meanp, corrIp);

    T **a = covIp;
    apply(DIV, ww, hh, a, covIp, varI);

    T **b = meanp;
    apply(SUBMUL, ww, hh, b, a, meanI, meanp);

    boxblur(a, a, r1, r1, ww, hh);
    boxblur(b, b, r1, r1, ww, hh);

    array2D<float> aa(W, H);
    rescaleBilinear(a, ww, hh, aa, W, H, multithread);
    
    array2D<float> bb(W, H);
    rescaleBilinear(b, ww, hh, bb, W, H, multithread);

    apply(ADDMUL, W, H, dst, aa, guide, bb);
}

} // namespace rtengine
