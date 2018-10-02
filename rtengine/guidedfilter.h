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

inline void guidedFilter(const array2D<float> &guide, const array2D<float> &src, array2D<float> &dst, int r, float epsilon, bool multithread, int subsampling=4)
{
    const int W = src.width();
    const int H = src.height();
    
    const auto dbgsave =
        [](array2D<float> &img, const char *out) -> void
        {
#if 0
            Imagefloat im(img.width(), img.height());
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for (int y = 0; y < im.getHeight(); ++y) {
                for (int x = 0; x < im.getWidth(); ++x) {
                    im.r(y, x) = im.g(y, x) = im.b(y, x) = img[y][x] * 65535.f;
                }
            }
            im.saveTIFF(out, 16);
#endif // if 0
        };
    
    enum Op { MUL, DIV, DIVEPSILON, ADD, SUB, ADDMUL, SUBMUL };

    const auto apply =
        [=](Op op, array2D<float> &res, const array2D<float> &a, const array2D<float> &b, const array2D<float> &c=array2D<float>()) -> void
        {
            const int w = res.width();
            const int h = res.height();
            
#ifdef _OPENMP
            #pragma omp parallel for if (multithread)
#endif
            for (int y = 0; y < h; ++y) {
                for (int x = 0; x < w; ++x) {
                    float r;
                    float aa = a[y][x];
                    float bb = b[y][x];
                    switch (op) {
                    case MUL:
                        r = aa * bb;
                        break;
                    case DIV:
                        r = aa / bb;
                    case DIVEPSILON:
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

    const int ww = W / subsampling;
    const int hh = H / subsampling;

    array2D<float> guide1(ww, hh);
    array2D<float> src1(ww, hh);

    rescaleBilinear(guide, guide1, multithread);
    rescaleBilinear(src, src1, multithread);

    float r1 = float(r) / subsampling;

    array2D<float> N(ww, hh);
#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < hh; ++y) {
        for (int x = 0; x < ww; ++x) {
            N[y][x] = 1.f;
        }
    }
    boxblur<float, float>(N, N, r1, r1, ww, hh);
    

    array2D<float> meanI(ww, hh);
    boxblur<float, float>(guide1, meanI, r1, r1, ww, hh);
    apply(DIV, meanI, meanI, N);
    dbgsave(meanI, "/tmp/meanI.tif");


    array2D<float> meanp(ww, hh);
    boxblur<float, float>(const_cast<array2D<float> &>(src), meanp, r1, r1, ww, hh);
    apply(DIV, meanp, meanp, N);
    dbgsave(meanp, "/tmp/meanp.tif");

    array2D<float> meanIp(ww, hh);
    apply(MUL, meanIp, guide1, src1);
    boxblur<float, float>(meanIp, meanIp, r1, r1, ww, hh);
    apply(DIV, meanIp, meanIp, N);
    dbgsave(meanIp, "/tmp/meanIp.tif");    

    array2D<float> covIp(ww, hh);
    apply(SUBMUL, covIp, meanI, meanp, meanIp);
    dbgsave(covIp, "/tmp/covIp.tif");

    array2D<float> meanII(ww, hh);
    apply(MUL, meanII, guide1, guide1);
    boxblur<float, float>(meanII, meanII, r1, r1, ww, hh);
    apply(DIV, meanII, meanII, N);
    dbgsave(meanII, "/tmp/meanII.tif");

    array2D<float> varI(ww, hh);
    apply(SUBMUL, varI, meanI, meanI, meanII);
    dbgsave(varI, "/tmp/varI.tif");
    
    array2D<float> a(ww, hh);
    apply(DIVEPSILON, a, covIp, varI);
    boxblur<float, float>(a, a, r1, r1, ww, hh);
    apply(DIV, a, a, N);
    dbgsave(a, "/tmp/a.tif");

    array2D<float> b(ww, hh);
    apply(SUBMUL, b, a, meanI, meanp);
    boxblur<float, float>(b, b, r1, r1, ww, hh);
    apply(DIV, b, b, N);
    dbgsave(b, "/tmp/b.tif");
    
    array2D<float> aa(W, H);
    rescaleBilinear(a, aa, multithread);
    dbgsave(aa, "/tmp/aa.tif");
    
    array2D<float> bb(W, H);
    rescaleBilinear(b, bb, multithread);
    dbgsave(bb, "/tmp/bb.tif");

    apply(ADDMUL, dst, aa, guide, bb);
}

} // namespace rtengine
