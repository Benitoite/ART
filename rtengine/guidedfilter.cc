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

#include "guidedfilter.h"
#include "boxblur.h"
#include "rescale.h"
#include "imagefloat.h"

namespace rtengine {

#if 0
#  define DEBUG_DUMP(arr)                                                 \
    do {                                                                \
        Imagefloat im(arr.width(), arr.height());                      \
        const char *out = "/tmp/" #arr ".tif";                     \
        for (int y = 0; y < im.getHeight(); ++y) {                      \
            for (int x = 0; x < im.getWidth(); ++x) {                   \
                im.r(y, x) = im.g(y, x) = im.b(y, x) = arr[y][x] * 65535.f; \
            }                                                           \
        }                                                               \
        im.saveTIFF(out, 16);                                           \
    } while (false)
#else
#  define DEBUG_DUMP(arr)
#endif


void guidedFilter(const array2D<float> &guide, const array2D<float> &src, array2D<float> &dst, int r, float epsilon, bool multithread, int subsampling)
{
    const int W = src.width();
    const int H = src.height();

    enum Op { MUL, DIV, ADD, SUB, ADDMUL, SUBMUL };

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

    DEBUG_DUMP(guide);
    DEBUG_DUMP(src);
    DEBUG_DUMP(guide1);
    DEBUG_DUMP(src1);

    float r1 = float(r) / subsampling;

    array2D<float> meanI(ww, hh);
    boxblur<float, float>(guide1, meanI, r1, r1, ww, hh);
    DEBUG_DUMP(meanI);


    array2D<float> meanp(ww, hh);
    boxblur<float, float>(src1, meanp, r1, r1, ww, hh);
    DEBUG_DUMP(meanp);

    array2D<float> meanIp(ww, hh);
    apply(MUL, meanIp, guide1, src1);
    boxblur<float, float>(meanIp, meanIp, r1, r1, ww, hh);
    DEBUG_DUMP(meanIp);

    array2D<float> covIp(ww, hh);
    apply(SUBMUL, covIp, meanI, meanp, meanIp);
    DEBUG_DUMP(covIp);

    array2D<float> meanII(ww, hh);
    apply(MUL, meanII, guide1, guide1);
    boxblur<float, float>(meanII, meanII, r1, r1, ww, hh);
    DEBUG_DUMP(meanII);

    array2D<float> varI(ww, hh);
    apply(SUBMUL, varI, meanI, meanI, meanII);
    DEBUG_DUMP(varI);
    
    array2D<float> a(ww, hh);
    apply(DIV, a, covIp, varI);
    boxblur<float, float>(a, a, r1, r1, ww, hh);
    DEBUG_DUMP(a);

    array2D<float> b(ww, hh);
    apply(SUBMUL, b, a, meanI, meanp);
    boxblur<float, float>(b, b, r1, r1, ww, hh);
    DEBUG_DUMP(b);
    
    array2D<float> aa(W, H);
    rescaleBilinear(a, aa, multithread);
    DEBUG_DUMP(aa);
    
    array2D<float> bb(W, H);
    rescaleBilinear(b, bb, multithread);
    DEBUG_DUMP(bb);

    apply(ADDMUL, dst, aa, guide, bb);
}

} // namespace rtengine
