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

#include "improcfun.h"
#include "guidedfilter.h"
#include "rt_math.h"
#include "rt_algo.h"
#include "curves.h"
#include <iostream>
#include <queue>

extern Options options;

namespace rtengine {

namespace {

void guided_smoothing(array2D<float> &R, array2D<float> &G, array2D<float> &B, const DenoiseParams &pp, const TMatrix &ws, const TMatrix &iws, double scale, bool multithread)
{
    if (pp.guidedRadius > 0) {
        const int W = R.width();
        const int H = R.height();
        int r = max(int(pp.guidedRadius / scale), 1);
        const float epsilon = std::pow(2.f, 2 * (pp.guidedEpsilon - 10.f)) / 2.f;
        array2D<float> iR(W, H, R, 0);
        array2D<float> iG(W, H, G, 0);
        array2D<float> iB(W, H, B, 0);
        
        for (int i = 0; i < pp.guidedIterations; ++i) {
            guidedFilter(R, R, R, r, epsilon, multithread);
            guidedFilter(G, G, G, r, epsilon, multithread);
            guidedFilter(B, B, B, r, epsilon, multithread);
        }

        float l_blend = LIM01(pp.guidedLumaBlend / 100.f);
        float ab_blend = LIM01(pp.guidedChromaBlend / 100.f);

#ifdef _OPENMP
        #pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                float rr = R[y][x];
                float gg = G[y][x];
                float bb = B[y][x];
                float ir = iR[y][x];
                float ig = iG[y][x];
                float ib = iB[y][x];
                float Y = Color::rgbLuminance(rr, gg, bb, ws);
                float iY = Color::rgbLuminance(ir, ig, ib, ws);
                float oY = intp(l_blend, Y, iY);
                rr = intp(ab_blend, rr - Y, ir - iY); 
                gg = intp(ab_blend, gg - Y, ig - iY); 
                bb = intp(ab_blend, bb - Y, ib - iY);
                R[y][x] = oY + rr;
                G[y][x] = oY + gg;
                B[y][x] = oY + bb;
            }
        }
    }
}


// void guided_decomposition(array2D<float> &R, array2D<float> &G, array2D<float> &B, const GuidedFilterParams &gf, double scale, bool multithread)
// {
//     if (gf.decompRadius > 0) {
//         LUTf curve(65536);
//         {
//             DiagonalCurve baseCurve1(gf.decompBaseCurve1);
//             DiagonalCurve baseCurve2(gf.decompBaseCurve2);
//             constexpr float base = 100.f;
//             for (int i = 0; i < 65536; ++i) {
//                 float x = float(i)/65535.f;
//                 x = std::log(x * (base - 1.0f) + 1.0f) / std::log(base);
//                 float y = baseCurve1.getVal(x);
//                 y = (std::pow(base, y) - 1.0f) / (base - 1.0f);
//                 x = y;
//                 y = baseCurve2.getVal(x);
//                 curve[i] = y * 65535.f;
//             }
//         }

//         const int W = R.width();
//         const int H = R.height();
//         array2D<float> tmp(W, H);
        
//         const int r = max(int(gf.decompRadius / scale), 1);
//         const float boost = gf.decompDetailBoost >= 0 ? 1.f + gf.decompDetailBoost : -1.f/gf.decompDetailBoost;
//         const float epsilon = std::pow(2.f, 2 * (gf.decompEpsilon - 10.f)) / 2.f;
        
//         const auto apply =
//             [&](array2D<float> &chan) -> void
//             {
//                 guidedFilter(chan, chan, tmp, r, epsilon, multithread);

// #ifdef _OPENMP
//                 #pragma omp parallel for if (multithread)
// #endif
//                 for (int y = 0; y < H; ++y) {
//                     for (int x = 0; x < W; ++x) {
//                         float base = tmp[y][x];
//                         float detail = chan[y][x] - base;
//                         base *= 65535.f;
//                         curves::setLutVal(curve, base);
//                         chan[y][x] = max(base / 65535.f + boost * detail, 1e-5f);
//                     }
//                 }
//             };
        
//         apply(R);
//         apply(G);
//         apply(B);
//     }
// }

} // namespace


void ImProcFunctions::guidedSmoothing(Imagefloat *rgb)
{
    if (!params->denoise.smoothingEnabled || params->denoise.smoothingMethod != procparams::DenoiseParams::SmoothingMethod::GUIDED) {
        return;
    }

    rgb->normalizeFloatTo1();

    const int W = rgb->getWidth();
    const int H = rgb->getHeight();
    array2D<float> R(W, H, rgb->r.ptrs, ARRAY2D_BYREFERENCE);
    array2D<float> G(W, H, rgb->g.ptrs, ARRAY2D_BYREFERENCE);
    array2D<float> B(W, H, rgb->b.ptrs, ARRAY2D_BYREFERENCE);

    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
    TMatrix iws = ICCStore::getInstance()->workingSpaceInverseMatrix(params->icm.workingProfile);

    guided_smoothing(R, G, B, params->denoise, ws, iws, scale, multiThread);
    rgb->normalizeFloatTo65535();
}

} // namespace rtengine
