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

void guided_smoothing(array2D<float> &R, array2D<float> &G, array2D<float> &B, const TMatrix &ws, bool luminance, int radius, int strength, double scale, bool multithread)
{
    if (radius > 0 && strength > 0) {
        const int W = R.width();
        const int H = R.height();
        int r = max(int(radius / scale), 1);
        const float epsilon = luminance ? 0.001f : 0.25f;
        array2D<float> iR(W, H, R, 0);
        array2D<float> iG(W, H, G, 0);
        array2D<float> iB(W, H, B, 0);
        
        guidedFilter(R, R, R, r, epsilon, multithread);
        guidedFilter(G, G, G, r, epsilon, multithread);
        guidedFilter(B, B, B, r, epsilon, multithread);

        float blend = LIM01(float(strength) / 100.f);

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
                float oY = luminance ? intp(blend, Y, iY) : iY;
                rr = luminance ? ir - iY : intp(blend, rr - Y, ir - iY); 
                gg = luminance ? ig - iY : intp(blend, gg - Y, ig - iY); 
                bb = luminance ? ib - iY : intp(blend, bb - Y, ib - iY);
                R[y][x] = oY + rr;
                G[y][x] = oY + gg;
                B[y][x] = oY + bb;
            }
        }
    }
}

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

    guided_smoothing(R, G, B, ws, false, params->denoise.guidedChromaRadius, params->denoise.guidedChromaStrength, scale, multiThread);
    guided_smoothing(R, G, B, ws, true, params->denoise.guidedLumaRadius, params->denoise.guidedLumaStrength, scale, multiThread);
    
    rgb->normalizeFloatTo65535();
}

} // namespace rtengine
