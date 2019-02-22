/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright 2018 Alberto Griggio <alberto.griggio@gmail.com>
 *
 *  
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

// adapted from the tone equalizer of darktable
/*
    This file is part of darktable,
    copyright (c) 2018 Aurelien Pierre.

    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/    

#include "improcfun.h"
#include "gauss.h"
#include "sleef.c"
#include "opthelper.h"
#include "guidedfilter.h"

namespace rtengine {

namespace {

void tone_eq(array2D<float> &R, array2D<float> &G, array2D<float> &B, const ToneEqualizerParams &pp, const Glib::ustring &workingProfile, double scale, bool multithread)
{
    const int W = R.width();
    const int H = R.height();
    array2D<float> Y(W, H);

    const int r = max(int(30 / scale), 1);
    const float epsilon = 0.001f;

    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(workingProfile);

#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            Y[y][x] = Color::rgbLuminance(R[y][x], G[y][x], B[y][x], ws);
        }
    }
    
    rtengine::guidedFilter(Y, Y, Y, r, epsilon, multithread);

    const auto log2 =
        [](float x) -> float
        {
            static const float l2 = xlogf(2);
            return xlogf(x) / l2;
        };

    const auto exp2 =
        [](float x) -> float
        {
            return pow_F(2.f, x);
        };

    const auto GAUSS = [](float b, float x) -> float
                       {
                           return xexpf((-SQR(x - b) / 4.0f));
                       };

     // Build the luma channels: band-pass filters with gaussian windows of
     // std 2 EV, spaced by 2 EV
    const float centers[12] = {
        -18.0f, -16.0f, -14.0f, -12.0f, -10.0f, -8.0f, -6.0f,
        -4.0f, -2.0f, 0.0f, 2.0f, 4.0f
    };

    const auto conv = [&](int v) -> float
                      {
                          return exp2(float(v) / 100.f * 2.f);
                      };

    const float factors[12] = {
        conv(pp.bands[0]), // -18 EV
        conv(pp.bands[0]), // -16 EV
        conv(pp.bands[0]), // -14 EV
        conv(pp.bands[0]), // -12 EV
        conv(pp.bands[0]), // -10 EV
        conv(pp.bands[0]), //  -8 EV
        conv(pp.bands[1]), //  -6 EV
        conv(pp.bands[2]), //  -4 EV
        conv(pp.bands[3]), //  -2 EV
        conv(pp.bands[4]), //   0 EV
        conv(pp.bands[4]), //   2 EV
        conv(pp.bands[4])  //   4 EV
    };

    // For every pixel luminance, the sum of the gaussian masks
    // evenly spaced by 2 EV with 2 EV std should be this constant
    float w_sum = 0.f;
    for (int i = 0; i < 12; ++i) {
        w_sum += GAUSS(centers[i], 0.f);
    }

    const auto process_pixel =
        [&](float y) -> float
        {
            // get the luminance of the pixel - just average channels
            const float luma = max(log2(max(y, 0.f)), -18.0f);

            // build the correction as the sum of the contribution of each
            // luminance channel to current pixel
            float correction = 0.0f;
            for (int c = 0; c < 12; ++c) {
                correction += GAUSS(centers[c], luma) * factors[c];
            }
            correction /= w_sum;

            return correction;
        };
        
#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            float oY = Color::rgbLuminance(R[y][x], G[y][x], B[y][x], ws);
            float cY = Y[y][x];
            float corr = process_pixel(cY);
            float corr2 = corr + pow_F(corr / 2.f, 4.f);
            float dY = oY - cY;
            if (oY > 0.f) {
                R[y][x] = R[y][x] / oY * (cY * corr + dY * corr2);
                G[y][x] = G[y][x] / oY * (cY * corr + dY * corr2);
                B[y][x] = B[y][x] / oY * (cY * corr + dY * corr2);
            }
        }
    }
}

} // namespace


void ImProcFunctions::toneEqualizer(Imagefloat *rgb)
{
    if (!params->toneEqualizer.enabled) {
        return;
    }

    rgb->normalizeFloatTo1();

    array2D<float> R(rgb->getWidth(), rgb->getHeight(), rgb->r.ptrs, ARRAY2D_BYREFERENCE);
    array2D<float> G(rgb->getWidth(), rgb->getHeight(), rgb->g.ptrs, ARRAY2D_BYREFERENCE);
    array2D<float> B(rgb->getWidth(), rgb->getHeight(), rgb->b.ptrs, ARRAY2D_BYREFERENCE);

    tone_eq(R, G, B, params->toneEqualizer, params->icm.workingProfile, scale, multiThread);
    
    rgb->normalizeFloatTo65535();
}

} // namespace rtengine
