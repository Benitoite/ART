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
#include "labmasks.h"
#include "sleef.c"
#include <iostream>
#include <queue>

extern Options options;

namespace rtengine {

namespace {

enum class Channel {
    L,
    C,
    LC
};


void guided_smoothing(array2D<float> &R, array2D<float> &G, array2D<float> &B, const TMatrix &ws, const TMatrix &iws, Channel chan, int radius, float epsilon, int strength, double scale, bool multithread)
{
    const auto rgb2xyY =
        [&](float r, float g, float b, float &x, float &y, float &Y) -> void
        {
            float X, Z;
            Color::rgbxyz(r, g, b, X, Y, Z, ws);
            Color::XYZ_to_xyY(X, Y, Z, x, y);
        };

    const auto xyY2rgb =
        [&](float x, float y, float Y, float &r, float &g, float &b) -> void
        {
            float X, Z;
            Color::xyY_to_XYZ(x, y, Y, X, Z);
            Color::xyz2rgb(X, Y, Z, r, g, b, iws);
        };
        
    if (radius > 0 && strength > 0) {
        const int W = R.width();
        const int H = R.height();
        int r = max(int(radius / scale), 1);

        array2D<float> iR(W, H, R, 0);
        array2D<float> iG(W, H, G, 0);
        array2D<float> iB(W, H, B, 0);

        rtengine::guidedFilterLog(10.f, R, r, epsilon, multithread);
        rtengine::guidedFilterLog(10.f, G, r, epsilon, multithread);
        rtengine::guidedFilterLog(10.f, B, r, epsilon, multithread);

        const float blend = LIM01(float(strength) / 100.f);
        const bool rgb = (chan == Channel::LC);
        const bool luminance = (chan == Channel::L);

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
                const float bf = blend;
                if (rgb) {
                    R[y][x] = intp(bf, rr, ir);
                    G[y][x] = intp(bf, gg, ig);
                    B[y][x] = intp(bf, bb, ib);
                } else {
                    // float Y = Color::rgbLuminance(rr, gg, bb, ws);
                    // float iY = Color::rgbLuminance(ir, ig, ib, ws);
                    // float oY = luminance ? intp(bf, Y, iY) : iY;
                    // rr = luminance ? ir - iY : intp(bf, rr - Y, ir - iY); 
                    // gg = luminance ? ig - iY : intp(bf, gg - Y, ig - iY); 
                    // bb = luminance ? ib - iY : intp(bf, bb - Y, ib - iY);
                    // R[y][x] = oY + rr;
                    // G[y][x] = oY + gg;
                    // B[y][x] = oY + bb;

                    float ix, iy, iY;
                    float ox, oy, oY;
                    rgb2xyY(ir, ig, ib, ix, iy, iY);
                    rgb2xyY(rr, gg, bb, ox, oy, oY);
                    if (luminance) {
                        oY = intp(bf, oY, iY);
                        ox = ix;
                        oy = iy;
                    } else {
                        ox = intp(bf, ox, ix);
                        oy = intp(bf, oy, iy);
                        oY = iY;
                    }
                    xyY2rgb(ox, oy, oY, R[y][x], G[y][x], B[y][x]);
                }
            }
        }
    }
}

} // namespace


void ImProcFunctions::denoiseGuidedSmoothing(Imagefloat *rgb)
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

    int r = max(int(params->denoise.guidedChromaRadius / scale), 1);
    const float c_eps = 0.1f / float(min(r, 10));
    const float l_eps = 5e-4f;

    guided_smoothing(R, G, B, ws, iws, Channel::C, params->denoise.guidedChromaRadius, c_eps, params->denoise.guidedChromaStrength, scale, multiThread);
    guided_smoothing(R, G, B, ws, iws, Channel::L, params->denoise.guidedLumaRadius, l_eps, params->denoise.guidedLumaStrength, scale, multiThread);
    
    rgb->normalizeFloatTo65535();
}


bool ImProcFunctions::guidedSmoothing(Imagefloat *rgb, int offset_x, int offset_y, int full_width, int full_height)
{
    PlanarWhateverData<float> *editWhatever = nullptr;
    EditUniqueID eid = pipetteBuffer ? pipetteBuffer->getEditID() : EUID_None;

    if ((eid == EUID_LabMasks_H3 || eid == EUID_LabMasks_C3 || eid == EUID_LabMasks_L3) && pipetteBuffer->getDataProvider()->getCurrSubscriber()->getPipetteBufferType() == BT_SINGLEPLANE_FLOAT) {
        editWhatever = pipetteBuffer->getSinglePlaneBuffer();
    }
    
    if (params->smoothing.enabled) {
        // LabImage tmplab(rgb->getWidth(), rgb->getHeight());
        // rgb2lab(*rgb, tmplab);
        // LabImage *lab = &tmplab;
        
        if (editWhatever) {
            LabMasksEditID id = static_cast<LabMasksEditID>(int(eid) - EUID_LabMasks_H3);
            fillPipetteLabMasks(rgb, editWhatever, id, multiThread);
        }
        
        int n = params->smoothing.regions.size();
        int show_mask_idx = params->smoothing.showMask;
        if (show_mask_idx >= n) {
            show_mask_idx = -1;
        }
        std::vector<array2D<float>> mask(n);
        if (!generateLabMasks(rgb, params->smoothing.labmasks, offset_x, offset_y, full_width, full_height, scale, multiThread, show_mask_idx, nullptr, &mask)) {
            // lab2rgb(*lab, *rgb);
            return true; // show mask is active, nothing more to do
        }

        const int W = rgb->getWidth();
        const int H = rgb->getHeight();

        Imagefloat working(W, H);
        rgb->normalizeFloatTo1();
        //rgb->assignColorSpace(params->icm.workingProfile);
        rgb->setMode(Imagefloat::Mode::RGB, multiThread);

        TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
        TMatrix iws = ICCStore::getInstance()->workingSpaceInverseMatrix(params->icm.workingProfile);
        
        for (int i = 0; i < n; ++i) {
            rgb->copyTo(&working);
            // lab2rgb(*lab, working, params->icm.workingProfile);
            // working.normalizeFloatTo1();
            
            auto &r = params->smoothing.regions[i];

            // const int W = working.getWidth();
            // const int H = working.getHeight();
            array2D<float> R(W, H, working.r.ptrs, ARRAY2D_BYREFERENCE);
            array2D<float> G(W, H, working.g.ptrs, ARRAY2D_BYREFERENCE);
            array2D<float> B(W, H, working.b.ptrs, ARRAY2D_BYREFERENCE);

            const float epsilon = 0.001f * std::pow(2, -r.epsilon);
            guided_smoothing(R, G, B, ws, iws, Channel(int(r.channel)), r.radius, epsilon, 100, scale, multiThread);
            
            const auto &blend = mask[i];
            
#ifdef _OPENMP
#           pragma omp parallel for if (multiThread)
#endif
            for (int y = 0; y < H; ++y) {
                for (int x = 0; x < W; ++x) {
                    // float ll, aa, bb;
                    // Color::rgb2lab(R[y][x] * 65535.f, G[y][x] * 65535.f, B[y][x] * 65535.f, ll, aa, bb, ws);
                    // lab->L[y][x] = intp(blend[y][x], ll, lab->L[y][x]);
                    // lab->a[y][x] = intp(blend[y][x], aa, lab->a[y][x]);
                    // lab->b[y][x] = intp(blend[y][x], bb, lab->b[y][x]);
                    rgb->r(y, x) = intp(blend[y][x], working.r(y, x), rgb->r(y, x));
                    rgb->g(y, x) = intp(blend[y][x], working.g(y, x), rgb->g(y, x));
                    rgb->b(y, x) = intp(blend[y][x], working.b(y, x), rgb->b(y, x));
                }
            }
        }

        rgb->normalizeFloatTo65535();
        // lab2rgb(*lab, *rgb);
    } else if (editWhatever) {
        editWhatever->fill(0.f);
    }

    return false;
}

} // namespace rtengine
