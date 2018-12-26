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
#include <iostream>
#include <queue>

extern Options options;

namespace rtengine {

namespace {

void guided_smoothing(array2D<float> &R, array2D<float> &G, array2D<float> &B, const TMatrix &ws, bool luminance, int radius, float epsilon, int strength, double scale, bool multithread)
{
    if (radius > 0 && strength > 0) {
        const int W = R.width();
        const int H = R.height();
        int r = max(int(radius / scale), 1);
        array2D<float> iR(W, H, R, 0);
        array2D<float> iG(W, H, G, 0);
        array2D<float> iB(W, H, B, 0);
        
        guidedFilter(R, R, R, r, epsilon, multithread);
        guidedFilter(G, G, G, r, epsilon, multithread);
        guidedFilter(B, B, B, r, epsilon, multithread);

        const float blend = LIM01(float(strength) / 100.f);

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
                const float bf = blend;
                float oY = luminance ? intp(bf, Y, iY) : iY;
                rr = luminance ? ir - iY : intp(bf, rr - Y, ir - iY); 
                gg = luminance ? ig - iY : intp(bf, gg - Y, ig - iY); 
                bb = luminance ? ib - iY : intp(bf, bb - Y, ib - iY);
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

    int r = max(int(params->denoise.guidedChromaRadius / scale), 1);
    guided_smoothing(R, G, B, ws, false, params->denoise.guidedChromaRadius, 0.1f / float(min(r, 10)), params->denoise.guidedChromaStrength, scale, multiThread);
    guided_smoothing(R, G, B, ws, true, params->denoise.guidedLumaRadius, 0.001f, params->denoise.guidedLumaStrength, scale, multiThread);
    
    rgb->normalizeFloatTo65535();
}


void ImProcFunctions::guidedSmoothing(LabImage *lab, int offset_x, int offset_y, int full_width, int full_height)
{
    PlanarWhateverData<float> *editWhatever = nullptr;
    EditUniqueID eid = pipetteBuffer ? pipetteBuffer->getEditID() : EUID_None;

    if ((eid == EUID_LabMasks_H3 || eid == EUID_LabMasks_C3 || eid == EUID_LabMasks_L3) && pipetteBuffer->getDataProvider()->getCurrSubscriber()->getPipetteBufferType() == BT_SINGLEPLANE_FLOAT) {
        editWhatever = pipetteBuffer->getSinglePlaneBuffer();
    }
    
    if (params->smoothing.enabled) {
        if (editWhatever) {
            LabMasksEditID id = static_cast<LabMasksEditID>(int(eid) - EUID_LabMasks_H3);
            fillPipetteLabMasks(lab, editWhatever, id, multiThread);
        }
        
        int n = params->smoothing.regions.size();
        int show_mask_idx = params->smoothing.showMask;
        if (show_mask_idx >= n) {
            show_mask_idx = -1;
        }
        std::vector<array2D<float>> Lmask(n);
        std::vector<array2D<float>> abmask(n);
        if (!generateLabMasks(lab, params->smoothing.labmasks, offset_x, offset_y, full_width, full_height, scale, multiThread, show_mask_idx, &Lmask, &abmask)) {
            return; // show mask is active, nothing more to do
        }

        Imagefloat rgb(lab->W, lab->H);

        TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
        
        for (int i = 0; i < n; ++i) {
            lab2rgb(*lab, rgb, params->icm.workingProfile);
            rgb.normalizeFloatTo1();
            
            auto &r = params->smoothing.regions[i];

            const int W = rgb.getWidth();
            const int H = rgb.getHeight();
            array2D<float> R(W, H, rgb.r.ptrs, ARRAY2D_BYREFERENCE);
            array2D<float> G(W, H, rgb.g.ptrs, ARRAY2D_BYREFERENCE);
            array2D<float> B(W, H, rgb.b.ptrs, ARRAY2D_BYREFERENCE);
            
            guided_smoothing(R, G, B, ws, false, r.chromaRadius, r.chromaStrength, scale, multiThread);
            guided_smoothing(R, G, B, ws, true, r.lumaRadius, r.lumaStrength, scale, multiThread);
            const auto &Lblend = Lmask[i];
            const auto &abblend = abmask[i];

            rgb.normalizeFloatTo65535();
            rgb2lab(rgb, *lab, params->icm.workingProfile);
            
#ifdef _OPENMP
#           pragma omp parallel for if (multiThread)
#endif
            for (int y = 0; y < lab->H; ++y) {
                for (int x = 0; x < lab->W; ++x) {
                    float l = lab->L[y][x];
                    float ll, aa, bb;
                    Color::rgb2lab(R[y][x], G[y][x], B[y][x], ll, aa, bb, ws);
                    lab->L[y][x] = intp(Lblend[y][x], ll, lab->L[y][x]);
                    lab->a[y][x] = intp(abblend[y][x], aa, lab->a[y][x]);
                    lab->b[y][x] = intp(abblend[y][x], bb, lab->b[y][x]);
                }
            }
        }
    } else if (editWhatever) {
        editWhatever->fill(0.f);
    }
}

} // namespace rtengine
