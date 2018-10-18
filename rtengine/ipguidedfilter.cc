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

inline void rgb2lab(float R, float G, float B, float &l, float &a, float &b, const TMatrix &ws)
{
    float x, y, z;
    Color::rgbxyz(R, G, B, x, y, z, ws);
    Color::XYZ2Lab(x, y, z, l, a, b);
}


inline void lab2rgb(float l, float a, float b, float &R, float &G, float &B, const TMatrix &iws)
{
    float x, y, z;
    Color::Lab2XYZ(l, a, b, x, y, z);
    Color::xyz2rgb(x, y, z, R, G, B, iws);
}


void decompose(LabImage *lab, array2D<float> &R, array2D<float> &G, array2D<float> &B, const TMatrix &iws, bool multithread)
{
    const int W = lab->W;
    const int H = lab->H;
    
#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            float r, g, b;
            lab2rgb(lab->L[y][x], lab->a[y][x], lab->b[y][x], r, g, b, iws);
            R[y][x] = r / 65535.f;
            G[y][x] = g / 65535.f;
            B[y][x] = b / 65535.f;
        }
    }
}


void recompose(const array2D<float> &R, const array2D<float> &G, const array2D<float> &B, LabImage *lab, const TMatrix &ws, bool multithread)
{
    const int W = lab->W;
    const int H = lab->H;
    
#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            rgb2lab(R[y][x] * 65535.f, G[y][x] * 65535.f, B[y][x] * 65535.f, lab->L[y][x], lab->a[y][x], lab->b[y][x], ws);
        }
    }
}


void guided_smoothing(LabImage *lab, array2D<float> &R, array2D<float> &G, array2D<float> &B, const GuidedFilterParams &gf, const TMatrix &ws, const TMatrix &iws, double scale, bool multithread)
{
    if (gf.smoothingRadius > 0) {
        const int W = lab->W;
        const int H = lab->H;
        int r = max(int(gf.smoothingRadius / scale), 1);
        int s = LIM(r / 2, 1, 4);
        for (int i = 0; i < gf.smoothingIterations; ++i) {
            guidedFilter(R, R, R, r, gf.smoothingEpsilon, multithread, s);
            guidedFilter(G, G, G, r, gf.smoothingEpsilon, multithread, s);
            guidedFilter(B, B, B, r, gf.smoothingEpsilon, multithread, s);
        }

        float l_blend = LIM01(gf.smoothingLumaBlend / 100.f);
        float ab_blend = LIM01(gf.smoothingChromaBlend / 100.f);

#ifdef _OPENMP
        #pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                float l, a, b;
                rgb2lab(R[y][x] * 65535.f, G[y][x] * 65535.f, B[y][x] * 65535.f, l, a, b, ws);
                float rr, gg, bb;
                lab2rgb(intp(l_blend, l, lab->L[y][x]), intp(ab_blend, a, lab->a[y][x]), intp(ab_blend, b, lab->b[y][x]), rr, gg, bb, iws);
                R[y][x] = rr / 65535.f;
                G[y][x] = gg / 65535.f;
                B[y][x] = bb / 65535.f;
            }
        }
    }
}


void guided_decomposition(LabImage *lab, array2D<float> &R, array2D<float> &G, array2D<float> &B, const GuidedFilterParams &gf, double scale, bool multithread)
{
    if (gf.decompRadius > 0) {
        StandardToneCurve curve;
        DiagonalCurve baseCurve(gf.decompBaseCurve);
        curve.Set(baseCurve, Color::sRGBGammaCurve);

        const int W = lab->W;
        const int H = lab->H;
        array2D<float> tmp(W, H);
        
        const int r = max(int(gf.decompRadius / scale), 1);
        const int s = LIM(r / 2, 1, 4);
        const float boost = gf.decompDetailBoost >= 0 ? 1.f + gf.decompDetailBoost : -1.f/gf.decompDetailBoost;

        const auto apply =
            [&](array2D<float> &chan) -> void
            {
                guidedFilter(chan, chan, tmp, r, gf.decompEpsilon, multithread, s);

#ifdef _OPENMP
                #pragma omp parallel for if (multithread)
#endif
                for (int y = 0; y < H; ++y) {
                    for (int x = 0; x < W; ++x) {
                        float base = tmp[y][x];
                        float detail = chan[y][x] - base;
                        base *= 65535.f;
                        curves::setLutVal(curve.lutToneCurve, base);
                        chan[y][x] = base / 65535.f + boost * detail;
                    }
                }
            };
        
        apply(R);
        apply(G);
        apply(B);
    }
}

} // namespace


void ImProcFunctions::guidedFilter(LabImage *lab)
{
    if (!params->guidedfilter.enabled) {
        return;
    }

    array2D<float> R(lab->W, lab->H);
    array2D<float> G(lab->W, lab->H);
    array2D<float> B(lab->W, lab->H);

    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
    TMatrix iws = ICCStore::getInstance()->workingSpaceInverseMatrix(params->icm.workingProfile);

    decompose(lab, R, G, B, iws, multiThread);
    guided_smoothing(lab, R, G, B, params->guidedfilter, ws, iws, scale, multiThread);
    guided_decomposition(lab, R, G, B, params->guidedfilter, scale, multiThread);
    recompose(R, G, B, lab, ws, multiThread);
}

} // namespace rtengine
