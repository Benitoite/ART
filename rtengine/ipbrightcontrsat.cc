/* -*- C++ -*-
 *
 *  This file is part of ART.
 *
 *  Copyright 2019 Alberto Griggio <alberto.griggio@gmail.com>
 *
 *  ART is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ART is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ART.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "improcfun.h"
#include "curves.h"
#include "color.h"
#include "rt_math.h"

namespace rtengine {

namespace {

float apply_vibrance(float x, float vib)
{
    static const float noise = pow_F(2.f, -16.f);
    float ax = std::abs(x / 65535.f);
    if (ax > noise) {
        return SGN(x) * pow_F(ax, vib) * 65535.f;
    } else {
        return x;
    }
}

} // namespace


void ImProcFunctions::brightnessContrastSaturation(Imagefloat *rgb)
{
    rgb->setMode(Imagefloat::Mode::RGB, multiThread);
    
    LUTf curve(65536);
    int bright = params->brightContrSat.enabled ? params->brightContrSat.brightness : 0;
    int contr = params->brightContrSat.enabled ? params->brightContrSat.contrast : 0;
    if (bright || contr) {
        LUTf curve1(65536);
        LUTf curve2(65536);
        LUTu dummy;
        LUTu hist16(65536);
        ToneCurve customToneCurve1, customToneCurve2;

        if (contr) {
            ImProcFunctions ipf(params, multiThread);
            ipf.firstAnalysis(rgb, *params, hist16);
        }
        CurveFactory::complexCurve(0, 0, 0, 0, 0, bright, contr,
                                   { DCT_Linear }, { DCT_Linear },
                                   hist16, curve1, curve2, curve, dummy,
                                   customToneCurve1, customToneCurve2, max(scale, 1.0));
    }

    const int W = rgb->getWidth();
    const int H = rgb->getHeight();

    if (!params->exposure.enabled || params->exposure.clampOOG) {
#ifdef _OPENMP
#       pragma omp parallel for if (multiThread)
#endif
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                float &r = rgb->r(i, j);
                float &g = rgb->g(i, j);
                float &b = rgb->b(i, j);
                if (OOG(r) || OOG(g) || OOG(b)) {
                    Color::filmlike_clip(&r, &g, &b);
                }
            }
        }
    }

    if (bright || contr) {
#ifdef _OPENMP
#       pragma omp parallel for if (multiThread)
#endif
        for (int i = 0; i < H; ++i) {
            int j = 0;
#ifdef __SSE2__
            vfloat tmpr;
            vfloat tmpg;
            vfloat tmpb;
            for (; j < W - 3; j += 4) {
                //brightness/contrast
                tmpr = curve(LVFU(rgb->r(i, j)));
                tmpg = curve(LVFU(rgb->g(i, j)));
                tmpb = curve(LVFU(rgb->b(i, j)));
                for (int k = 0; k < 4; ++k) {
                    setUnlessOOG(rgb->r(i, j+k), rgb->g(i, j+k), rgb->b(i, j+k), tmpr[k], tmpg[k], tmpb[k]);
                }
            }
#endif
            for (; j < W; ++j) {
                //brightness/contrast
                setUnlessOOG(rgb->r(i, j), rgb->g(i, j), rgb->b(i, j), curve[rgb->r(i, j)], curve[rgb->g(i, j)], curve[rgb->b(i, j)]);
            }
        }
    }

    if (params->brightContrSat.enabled &&
        (params->brightContrSat.saturation || params->brightContrSat.vibrance)) {
        const float saturation = 1.f + params->brightContrSat.saturation / 100.f;
        const float vibrance = 1.f - params->brightContrSat.vibrance / 1000.f;
        TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
        const float noise = pow_F(2.f, -16.f);
        const bool vib = params->brightContrSat.vibrance;
        
#ifdef _OPENMP
#       pragma omp parallel for if (multiThread)
#endif
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                float &r = rgb->r(i, j);
                float &g = rgb->g(i, j);
                float &b = rgb->b(i, j);
                float l = Color::rgbLuminance(r, g, b, ws);
                float rl = r - l;
                float gl = g - l;
                float bl = b - l;
                if (vib) {
                    rl = apply_vibrance(rl, vibrance);
                    gl = apply_vibrance(gl, vibrance);
                    bl = apply_vibrance(bl, vibrance);
                    assert(rl == rl);
                    assert(gl == gl);
                    assert(bl == bl);
                }
                r = max(l + saturation * rl, noise);
                g = max(l + saturation * gl, noise);
                b = max(l + saturation * bl, noise);
            }
        }
    }
}

} // namespace rtengine
