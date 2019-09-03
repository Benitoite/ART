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
// extracted and datapted from ImProcFunctions::rgbProc (improcfun.cc) of
// RawTherapee

#include "improcfun.h"
#include "curves.h"
#include "color.h"

namespace rtengine {

namespace {

void shadowToneCurve(const LUTf &shtonecurve, Imagefloat *rgb, bool multiThread)
{
    float **rtemp = rgb->r.ptrs;
    float **gtemp = rgb->g.ptrs;
    float **btemp = rgb->b.ptrs;
    int W = rgb->getWidth();
    int H = rgb->getHeight();

#if defined( __SSE2__ ) && defined( __x86_64__ )
    vfloat cr = F2V(0.299f);
    vfloat cg = F2V(0.587f);
    vfloat cb = F2V(0.114f);
#endif

#ifdef _OPENMP
#   pragma omp parallel for if (multiThread)
#endif
    for (int y = 0; y < H; ++y) {
        int x = 0;
#if defined( __SSE2__ ) && defined( __x86_64__ )
        for (; x < W - 3; x += 4) {
            vfloat rv = LVF(rtemp[y][x]);
            vfloat gv = LVF(gtemp[y][x]);
            vfloat bv = LVF(btemp[y][x]);

            //shadow tone curve
            vfloat Yv = cr * rv + cg * gv + cb * bv;
            vfloat tonefactorv = shtonecurve(Yv);
            STVF(rtemp[y][x], rv * tonefactorv);
            STVF(gtemp[y][x], gv * tonefactorv);
            STVF(btemp[y][x], bv * tonefactorv);
        }
#endif

        for (; x < W; ++x) {
            float r = rtemp[y][x];
            float g = gtemp[y][x];
            float b = btemp[y][x];

            //shadow tone curve
            float Y = (0.299f * r + 0.587f * g + 0.114f * b);
            float tonefactor = shtonecurve[Y];
            rtemp[y][x] = rtemp[y][x] * tonefactor;
            gtemp[y][x] = gtemp[y][x] * tonefactor;
            btemp[y][x] = btemp[y][x] * tonefactor;
        }
    }
}


void highlightToneCurve(const LUTf &hltonecurve, Imagefloat *rgb, float exp_scale, float comp, float hlrange, bool multiThread)
{
    float **rtemp = rgb->r.ptrs;
    float **gtemp = rgb->g.ptrs;
    float **btemp = rgb->b.ptrs;
    int W = rgb->getWidth();
    int H = rgb->getHeight();

#if defined( __SSE2__ ) && defined( __x86_64__ )
    vfloat threev = F2V(3.f);
    vfloat maxvalfv = F2V(MAXVALF);
#endif

#ifdef _OPENMP
#   pragma omp parallel for if (multiThread)
#endif
    for (int y = 0; y < H; ++y) {
        int x = 0;
#if defined( __SSE2__ ) && defined( __x86_64__ )
        for (; x < W - 3; x += 4) {

            vfloat rv = LVF(rtemp[y][x]);
            vfloat gv = LVF(gtemp[y][x]);
            vfloat bv = LVF(btemp[y][x]);

            //TODO: proper treatment of out-of-gamut colors
            //float tonefactor = hltonecurve[(0.299f*r+0.587f*g+0.114f*b)];
            vmask maxMask = vmaskf_ge(vmaxf(rv, vmaxf(gv, bv)), maxvalfv);

            if (_mm_movemask_ps((vfloat)maxMask)) {
                for (int k = 0; k < 4; ++k) {
                    float r = rtemp[y][x + k];
                    float g = gtemp[y][x + k];
                    float b = btemp[y][x + k];
                    float tonefactor = ((r < MAXVALF ? hltonecurve[r] : CurveFactory::hlcurve(exp_scale, comp, hlrange, r)) +
                                        (g < MAXVALF ? hltonecurve[g] : CurveFactory::hlcurve(exp_scale, comp, hlrange, g)) +
                                        (b < MAXVALF ? hltonecurve[b] : CurveFactory::hlcurve(exp_scale, comp, hlrange, b))) / 3.0;

                    // note: tonefactor includes exposure scaling, that is here exposure slider and highlight compression takes place
                    rtemp[y][x + k] = r * tonefactor;
                    gtemp[y][x + k] = g * tonefactor;
                    btemp[y][x + k] = b * tonefactor;
                }
            } else {
                vfloat tonefactorv = (hltonecurve.cb(rv) + hltonecurve.cb(gv) + hltonecurve.cb(bv)) / threev;
                // note: tonefactor includes exposure scaling, that is here exposure slider and highlight compression takes place
                STVF(rtemp[y][x], rv * tonefactorv);
                STVF(gtemp[y][x], gv * tonefactorv);
                STVF(btemp[y][x], bv * tonefactorv);
            }
        }

#endif

        for (; x < W; ++x) {
            float r = rtemp[y][x];
            float g = gtemp[y][x];
            float b = btemp[y][x];

            //TODO: proper treatment of out-of-gamut colors
            //float tonefactor = hltonecurve[(0.299f*r+0.587f*g+0.114f*b)];
            float tonefactor = ((r < MAXVALF ? hltonecurve[r] : CurveFactory::hlcurve(exp_scale, comp, hlrange, r)) +
                                (g < MAXVALF ? hltonecurve[g] : CurveFactory::hlcurve(exp_scale, comp, hlrange, g)) +
                                (b < MAXVALF ? hltonecurve[b] : CurveFactory::hlcurve(exp_scale, comp, hlrange, b))) / 3.0;

            // note: tonefactor includes exposure scaling, that is here exposure slider and highlight compression takes place
            rtemp[y][x] = r * tonefactor;
            gtemp[y][x] = g * tonefactor;
            btemp[y][x] = b * tonefactor;
        }
    }
}

} // namespace


void ImProcFunctions::exposure(Imagefloat *img)
{
    if (!params->exposure.enabled) {
        return;
    }
    
    img->setMode(Imagefloat::Mode::RGB, multiThread);

    LUTf hltonecurve(65536);
    LUTf shtonecurve(65536);

    const double expcomp = params->exposure.expcomp;
    const int hlcompr = params->exposure.hlcompr;
    const int hlcomprthresh = params->exposure.hlcomprthresh;

    {
        const int black = params->exposure.black;
        const int shcompr = params->exposure.shcompr;
        
        LUTf tonecurve(65536);
        LUTu vhist16(65536), histToneCurve(256);
        ToneCurve customToneCurve1, customToneCurve2;
        
        CurveFactory::complexCurve(expcomp, black / 65535.0,
                                   hlcompr, hlcomprthresh,
                                   shcompr, 0, 0, 
                                   { DCT_Linear }, { DCT_Linear },
                                   vhist16, hltonecurve, shtonecurve, tonecurve,
                                   histToneCurve, customToneCurve1,
                                   customToneCurve2, scale);
    }

    const float exp_scale = pow(2.0, expcomp);
    const float comp = (max(0.0, expcomp) + 1.0) * hlcompr / 100.0;
    const float shoulder = ((65536.0 / max(1.0f, exp_scale)) * (hlcomprthresh / 200.0)) + 0.1;
    const float hlrange = 65536.0 - shoulder;

    highlightToneCurve(hltonecurve, img, exp_scale, comp, hlrange, multiThread);
    shadowToneCurve(shtonecurve, img, multiThread);    
}

} // namespace rtengine
