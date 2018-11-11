/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright 2018 Alberto Griggio <alberto.griggio@gmail.com>
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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "improcfun.h"
#include "sleef.c"
#include "imagesource.h"
#include "rt_algo.h"
#include "curves.h"

namespace rtengine {

namespace {

template <class Curve>
inline void apply_batch(const Curve &c, Imagefloat *rgb, int W, int H, bool multithread)
{
#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        c.BatchApply(0, W, rgb->r.ptrs[y], rgb->g.ptrs[y], rgb->b.ptrs[y]);
    }
}


template <class Curve>
inline void apply(const Curve &c, Imagefloat *rgb, int W, int H, bool multithread)
{
#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            c.Apply(rgb->r(y, x), rgb->g(y, x), rgb->b(y, x));
        }
    }
}


void apply_tc(Imagefloat *rgb, const ToneCurve &tc, ToneCurveParams::TcMode curveMode, const Glib::ustring &working_profile, bool multithread)
{
    const int W = rgb->getWidth();
    const int H = rgb->getHeight();
    
    if (curveMode == ToneCurveParams::TcMode::PERCEPTUAL) {
        const PerceptualToneCurve &c = static_cast<const PerceptualToneCurve&>(tc);
        PerceptualToneCurveState state;
        c.initApplyState(state, working_profile);

#ifdef _OPENMP
        #pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < H; ++y) {
            c.BatchApply(0, W, rgb->r.ptrs[y], rgb->g.ptrs[y], rgb->b.ptrs[y], state);
        }
    } else if (curveMode == ToneCurveParams::TcMode::STD) {
        const StandardToneCurve &c = static_cast<const StandardToneCurve &>(tc);
        apply_batch(c, rgb, W, H, multithread);
    } else if (curveMode == ToneCurveParams::TcMode::WEIGHTEDSTD) {
        const WeightedStdToneCurve &c = static_cast<const WeightedStdToneCurve &>(tc);
        apply_batch(c, rgb, W, H, multithread);
    } else if (curveMode == ToneCurveParams::TcMode::FILMLIKE) {
        const AdobeToneCurve &c = static_cast<const AdobeToneCurve &>(tc);
        apply(c, rgb, W, H, multithread);
    } else if (curveMode == ToneCurveParams::TcMode::SATANDVALBLENDING) {
        const SatAndValueBlendingToneCurve &c = static_cast<const SatAndValueBlendingToneCurve &>(tc);
        apply(c, rgb, W, H, multithread);
    } else if (curveMode == ToneCurveParams::TcMode::LUMINANCE) {
        const LuminanceToneCurve &c = static_cast<const LuminanceToneCurve &>(tc);
        apply(c, rgb, W, H, multithread);
    }
}


float find_brightness(float source_gray, float target_gray)
{
    // find a base such that log2lin(base, source_gray) = target_gray
    // log2lin is (base^source_gray - 1) / (base - 1), so we solve
    //
    //  (base^source_gray - 1) / (base - 1) = target_gray, that is
    //
    //  base^source_gray - 1 - base * target_gray - target_gray = 0
    //
    // use a bisection method (maybe later change to Netwon)

    const auto f =
        [=](float x) -> float
        {
            return pow_F(x, source_gray) - 1 - target_gray * x - target_gray;
        };

    // first find the interval we are interested in
    float lo = 1.f;
    while (f(lo) <= 0.f) {
        lo *= 2.f;
    }

    float hi = lo * 2.f;
    while (f(hi) >= 0.f) {
        hi *= 2.f;
    }

    if (std::isinf(hi)) {
        return 0.f;
    }

    // now search for a zero
    for (int iter = 0; iter < 100; ++iter) {
        float mid = lo + (hi - lo) / 2.f;
        float v = f(mid);
        if (std::abs(v) < 1e-4f || (hi - lo) / lo <= 1e-4f) {
            return mid;
        }
        if (v > 0.f) {
            lo = mid;
        } else {
            hi = mid;
        }
    }

    return 0.f; // not found
}


} // namespace


// taken from darktable (src/iop/profile_gamma.c)
/*
   copyright (c) 2009--2010 johannes hanika.
   copyright (c) 2014 LebedevRI.
   copyright (c) 2018 Aurélien Pierre.

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
void ImProcFunctions::logEncoding(float *r, float *g, float *b, int istart, int jstart, int tW, int tH, int TS)
{
    if (!params->logenc.enabled) {
        return;
    }
    
    const float gray = params->logenc.grayPoint / 100.f;
    const float shadows_range = params->logenc.blackEv;
    const float dynamic_range = params->logenc.whiteEv - params->logenc.blackEv;
    const float noise = pow_F(2.f, -16.f);
    const float log2 = xlogf(2.f);
    const bool brightness_enabled = params->toneCurve.brightness;
    const float brightness = 1.f + params->toneCurve.brightness / 100.f; //params->toneCurve.brightness <= 0 ? -params->toneCurve.brightness / 10.f : 1.f / (1.f + params->toneCurve.brightness);// / 10.f);
//    const float base = pow_F(2.f, brightness);
    const bool contrast_enabled = params->toneCurve.contrast;
    const float contrast = 1.f + params->toneCurve.contrast / 100.f; //params->toneCurve.contrast >= 0 ? 1.f + params->toneCurve.contrast / 100.f : 1.f/(1.f - params->toneCurve.contrast / 50.f);
    const bool saturation_enabled = params->toneCurve.saturation;
    const float saturation = 1.f + params->toneCurve.saturation / 100.f;
    const float norm = params->logenc.base > 0 ? pow_F(2.f, params->logenc.base) : 0;
    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);

    const auto apply =
        [=](float x) -> float
        {
            x /= 65535.f;
            x = max(x, noise);
            if (brightness_enabled) {
                x *= brightness;
            }
            if (contrast_enabled) {
                x = pow_F(x / gray, contrast) * gray;
            }
            x = max(x / gray, noise);
            x = max((xlogf(x)/log2 - shadows_range) / dynamic_range, noise);
            assert(x == x);
            if (norm > 0.f) {
                x = xlog2lin(x, norm);
            }
            // if (brightness_enabled) {
            //     x = xlog2lin(x, base);
            // }
            return x * 65535.f;
        };

    for (int i = istart, ti = 0; i < tH; i++, ti++) {
        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
            int idx = ti * TS + tj;
            r[idx] = apply(r[idx]);
            g[idx] = apply(g[idx]);
            b[idx] = apply(b[idx]);
            if (saturation_enabled) {
                float l = Color::rgbLuminance(r[idx], g[idx], b[idx], ws);
                r[idx] = l + saturation * (r[idx] - l);
                g[idx] = l + saturation * (g[idx] - l);
                b[idx] = l + saturation * (b[idx] - l);
            }
        }
    }
}


void ImProcFunctions::getAutoLog(ImageSource *imgsrc, LogEncodingParams &lparams)
{
    int fw, fh, tr = TR_NONE;
    imgsrc->getFullSize(fw, fh, tr);
    PreviewProps pp(0, 0, fw, fh, 10);
    Imagefloat img(int(fw / 10 + 0.5), int(fh / 10 + 0.5));
    imgsrc->getImage(imgsrc->getWB(), tr, &img, pp, params->toneCurve, params->raw);
    imgsrc->convertColorSpace(&img, params->icm, imgsrc->getWB());

    std::vector<float> data;
    data.reserve(fw/10 * fh/10 * 2);

    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
    double tot = 0.0;

    for (int y = 0, h = fh / 10; y < h; ++y) {
        for (int x = 0, w = fw / 10; x < w; ++x) {
            float l = Color::rgbLuminance(img.r(y, x), img.g(y, x), img.b(y, x), ws);
            if (l > 0.f) {
                data.push_back(l);
                tot += l;
            }
        }
    }

    std::sort(data.begin(), data.end());
    int n = data.size();
    if (!n) {
        return;
    }
    
    float vmin = data[0];
    float vmax = data[n-1];
    float vmid = tot / n;

    constexpr float dr_headroom = 1.f;
    constexpr float black_headroom = 0.5f;
    
    if (vmax > vmin) {
        const float log2 = xlogf(2.f);
        const float noise = pow_F(2.f, -16.f);
        lparams.grayPoint = int(vmid / 65535.f * 100.f);
        const float gray = float(lparams.grayPoint) / 100.f;
        vmin = max(vmin / vmax, noise);
        float lmin = xlogf(vmin) / log2;
        lparams.blackEv = xlogf(vmin/gray) / log2 - black_headroom;
        float dynamic_range = std::abs(lmin) + dr_headroom;
        lparams.whiteEv = dynamic_range + lparams.blackEv;
        float b = find_brightness(std::abs(lparams.blackEv) / dynamic_range, 0.18);
        if (b > 0.f) {
            lparams.base = std::log(b) / log2;
        } else {
            lparams.base = 0;
        }
    }
}


void ImProcFunctions::logEncodingCurves(LabImage *lab)
{
    const bool logenc_enabled = params->logenc.enabled && ((!params->toneCurve.curve.empty() && params->toneCurve.curve[0] != DCT_Linear) || (!params->toneCurve.curve2.empty() && params->toneCurve.curve2[0] != DCT_Linear));

    if (!logenc_enabled) {
        return;
    }

    Imagefloat working(lab->W, lab->H);
    lab2rgb(*lab, working, params->icm.workingProfile);

    ToneCurve tc;
    const DiagonalCurve tcurve1(params->toneCurve.curve, CURVES_MIN_POLY_POINTS / max(int(scale), 1));

    if (!tcurve1.isIdentity()) {
        tc.Set(tcurve1, Color::sRGBGammaCurve);
        apply_tc(&working, tc, params->toneCurve.curveMode, params->icm.workingProfile, multiThread);
    }

    const DiagonalCurve tcurve2(params->toneCurve.curve2, CURVES_MIN_POLY_POINTS / max(int(scale), 1));

    if (!tcurve2.isIdentity()) {
        tc.Set(tcurve2, Color::sRGBGammaCurve);
        apply_tc(&working, tc, params->toneCurve.curveMode2, params->icm.workingProfile, multiThread);
    }

    rgb2lab(working, *lab, params->icm.workingProfile);
}


} // namespace rtengine
