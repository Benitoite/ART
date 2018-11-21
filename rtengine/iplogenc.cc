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
    //  base^source_gray - 1 - base * target_gray + target_gray = 0
    //
    // use a bisection method (maybe later change to Netwon)

    const auto f =
        [=](float x) -> float
        {
            return std::pow(x, source_gray) - 1 - target_gray * x + target_gray;
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


// basic log encoding taken from ACESutil.Lin_to_Log2, from
// https://github.com/ampas/aces-dev
// (as seen on pixls.us)
void log_encode(Imagefloat *rgb, const ProcParams *params, bool multithread)
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

#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < rgb->getHeight(); ++y) {
        for (int x = 0; x < rgb->getWidth(); ++x) {
            float r = rgb->r(y, x);
            float g = rgb->g(y, x);
            float b = rgb->b(y, x);
            r = apply(r);
            g = apply(g);
            b = apply(b);
            if (saturation_enabled) {
                float l = Color::rgbLuminance(r, g, b, ws);
                r = max(l + saturation * (r - l), noise);
                g = max(l + saturation * (g - l), noise);
                b = max(l + saturation * (b - l), noise);
            }

            if (OOG(r) || OOG(g) || OOG(b)) {
                Color::filmlike_clip(&r, &g, &b);
            }
            
            rgb->r(y, x) = r;
            rgb->g(y, x) = g;
            rgb->b(y, x) = b;
        }
    }
}

} // namespace


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

    for (int y = 0, h = fh / 10; y < h; ++y) {
        for (int x = 0, w = fw / 10; x < w; ++x) {
            float l = Color::rgbLuminance(img.r(y, x), img.g(y, x), img.b(y, x), ws);
            if (l > 0.f) {
                data.push_back(l);
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

    constexpr float dr_headroom = 1.f;
    constexpr float black_headroom = 0.5f;
    
    if (vmax > vmin) {
        const float log2 = xlogf(2.f);
        const float noise = pow_F(2.f, -16.f);
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


void ImProcFunctions::logEncoding(LabImage *lab)
{
    if (!params->logenc.enabled) {
        return;
    }

    Imagefloat working(lab->W, lab->H);
    lab2rgb(*lab, working, params->icm.workingProfile);

    log_encode(&working, params, multiThread);

    if (dcpProf && dcpApplyState) {
        for (int y = 0; y < lab->H; ++y) {
            float *r = working.r.ptrs[y];
            float *g = working.g.ptrs[y];
            float *b = working.b.ptrs[y];
            dcpProf->step2ApplyTile(r, g, b, lab->W, 1, 1, *dcpApplyState);
        }
    }    

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
