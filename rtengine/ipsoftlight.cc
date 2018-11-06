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
#include "curves.h"

namespace rtengine {

namespace {

inline float sl(float blend, float x)
{
    if (!OOG(x)) {
        const float orig = 1.f - blend;
        float v = Color::gamma_srgb(x) / MAXVALF;
        // Pegtop's formula from
        // https://en.wikipedia.org/wiki/Blend_modes#Soft_Light
        float v2 = v * v;
        float v22 = v2 * 2.f;
        v = v2 + v22 - v22 * v;
        x = blend * Color::igamma_srgb(v * MAXVALF) + orig * x;
    }
    return x;
}


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


} // namespace


void ImProcFunctions::softLight(LabImage *lab)
{
    const bool sl_enabled = params->softlight.enabled && params->softlight.strength > 0;
    const bool logenc_enabled = params->logenc.enabled && (!params->toneCurve.curve2.empty() && params->toneCurve.curve2[0] != DCT_Linear);
                                                           
    if (!sl_enabled && !logenc_enabled) {
        return;
    }

    Imagefloat working(lab->W, lab->H);
    lab2rgb(*lab, working, params->icm.workingProfile);
    
    if (logenc_enabled) {
        ToneCurve tc;
        const DiagonalCurve tcurve(params->toneCurve.curve2, CURVES_MIN_POLY_POINTS / max(int(scale), 1));

        if (!tcurve.isIdentity()) {
            tc.Set(tcurve, Color::sRGBGammaCurve);
            apply_tc(&working, tc, params->toneCurve.curveMode2, params->icm.workingProfile, multiThread);
        }
    }
    
    if (sl_enabled) {
        const float blend = params->softlight.strength / 100.f;

#ifdef _OPENMP
        #pragma omp parallel for
#endif
        for (int y = 0; y < working.getHeight(); ++y) {
            for (int x = 0; x < working.getWidth(); ++x) {
                working.r(y, x) = sl(blend, working.r(y, x));
                working.g(y, x) = sl(blend, working.g(y, x));
                working.b(y, x) = sl(blend, working.b(y, x));
            }
        }
    }

    rgb2lab(working, *lab, params->icm.workingProfile);
}

} // namespace rtengine
