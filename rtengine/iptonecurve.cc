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
#include "sleef.h"
#include "curves.h"

namespace rtengine {

namespace {

void apply_contrast(Imagefloat *rgb, const ImProcData &im, int contrast)
{
    const int W = rgb->getWidth();
    const int H = rgb->getHeight();
    
#ifdef _OPENMP
#   pragma omp parallel for if (im.multiThread)
#endif
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            float &r = rgb->r(i, j);
            float &g = rgb->g(i, j);
            float &b = rgb->b(i, j);
            Color::filmlike_clip(&r, &g, &b);
        }
    }
    
    if (!contrast) {
        return;
    }

    LUTf curve(65536);
    LUTf curve1(65536);
    LUTf curve2(65536);
    LUTu dummy;
    LUTu hist16(65536);
    ToneCurve customToneCurve1, customToneCurve2;

    ImProcFunctions ipf(im.params, im.multiThread);
    ipf.firstAnalysis(rgb, *im.params, hist16);
    CurveFactory::complexCurve(0, 0, 0, 0, 0, 0, contrast,
                               { DCT_Linear }, { DCT_Linear },
                               hist16, curve1, curve2, curve, dummy,
                               customToneCurve1, customToneCurve2,
                               max(im.scale, 1.0));
    

#ifdef _OPENMP
#   pragma omp parallel for if (im.multiThread)
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
            STVFU(rgb->r(i, j), tmpr);
            STVFU(rgb->g(i, j), tmpg);
            STVFU(rgb->b(i, j), tmpb);
        }
#endif
        for (; j < W; ++j) {
            //brightness/contrast
            rgb->r(i, j) = curve[rgb->r(i, j)];
            rgb->g(i, j) = curve[rgb->g(i, j)];
            rgb->b(i, j) = curve[rgb->b(i, j)];
        }
    }
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
        TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(working_profile);
        const LuminanceToneCurve &c = static_cast<const LuminanceToneCurve &>(tc);
//        apply(c, rgb, W, H, multithread);
#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                c.Apply(rgb->r(y, x), rgb->g(y, x), rgb->b(y, x), ws);
            }
        }
    }
}


void satcurve_lut(const FlatCurve &curve, LUTf &sat)
{
    sat(65536, LUT_CLIP_BELOW);
    sat[0] = curve.getVal(0) * 2.f;
    for (int i = 1; i < 65536; ++i) {
        float v = curve.getVal(Color::gamma2curve[i] / 65535.f);//pow_F(i / 65535.f, 1.f/2.2f));
        sat[i] = v * 2.f;
    }
}


void apply_satcurve(Imagefloat *rgb, const FlatCurve &curve, const Glib::ustring &working_profile, bool multithread)
{
    LUTf sat;
    satcurve_lut(curve, sat);

    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(working_profile);

#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < rgb->getHeight(); ++y) {
        for (int x = 0; x < rgb->getWidth(); ++x) {
            float r = rgb->r(y, x), g = rgb->g(y, x), b = rgb->b(y, x);
            float Y = Color::rgbLuminance(r, g, b, ws);
            float s = sat[Y];
            rgb->r(y, x) = Y + s * (r - Y);
            rgb->g(y, x) = Y + s * (g - Y);
            rgb->b(y, x) = Y + s * (b - Y);
        }
    }
}


void fill_satcurve_pipette(Imagefloat *rgb, PlanarWhateverData<float>* editWhatever, const Glib::ustring &working_profile, bool multithread)
{
    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(working_profile);

#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < rgb->getHeight(); ++y) {
        for (int x = 0; x < rgb->getWidth(); ++x) {
            float r = rgb->r(y, x), g = rgb->g(y, x), b = rgb->b(y, x);
            float Y = Color::rgbLuminance(r, g, b, ws);
            float s = Color::gamma2curve[Y];
            editWhatever->v(y, x) = LIM01(s / 65535.f);
        }
    }
}


void update_tone_curve_histogram(Imagefloat *img, LUTu &hist, const Glib::ustring &profile, bool multithread)
{
    hist.clear();
    const int compression = log2(65536 / hist.getSize());

#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < img->getHeight(); ++y) {
        for (int x = 0; x < img->getWidth(); ++x) {
            float r = CLIP(img->r(y, x));
            float g = CLIP(img->g(y, x));
            float b = CLIP(img->b(y, x));

            int y = CLIP<int>(Color::gamma2curve[max(r, g, b)]);
            hist[y >> compression]++;
        }
    }
}

void fill_pipette(Imagefloat *img, Imagefloat *pipette, bool multithread)
{
    const int W = img->getWidth();
    const int H = img->getHeight();
    
#ifdef _OPENMP
#    pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            pipette->r(y, x) = Color::gamma2curve[CLIP(img->r(y, x))] / 65535.f;
            pipette->g(y, x) = Color::gamma2curve[CLIP(img->g(y, x))] / 65535.f;
            pipette->b(y, x) = Color::gamma2curve[CLIP(img->b(y, x))] / 65535.f;
        }
    }
}

} // namespace


void ImProcFunctions::toneCurve(Imagefloat *img)
{
    if (histToneCurve && *histToneCurve) {
        img->setMode(Imagefloat::Mode::RGB, multiThread);
        update_tone_curve_histogram(img, *histToneCurve, params->icm.workingProfile, multiThread);
    }

    Imagefloat *editImgFloat = nullptr;
    PlanarWhateverData<float> *editWhatever = nullptr;
    EditUniqueID editID = pipetteBuffer ? pipetteBuffer->getEditID() : EUID_None;

    if ((editID == EUID_ToneCurve1 || editID == EUID_ToneCurve2) && pipetteBuffer->getDataProvider()->getCurrSubscriber()->getPipetteBufferType() == BT_IMAGEFLOAT) {
        editImgFloat = pipetteBuffer->getImgFloatBuffer();
    } else if (editID == EUID_ToneCurveSaturation && pipetteBuffer->getDataProvider()->getCurrSubscriber()->getPipetteBufferType() == BT_SINGLEPLANE_FLOAT) {
        editWhatever = pipetteBuffer->getSinglePlaneBuffer();
    }

    if (params->toneCurve.enabled) {
        img->setMode(Imagefloat::Mode::RGB, multiThread);

        ImProcData im(params, scale, multiThread);
        apply_contrast(img, im, params->toneCurve.contrast);

        if (editImgFloat && editID == EUID_ToneCurve1) {
            fill_pipette(img, editImgFloat, multiThread);
        }
        
        ToneCurve tc;
        const DiagonalCurve tcurve1(params->toneCurve.curve, CURVES_MIN_POLY_POINTS / max(int(scale), 1));

        if (!tcurve1.isIdentity()) {
            tc.Set(tcurve1, Color::sRGBGammaCurve);
            apply_tc(img, tc, params->toneCurve.curveMode, params->icm.workingProfile, multiThread);
        }

        if (editImgFloat && editID == EUID_ToneCurve2) {
            fill_pipette(img, editImgFloat, multiThread);
        }

        const DiagonalCurve tcurve2(params->toneCurve.curve2, CURVES_MIN_POLY_POINTS / max(int(scale), 1));

        if (!tcurve2.isIdentity()) {
            tc.Set(tcurve2, Color::sRGBGammaCurve);
            apply_tc(img, tc, params->toneCurve.curveMode2, params->icm.workingProfile, multiThread);
        }

        if (editWhatever) {
            fill_satcurve_pipette(img, editWhatever, params->icm.workingProfile, multiThread);
        }

        const FlatCurve satcurve(params->toneCurve.saturation, false, CURVES_MIN_POLY_POINTS / max(int(scale), 1));
        if (!satcurve.isIdentity()) {
            apply_satcurve(img, satcurve, params->icm.workingProfile, multiThread);
        }
    } else if (editImgFloat) {
        const int W = img->getWidth();
        const int H = img->getHeight();

#ifdef _OPENMP
#       pragma omp parallel for if (multiThread)
#endif
        for (int y = 0; y < H; ++y) {
            std::fill(editImgFloat->r(y), editImgFloat->r(y)+W, 0.f);
            std::fill(editImgFloat->g(y), editImgFloat->g(y)+W, 0.f);
            std::fill(editImgFloat->b(y), editImgFloat->b(y)+W, 0.f);
        }
    } else if (editWhatever) {
        editWhatever->fill(0.f);
    }
}

} // namespace rtengine
