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
#include "clutstore.h"
#include "guidedfilter.h"

namespace rtengine {

extern const Settings *settings;

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


float find_gray(float source_gray, float target_gray)
{
    // find a base such that log2lin(base, source_gray) = target_gray
    // log2lin is (base^source_gray - 1) / (base - 1), so we solve
    //
    //  (base^source_gray - 1) / (base - 1) = target_gray, that is
    //
    //  base^source_gray - 1 - base * target_gray + target_gray = 0
    //
    // use a bisection method (maybe later change to Netwon)

    if (source_gray <= 0.f) {
        return 0.f;
    }

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
void log_encode(Imagefloat *rgb, const ProcParams *params, float scale, bool multithread)
{
    if (!params->logenc.enabled) {
        return;
    }

    const float gray = params->logenc.sourceGray / 100.f;
    const float shadows_range = params->logenc.blackEv;
    const float dynamic_range = params->logenc.whiteEv - params->logenc.blackEv;
    const float noise = pow_F(2.f, -16.f);
    const float log2 = xlogf(2.f);
    const bool brightness_enabled = params->brightContrSat.enabled && params->brightContrSat.brightness;
    const float brightness = 1.f + params->brightContrSat.brightness / 100.f;
    const bool contrast_enabled = params->brightContrSat.enabled && params->brightContrSat.contrast;
    const float contrast = 1.f + params->brightContrSat.contrast / 100.f;
    const bool saturation_enabled = params->brightContrSat.enabled && params->brightContrSat.saturation;
    const float saturation = 1.f + params->brightContrSat.saturation / 100.f;
    const float b = params->logenc.targetGray > 1 && params->logenc.targetGray < 100 && dynamic_range > 0 ? find_gray(std::abs(params->logenc.blackEv) / dynamic_range, params->logenc.targetGray / 100.f) : 0.f;
    const float linbase = max(b, 0.f);
    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);

    const auto apply =
        [=](float x, bool scale=true) -> float
        {
            if (scale) {
                x /= 65535.f;
            }
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
            if (linbase > 0.f) {
                x = xlog2lin(x, linbase);
            }
            if (scale) {
                return x * 65535.f;
            } else {
                return x;
            }
        };

    const auto norm =
        [&](float r, float g, float b) -> float
        {
            return Color::rgbLuminance(r, g, b, ws);

            // other possible alternatives (so far, luminance seems to work
            // fine though). See also
            // https://discuss.pixls.us/t/finding-a-norm-to-preserve-ratios-across-non-linear-operations
            //
            // MAX
            //return max(r, g, b);
            //
            // Euclidean
            //return std::sqrt(SQR(r) + SQR(g) + SQR(b));
            
            // weighted yellow power norm from https://youtu.be/Z0DS7cnAYPk
            // float rr = 1.22f * r / 65535.f;
            // float gg = 1.20f * g / 65535.f;
            // float bb = 0.58f * b / 65535.f;
            // float rr4 = SQR(rr) * SQR(rr);
            // float gg4 = SQR(gg) * SQR(gg);
            // float bb4 = SQR(bb) * SQR(bb);
            // float den = (rr4 + gg4 + bb4);
            // if (den > 0.f) {
            //     return 0.8374319f * ((rr4 * rr + gg4 * gg + bb4 * bb) / den) * 65535.f;
            // } else {
            //     return 0.f;
            // }
        };

    const int detail = float(max(params->logenc.detail, 0)) / scale + 0.5f;
    if (detail == 0) {
#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < rgb->getHeight(); ++y) {
            for (int x = 0; x < rgb->getWidth(); ++x) {
                float r = rgb->r(y, x);
                float g = rgb->g(y, x);
                float b = rgb->b(y, x);
                float m = norm(r, g, b);
                if (m > noise) {
                    float mm = apply(m);
                    float f = mm / m;
                    r *= f;
                    b *= f;
                    g *= f;
                    if (saturation_enabled) {
                        float l = Color::rgbLuminance(r, g, b, ws);
                        r = max(l + saturation * (r - l), noise);
                        g = max(l + saturation * (g - l), noise);
                        b = max(l + saturation * (b - l), noise);
                    }

                    // if (OOG(r) || OOG(g) || OOG(b)) {
                    //     Color::filmlike_clip(&r, &g, &b);
                    // }
                }
            
                assert(r == r);
                assert(g == g);
                assert(b == b);

                rgb->r(y, x) = r;
                rgb->g(y, x) = g;
                rgb->b(y, x) = b;
            }
        }
    } else {
        const int W = rgb->getWidth(), H = rgb->getHeight();
        array2D<float> tmp(W, H);
#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                tmp[y][x] = norm(rgb->r(y, x), rgb->g(y, x), rgb->b(y, x)) / 65535.f;
            }
        }
            
        const float epsilon = 0.01f;
        guidedFilter(tmp, tmp, tmp, detail, epsilon, multithread);

        // array2D<float> blend(W, H);
        // float contrast = 0.2f;
        // buildBlendMask(tmp, blend, W, H, contrast);
               
#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                float &r = rgb->r(y, x);
                float &g = rgb->g(y, x);
                float &b = rgb->b(y, x);
                float m = norm(r, g, b);
                float t = intp(0.33f// blend[y][x]
                               , m, tmp[y][x] * 65535.f);
                if (t > noise) {
                    float c = apply(t);
                    float f = c / t;
                    r *= f;
                    g *= f;
                    b *= f;
                    if (saturation_enabled) {
                        float l = Color::rgbLuminance(r, g, b, ws);
                        r = max(l + saturation * (r - l), noise);
                        g = max(l + saturation * (g - l), noise);
                        b = max(l + saturation * (b - l), noise);
                    }
                    // if (OOG(r) || OOG(g) || OOG(b)) {
                    //     Color::filmlike_clip(&r, &g, &b);
                    // }
                }
            }
        }

        // rgb->normalizeFloatTo65535();
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


void brightness_contrast_saturation(Imagefloat *rgb, const ProcParams *params, float scale, bool multithread)
{
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
            ImProcFunctions ipf(params, multithread);
            ipf.firstAnalysis(rgb, *params, hist16);
        }
        CurveFactory::complexCurve(0, 0, 0, 0, 0, bright, contr,
                                   { DCT_Linear }, { DCT_Linear },
                                   hist16, curve1, curve2, curve, dummy,
                                   customToneCurve1, customToneCurve2, max(scale, 1.f));
    }

    const int W = rgb->getWidth();
    const int H = rgb->getHeight();

    if (!params->exposure.enabled || params->exposure.clampOOG) {
#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
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
#       pragma omp parallel for if (multithread)
#endif
        for (int i = 0; i < H; ++i) {
            int j = 0;
#ifdef __SSE2__
            vfloat tmpr;
            vfloat tmpg;
            vfloat tmpb;
            for (; j < W - 3; j += 4) {
                //brightness/contrast
                STVF(tmpr[0], curve(LVF(rgb->r(i, j))));
                STVF(tmpg[0], curve(LVF(rgb->g(i, j))));
                STVF(tmpb[0], curve(LVF(rgb->b(i, j))));
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

    if (params->brightContrSat.enabled && params->brightContrSat.saturation) {
        const float saturation = 1.f + params->brightContrSat.saturation / 100.f;
        TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
        const float noise = pow_F(2.f, -16.f);
        
#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                float &r = rgb->r(i, j);
                float &g = rgb->g(i, j);
                float &b = rgb->b(i, j);
                float l = Color::rgbLuminance(r, g, b, ws);
                r = max(l + saturation * (r - l), noise);
                g = max(l + saturation * (g - l), noise);
                b = max(l + saturation * (b - l), noise);
            }
        }
    }
}

} // namespace


void ImProcFunctions::getAutoLog(ImageSource *imgsrc, LogEncodingParams &lparams)
{
    constexpr int SCALE = 10;
    int fw, fh, tr = TR_NONE;
    imgsrc->getFullSize(fw, fh, tr);
    PreviewProps pp(0, 0, fw, fh, SCALE);
    Imagefloat img(int(fw / SCALE + 0.5), int(fh / SCALE + 0.5));
    ProcParams neutral;
    neutral.exposure.enabled = true;
    neutral.exposure.clampOOG = false;
    imgsrc->getImage(imgsrc->getWB(), tr, &img, pp, neutral.exposure, neutral.raw);
    imgsrc->convertColorSpace(&img, params->icm, imgsrc->getWB());

    float vmin = RT_INFINITY;
    float vmax = -RT_INFINITY;
    const float ec = params->exposure.enabled ? std::pow(2.f, params->exposure.expcomp) : 1.f;

    constexpr float noise = 1e-5;

    for (int y = 0, h = fh / SCALE; y < h; ++y) {
        for (int x = 0, w = fw / SCALE; x < w; ++x) {
            float r = img.r(y, x), g = img.g(y, x), b = img.b(y, x);
            float m = max(0.f, r, g, b) / 65535.f * ec;
            if (m > noise) {
                float l = min(r, g, b) / 65535.f * ec;
                vmin = min(vmin, l > noise ? l : m);
                vmax = max(vmax, m);
            }
        }
    }

    if (vmax > vmin) {
        const float log2 = xlogf(2.f);
        float dynamic_range = -xlogf(vmin / vmax) / log2;
        if (settings->verbose) {
            std::cout << "AutoLog: min = " << vmin << ", max = " << vmax
                      << ", DR = " << dynamic_range << std::endl;
        }

        if (lparams.autogray) {
            double tot = 0.f;
            int n = 0;
            float gmax = std::min(vmax / 2.f, 0.25f);
            float gmin = std::max(vmin * std::pow(2.f, std::max((dynamic_range - 1.f) / 2.f, 1.f)), 0.05f);
            if (settings->verbose) {
                std::cout << "         gray boundaries: " << gmin << ", " << gmax << std::endl;
            }
            for (int y = 0, h = fh / SCALE; y < h; ++y) {
                for (int x = 0, w = fw / SCALE; x < w; ++x) {
                    float l = img.g(y, x) / 65535.f;
                    if (l >= gmin && l <= gmax) {
                        tot += l;
                        ++n;
                    }
                }
            }
            if (n > 0) {
                lparams.sourceGray = tot / n * 100.f;
                if (settings->verbose) {
                    std::cout << "         computed gray point from " << n << " samples: " << lparams.sourceGray << std::endl;
                }
            } else if (settings->verbose) {
                std::cout << "         no samples found in range, resorting to default gray point value" << std::endl;
                lparams.sourceGray = LogEncodingParams().sourceGray;
            }
        }
        
        float gray = float(lparams.sourceGray) / 100.f;
        lparams.whiteEv = xlogf(vmax / gray) / log2;
        lparams.blackEv = lparams.whiteEv - dynamic_range;
    }
}


void ImProcFunctions::logEncoding(LabImage *lab, LUTu *histToneCurve)
{
    Imagefloat working(lab->W, lab->H);
    lab2rgb(*lab, working, params->icm.workingProfile);

    if (params->logenc.enabled) {
        log_encode(&working, params, scale, multiThread);
    } else {
        brightness_contrast_saturation(&working, params, scale, multiThread);
    }

    if (dcpProf && dcpApplyState) {
        for (int y = 0; y < lab->H; ++y) {
            float *r = working.r.ptrs[y];
            float *g = working.g.ptrs[y];
            float *b = working.b.ptrs[y];
            dcpProf->step2ApplyTile(r, g, b, lab->W, 1, 1, *dcpApplyState);
        }
    }

//    shadowsHighlights(&working);
    
    HaldCLUTApplication hald_clut(params->filmSimulation.clutFilename, params->icm.workingProfile);
    if (params->filmSimulation.enabled) {
        constexpr int TS = 112;
        hald_clut.init(float(params->filmSimulation.strength)/100.f, TS);

        if (hald_clut) {
#ifdef _OPENMP
            #pragma omp parallel for if (multiThread)
#endif
            for (int y = 0; y < working.getHeight(); ++y) {
                for (int jj = 0; jj < working.getWidth(); jj += TS) {
                    int jstart = jj;
                    float *r = working.r.ptrs[y]+jstart;
                    float *g = working.g.ptrs[y]+jstart;
                    float *b = working.b.ptrs[y]+jstart;
                    int tW = min(jj + TS, working.getWidth());
                    hald_clut(r, g, b, 0, jstart, tW, 1);
                }
            }
        }
    }

    if (histToneCurve && *histToneCurve) {
        update_tone_curve_histogram(&working, *histToneCurve, params->icm.workingProfile, multiThread);
    }

    if (params->toneCurve.enabled) {
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
    }

    shadowsHighlights(&working);

    rgb2lab(working, *lab, params->icm.workingProfile);
}


} // namespace rtengine
