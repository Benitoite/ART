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
#include "linalgebra.h"

#include <fstream>

namespace rtengine {

namespace {

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


void apply_tc(Imagefloat *rgb, const ToneCurve &tc, ToneCurveParams::TcMode curveMode, const Glib::ustring &working_profile, const Glib::ustring &outprofile, int perceptual_strength, float whitept, bool multithread)
{
    const int W = rgb->getWidth();
    const int H = rgb->getHeight();
    
    if (curveMode == ToneCurveParams::TcMode::PERCEPTUAL) {
        const PerceptualToneCurve &c = static_cast<const PerceptualToneCurve&>(tc);
        PerceptualToneCurveState state;
        c.initApplyState(state, working_profile);
        state.strength = LIM01(float(perceptual_strength) / 100.f);

#ifdef _OPENMP
        #pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < H; ++y) {
            c.BatchApply(0, W, rgb->r.ptrs[y], rgb->g.ptrs[y], rgb->b.ptrs[y], state);
        }
    } else if (curveMode == ToneCurveParams::TcMode::STD) {
        const StandardToneCurve &c = static_cast<const StandardToneCurve &>(tc);
        apply(c, rgb, W, H, multithread);
    } else if (curveMode == ToneCurveParams::TcMode::WEIGHTEDSTD) {
        const WeightedStdToneCurve &c = static_cast<const WeightedStdToneCurve &>(tc);
        apply(c, rgb, W, H, multithread);
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
    } else if (curveMode == ToneCurveParams::TcMode::NEUTRAL) {
        const NeutralToneCurve &c = static_cast<const NeutralToneCurve &>(tc);
        NeutralToneCurve::ApplyState state(working_profile, outprofile);

#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < H; ++y) {
            c.BatchApply(0, W, rgb->r.ptrs[y], rgb->g.ptrs[y], rgb->b.ptrs[y], state);
        }
    }
}


class ContrastCurve: public Curve {
public:
    ContrastCurve(double a, double b, double w): a_(a), b_(b), w_(w) {}
    void getVal(const std::vector<double>& t, std::vector<double>& res) const {}
    bool isIdentity () const { return false; }
    
    double getVal(double x) const
    {
        double res = lin2log(std::pow(x/w_, a_), b_)*w_;
        return res;
    }

private:
    double a_;
    double b_;
    double w_;
};


void filmlike_clip(Imagefloat *rgb, float whitept, bool multithread)
{
    const int W = rgb->getWidth();
    const int H = rgb->getHeight();
    const float Lmax = 65535.f * whitept;

#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int i = 0; i < H; ++i) {
        for (int j = 0; j < W; ++j) {
            float &r = rgb->r(i, j);
            float &g = rgb->g(i, j);
            float &b = rgb->b(i, j);
            Color::filmlike_clip(&r, &g, &b, Lmax);
        }
    }
}


void legacy_contrast(Imagefloat *rgb, const ImProcData &im, int contrast, const Glib::ustring &working_profile, float whitept)
{
    if (contrast) {
        ToneCurve tc;
        auto &curve = tc.lutToneCurve;
        curve(65536);

        tc.Set(DiagonalCurve({DCT_Empty}));
        
        LUTu hist16(65536);
        ImProcFunctions ipf(im.params, im.multiThread);
        ipf.firstAnalysis(rgb, *im.params, hist16);

        CurveFactory::contrastCurve(contrast, hist16, curve, max(im.scale, 1.0));
        apply_tc(rgb, tc, ToneCurveParams::TcMode::STD, working_profile, im.params->icm.outputProfile, 100, whitept, im.multiThread);
    }
}


std::unique_ptr<Curve> get_contrast_curve(Imagefloat *rgb, const ImProcData &im, int contrast, float whitept)
{
    std::unique_ptr<Curve> ccurve;
    
    if (contrast) {
        const double pivot = (im.params->logenc.enabled ? im.params->logenc.targetGray / 100.0 : 0.18) / whitept;
        const double c = std::pow(std::abs(contrast) / 100.0, 1.5) * 16.0;
        const double b = contrast > 0 ? (1 + c) : 1.0 / (1 + c);
        const double a = std::log((std::exp(std::log(b) * pivot) - 1) / (b - 1)) / std::log(pivot);

        ccurve.reset(new ContrastCurve(a, b, whitept));
    }
    return ccurve;
}


float expand_range(float whitept, float x)
{
    if (whitept <= 1.001f) {
        return x;
    }

    float f = (pow_F(whitept, x) - 1) / (whitept - 1);
    float g = rtengine::intp(SQR(x)*x, f * whitept, x);
    return g;
}


void satcurve_lut(const FlatCurve &curve, LUTf &sat, float whitept)
{
    sat(65536, LUT_CLIP_BELOW);
    sat[0] = curve.getVal(0) * 2.f;
    for (int i = 1; i < 65536; ++i) {
        float x = Color::gamma2curve[i] / 65535.f;
        float v = curve.getVal(x);
        sat[i] = v * 2.f;
    }
}


class SatCurveRemap {
public:
    SatCurveRemap(float whitept):
        whitept_(whitept),
        remapcurve_({
            DCT_CatmullRom,
                0.0, 0.0,
                0.4, 0.4,
                whitept, 1.0
            })
    {}

    float operator()(float x) const
    {
        return Color::gamma2curve[LIM01((whitept_ == 1.f) ? x : remapcurve_.getVal(x)) * 65535.f] / 65535.f;
    }

private:
    float whitept_;
    DiagonalCurve remapcurve_;
};


void apply_satcurve(Imagefloat *rgb, const FlatCurve &curve, const DiagonalCurve &curve2, const Glib::ustring &working_profile, float whitept, bool multithread)
{
    LUTf sat;
    const bool use_lut = (whitept == 1.f);
    if (use_lut) {
        satcurve_lut(curve, sat, whitept);
    }

    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(working_profile);
    TMatrix iws = ICCStore::getInstance()->workingSpaceInverseMatrix(working_profile);

    const bool use_c2 = !curve2.isIdentity();
    SatCurveRemap remap(whitept);

#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < rgb->getHeight(); ++y) {
        float X, Y, Z;
        float Jz, az, bz;
        float cz, hz;
        for (int x = 0; x < rgb->getWidth(); ++x) {
            float &R = rgb->r(y, x);
            float &G = rgb->g(y, x);
            float &B = rgb->b(y, x);
            Color::rgbxyz(R/65535.f, G/65535.f, B/65535.f, X, Y, Z, ws);
            Color::xyz2jzazbz(X, Y, Z, Jz, az, bz);
            Color::jzazbz2jzch(az, bz, cz, hz);
            if (use_c2) {
                cz = curve2.getVal(cz * 50.f) / 50.f;
            }
            float s = use_lut ? sat[Y * 65535.f] : curve.getVal(remap(Y)) * 2.f;
            cz *= s;
            Color::jzczhz2rgb(Jz, cz, hz, R, G, B, iws);
            R *= 65535.f;
            G *= 65535.f;
            B *= 65535.f;
        }
    }
}


void fill_satcurve_pipette(Imagefloat *rgb, EditUniqueID editID, PlanarWhateverData<float>* editWhatever, const Glib::ustring &working_profile, float whitept, bool multithread)
{
    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(working_profile);

    if (editID == EUID_ToneCurveSaturation) {
        SatCurveRemap remap(whitept);

#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < rgb->getHeight(); ++y) {
            for (int x = 0; x < rgb->getWidth(); ++x) {
                float r = rgb->r(y, x), g = rgb->g(y, x), b = rgb->b(y, x);
                float Y = Color::rgbLuminance(r, g, b, ws);
                float s = remap(Y / 65535.f);
                editWhatever->v(y, x) = LIM01(s);
            }
        }
    } else if (editID == EUID_ToneCurveSaturation2) {
#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < rgb->getHeight(); ++y) {
            for (int x = 0; x < rgb->getWidth(); ++x) {
                float r = rgb->r(y, x), g = rgb->g(y, x), b = rgb->b(y, x);
                float Jz, cz, hz;
                Color::rgb2jzczhz(r / 65535.f, g / 65535.f, b / 65535.f, Jz, cz, hz, ws);
                float v = 0.f;
                v = Jz > 1e-7f ? cz * 50 : 0.f;
                editWhatever->v(y, x) = LIM01(v);
            }
        }
    }
}


void update_tone_curve_histogram(Imagefloat *img, LUTu &hist, const Glib::ustring &profile, bool multithread)
{
    hist.clear();
    const int compression = log2(65536 / hist.getSize());

    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(profile);

#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < img->getHeight(); ++y) {
        for (int x = 0; x < img->getWidth(); ++x) {
            float r = CLIP(img->r(y, x));
            float g = CLIP(img->g(y, x));
            float b = CLIP(img->b(y, x));

            int y = CLIP<int>(Color::gamma2curve[Color::rgbLuminance(r, g, b, ws)]);//max(r, g, b)]);
            hist[y >> compression]++;
        }
    }

    // we make this log encoded
    int n = hist.getSize();
    float f = float(n);
    for (int i = 0; i < n; ++i) {
        hist[i] = xlin2log(float(hist[i]) / f, 2.f) * f;
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


class DoubleCurve: public Curve {
public:
    DoubleCurve(const Curve &c1, const Curve &c2):
        c1_(c1), c2_(c2) {}

    double getVal(double t) const override
    {
        return c2_.getVal(c1_.getVal(t));
    }
    
    void getVal(const std::vector<double>& t, std::vector<double>& res) const override
    {
        c1_.getVal(t, res);
        c2_.getVal(res, res);
    }

    bool isIdentity() const override
    {
        return c1_.isIdentity() && c2_.isIdentity();
    }
    
private:
    const Curve &c1_;
    const Curve &c2_;
};


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
    } else if ((editID == EUID_ToneCurveSaturation || editID == EUID_ToneCurveSaturation2) && pipetteBuffer->getDataProvider()->getCurrSubscriber()->getPipetteBufferType() == BT_SINGLEPLANE_FLOAT) {
        editWhatever = pipetteBuffer->getSinglePlaneBuffer();
    }

    if (params->toneCurve.enabled) {
        img->setMode(Imagefloat::Mode::RGB, multiThread);


        const float whitept = params->toneCurve.hasWhitePoint() ? params->toneCurve.whitePoint : 1.f;

        const bool single_curve = params->toneCurve.curveMode == params->toneCurve.curveMode2;
        
        ImProcData im(params, scale, multiThread);
        if (!(single_curve && params->toneCurve.curveMode == ToneCurveParams::TcMode::NEUTRAL)) {
            filmlike_clip(img, whitept, im.multiThread);
        }

        std::unique_ptr<Curve> ccurve;
        if (params->toneCurve.contrastLegacyMode) {
            legacy_contrast(img, im, params->toneCurve.contrast, params->icm.workingProfile, whitept);
        } else {
            ccurve = get_contrast_curve(img, im, params->toneCurve.contrast, whitept);
        }

        const auto expand =
            [whitept](double x) -> double { return expand_range(whitept, x); };
        
        const auto adjust =
            [whitept,&expand](std::vector<double> c) -> std::vector<double>
            {
                std::map<double, double> m;
                DiagonalCurveType tp = DiagonalCurveType(c[0]);
                bool add_c = (tp == DCT_CatmullRom || tp == DCT_Spline);
                DiagonalCurve curve(c);
                for (int i = 0; i < 25; ++i) {
                    double x = double(i)/100.0;
                    double v = Color::gammatab_srgb[x * 65535.0] / 65535.0;
                    double y = curve.getVal(v);
                    y = Color::igammatab_srgb[y * 65535.0] / 65535.0;
                    m[expand(x)] = expand(y);
                }
                for (int i = 25, j = 2; i < 100; ) {
                    double x = double(i)/100.0;
                    double v = Color::gammatab_srgb[x * 65535.0] / 65535.0;
                    double y = curve.getVal(v);
                    y = Color::igammatab_srgb[y * 65535.0] / 65535.0;
                    m[expand(x)] = expand(y);
                    i += j;
                    j *= 2;
                }
                if (add_c) {
                    for (size_t i = 0; i < (c.size()-1)/2; ++i) {
                        double x = c[2*i+1];
                        double v = Color::gammatab_srgb[x * 65535.0] / 65535.0;
                        double y = curve.getVal(v);
                        y = Color::igammatab_srgb[y * 65535.0] / 65535.0;
                        m[expand(x)] = expand(y);
                    }
                } else {
                    m[expand(1.0)] = expand(curve.getVal(1.0));
                }
                c = { DCT_CatmullRom };
                for (auto &p : m) {
                    c.push_back(p.first);
                    c.push_back(p.second);
                }
#if 0
                auto name = "/tmp/curve" + std::to_string(int(whitept)) + ".txt";
                std::ofstream out(name.c_str());
                for (size_t i = 1; i < c.size(); i += 2) {
                    out << c[i] << " " << c[i+1] << "\n";
                }
#endif // if 0
                return c;
            };

        ToneCurve tc;
        DiagonalCurve tcurve2(adjust(params->toneCurve.curve2), CURVES_MIN_POLY_POINTS / max(int(scale), 1));
        DiagonalCurve tcurve1(adjust(params->toneCurve.curve), CURVES_MIN_POLY_POINTS / max(int(scale), 1));
        DoubleCurve dcurve(tcurve1, tcurve2);
        std::unique_ptr<Curve> dccurve;
        Curve *tcurve = &dcurve;
        if (ccurve) {
            dccurve.reset(new DoubleCurve(*ccurve, dcurve));
            tcurve = dccurve.get();
        }

        if (single_curve) {
            if (editImgFloat && (editID == EUID_ToneCurve1 || editID == EUID_ToneCurve2)) {
                fill_pipette(img, editImgFloat, multiThread);
            }
            tc.Set(*tcurve, whitept);
            apply_tc(img, tc, params->toneCurve.curveMode, params->icm.workingProfile, params->icm.outputProfile, params->toneCurve.perceptualStrength, whitept, multiThread);
        } else {
            if (ccurve) {
                tc.Set(*ccurve, whitept);
                apply_tc(img, tc, params->toneCurve.curveMode, params->icm.workingProfile, params->icm.outputProfile, 100, whitept, multiThread);
            }
            
            if (editImgFloat && editID == EUID_ToneCurve1) {
                fill_pipette(img, editImgFloat, multiThread);
            }
        
            if (!tcurve1.isIdentity()) {
                tc.Set(tcurve1, whitept);
                apply_tc(img, tc, params->toneCurve.curveMode, params->icm.workingProfile, params->icm.outputProfile, params->toneCurve.perceptualStrength, whitept, multiThread);
            }

            if (editImgFloat && editID == EUID_ToneCurve2) {
                fill_pipette(img, editImgFloat, multiThread);
            }

            if (!tcurve2.isIdentity()) {
                tc.Set(tcurve2, whitept);
                apply_tc(img, tc, params->toneCurve.curveMode2, params->icm.workingProfile, params->icm.outputProfile, params->toneCurve.perceptualStrength, whitept, multiThread);
            }
        }

        if (editWhatever) {
            fill_satcurve_pipette(img, editID, editWhatever, params->icm.workingProfile, whitept, multiThread);
        }

        auto satcurve_pts = params->toneCurve.saturation;
        const FlatCurve satlcurve(satcurve_pts, false, CURVES_MIN_POLY_POINTS / max(int(scale), 1));
        const DiagonalCurve satccurve(params->toneCurve.saturation2);
        if (!satlcurve.isIdentity() || !satccurve.isIdentity()) {
            apply_satcurve(img, satlcurve, satccurve, params->icm.workingProfile, whitept, multiThread);
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
