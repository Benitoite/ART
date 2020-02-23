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

#include "labmasks.h"
#include "guidedfilter.h"
#include "sleef.h"
#include "coord.h"
#include "gauss.h"
#include "color.h"
#include "curves.h"
#include "iccstore.h"
#include "rt_algo.h"
#include "opthelper.h"

namespace rtengine {

using procparams::AreaMask;
using procparams::Mask;
using procparams::DrawnMask;

namespace {

#ifdef __SSE2__
void fastlin2log(float *x, float factor, float base, int w)
{
    float baseLog = 1.f / xlogf(base);
    vfloat baseLogv = F2V(baseLog);
    factor = factor * (base - 1.f);
    vfloat factorv = F2V(factor);
    vfloat onev = F2V(1.f);
    int i = 0;
    for (; i < w - 3; i += 4) {
        STVFU(x[i], xlogf(LVFU(x[i]) * factorv + onev) * baseLogv);
    }
    for (; i < w; ++i) {
        x[i] = xlogf(x[i] * factor + 1.f) * baseLog;
    }
}
#endif


bool generate_area_mask(int ox, int oy, int width, int height, const array2D<float> &guide, const AreaMask &areaMask, float scale, bool multithread, array2D<float> &mask)
{
    bool enabled = areaMask.enabled;
    float blur = areaMask.blur / scale;
    
    if (!enabled || areaMask.shapes.empty() || (areaMask.isTrivial() && blur <= 0.f)) {
        return false;
    }

    float w2 = float(width) / 2;
    float h2 = float(height) / 2;

    Coord origin(ox, oy);

    const auto inside =
        [&](int x, int y) -> bool
        {
            return (x >= 0 && x < mask.width() &&
                    y >= 0 && y < mask.height());
        };

    const int dir[] = { 1, 1, 1, -1, -1, 1, -1, -1 };

    constexpr float bgcolor = 1.f;
    constexpr float fgcolor = 1.f - bgcolor;

    // first fill with background
    mask(guide.width(), guide.height());
    float *maskdata = mask;
    std::fill(maskdata, maskdata + (mask.width() * mask.height()), bgcolor);
    array2D<float> intersect;

    float min_feather = RT_INFINITY;

    for (const auto &area : areaMask.shapes) {
        Coord center(w2 + area.x / 100.0 * w2, h2 + area.y / 100.0 * h2);
        float area_w = area.width / 100.0 * width;
        float area_h = area.height / 100.0 * height;
    
        float a_min = area_w / 2;
        float b_min = area_h / 2;
        float r = b_min / a_min;
        float a_max = std::sqrt(2) * a_min;
        float a = a_max - area.roundness / 100.0 * (a_max - a_min);

        min_feather = std::min(a_min, b_min);

        const auto get =
            [&](int x, int y) -> Coord
            {
                PolarCoord p(Coord(x, y));
                double r, a;
                p.get(r, a);
                p.set(r, a - area.angle);
                Coord ret(p);
                ret += center;
                ret -= origin;
                return ret;
            };

        float **marr = mask;
        if (area.mode == AreaMask::Shape::INTERSECT) {
            intersect(mask.width(), mask.height());
            marr = intersect;
            float *p = intersect;
            std::fill(p, p + (mask.width() * mask.height()), bgcolor);
        }

        // draw the (bounded) ellipse
        for (int x = 0, n = int(a_min); x < n; ++x) {
            int yy = r * std::sqrt(a*a - float(x*x));
            for (int y = 0, m = std::min(yy, int(b_min)); y < m; ++y) {
                for (int d = 0; d < 4; ++d) {
                    int dx = dir[2*d], dy = dir[2*d+1];
                    Coord point = get(dx * x, dy * y);
                    for (int i = -1; i < 2; ++i) {
                        for (int j = -1; j < 2; ++j) {
                            if (inside(point.x+i, point.y+j)) {
                                switch (area.mode) {
                                case AreaMask::Shape::ADD:
                                case AreaMask::Shape::INTERSECT:
                                    marr[point.y+j][point.x+i] = fgcolor;
                                    break;
                                case AreaMask::Shape::SUBTRACT:
                                    marr[point.y+j][point.x+i] = bgcolor;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }

        if (area.mode == AreaMask::Shape::INTERSECT) {
#ifdef _OPENMP
#           pragma omp parallel for if (multithread)
#endif
            for (int y = 0; y < mask.height(); ++y) {
                for (int x = 0; x < mask.width(); ++x) {
                    if (mask[y][x] == fgcolor && intersect[y][x] != fgcolor) {
                        mask[y][x] = bgcolor;
                    }
                }
            }
        }
    }

    // guided feathering and contrast
    int radius = std::max(int(areaMask.feather / 100.0 * min_feather), 1);
    guidedFilter(guide, mask, mask, radius, 1e-7, multithread);

    DiagonalCurve curve(areaMask.contrast);
    const auto contrast =
        [&curve](float x) -> float
        {
            return curve.getVal(x);
        };
#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < mask.height(); ++y) {
        for (int x = 0; x < mask.width(); ++x) {
            float v = 1.f - LIM01(mask[y][x]);
            // if (!areaMask.inverted) {
            //     v = 1.f - v;
            // }
            mask[y][x] = contrast(v);
            assert(mask[y][x] == mask[y][x]);
        }
    }

    // and blur
    if (blur > 0.f) {
#ifdef _OPENMP
#       pragma omp parallel if (multithread)
#endif
        gaussianBlur(mask, mask, mask.width(), mask.height(), blur);
    }

    return true;
}


bool generate_drawn_mask(int ox, int oy, int width, int height, const DrawnMask &drawnMask, const array2D<float> &guide, bool multithread, array2D<float> &mask)
{
    if (drawnMask.isTrivial()) {
        return false;
    }

    const int mask_w = guide.width();
    const int mask_h = guide.height();

    const bool add = drawnMask.addmode;

    mask(mask_w, mask_h);
    float *maskdata = mask;
    const float bgcolor = 0.f;
    std::fill(maskdata, maskdata + (mask.width() * mask.height()), bgcolor);

    struct StrokeEval {
        StrokeEval(int ox, int oy, int width, int height,
                   int mask_w, int mask_h):
            ox_(ox), oy_(oy),
            width_(width), height_(height),
            mask_w_(mask_w), mask_h_(mask_h),
            radius_(-1), neg_(false)
        {
            curflag_ = 1;
            flag_.resize(mask_w_ * mask_h_, 0);
        }

        void update(const DrawnMask::Stroke &s)
        {
            int r = std::min(width_, height_) * s.radius * 0.25;

            if (r != radius_ || neg_ != s.erase || hardness_ != s.hardness) {
                radius_ = r;
                neg_ = s.erase;
                hardness_ = s.hardness;

                int w = 2*radius_ + 1;
                int h = 2*radius_ + 1;
                const float f = LIM01(hardness_);
                const float val = (neg_ ? -1.f : 1.f) + (1.f - f) * (neg_ ? 0.99f : -0.99f);
                
                buf_(w, h);
                for (int y = 0; y < h; ++y) {
                    for (int x = 0; x < w; ++x) {
                        float d = std::sqrt(SQR(x - radius_) + SQR(y - radius_));
                        if (d <= radius_) {
                            buf_[y][x] = val;
                        } else {
                            buf_[y][x] = RT_NAN;
                        }
                    }
                }

                ++curflag_;
            }

            int cx = width_ * s.x - ox_;
            int cy = height_ * s.y - oy_;
            lx = std::max(cx - radius_, 0);
            ly = std::max(cy - radius_, 0);
            ux = std::min(cx + radius_, mask_w_-1);
            uy = std::min(cy + radius_, mask_h_-1);
        }

        bool operator()(int x, int y, float &val)
        {
            val = buf_[y - ly][x - lx];
            size_t idx = y * mask_w_ + x;
            if (!xisnanf(val) && flag_[idx] != curflag_) {
                flag_[idx] = curflag_;
                // flagmods_.push_back(idx);
                return true;
            }
            return false;
        }

        bool is_bg(int x, int y) const
        {
            return flag_[y * mask_w_ + x] == 0;
        }

        int lx;
        int ly;
        int ux;
        int uy;

    private:
        int ox_;
        int oy_;
        int width_;
        int height_;
        int mask_w_;
        int mask_h_;
        int radius_;
        bool neg_;
        double hardness_;
        
        array2D<float> buf_;
        std::vector<uint16_t> flag_;
        std::vector<size_t> flagmods_;
        uint16_t curflag_;
    };

    double maxradius = 0.0;
    StrokeEval se(ox, oy, width, height, mask_w, mask_h);

    for (size_t i = 0; i < drawnMask.strokes.size(); ++i) {
        const auto &s = drawnMask.strokes[i];
        se.update(s);
        maxradius = std::max(maxradius, s.radius);
        for (int y = se.ly; y <= se.uy; ++y) {
            for (int x = se.lx; x <= se.ux; ++x) {
                float v;
                if (se(x, y, v)) {
                    //mask[y][x] = LIM(mask[y][x] + v, -1.f, 1.f);
                    if (add) {
                        if (signf(mask[y][x]) == signf(v)) {
                            mask[y][x] = LIM(mask[y][x] + v, -1.f, 1.f);
                        } else {
                            mask[y][x] = LIM(LIM01(mask[y][x]) + v, -1.f, 1.f);
                        }
                    } else {
                        mask[y][x] = LIM01(mask[y][x] + v);
                    }
                }
            }
        }
    }

    DiagonalCurve ccurve(drawnMask.contrast);
    const bool needscale = add && (drawnMask.smoothness > 0.f || drawnMask.feather > 0 || !ccurve.isIdentity());
    
    if (needscale) {
#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < mask_h; ++y) {
            for (int x = 0; x < mask_w; ++x) {
                if (se.is_bg(x, y)) {
                    mask[y][x] = 0.5f;
                } else {
                    mask[y][x] = (mask[y][x] + 1.f) / 2.f;
                }
            }
        }
    }

    if (drawnMask.smoothness > 0.f) {
        float sigma = std::min(width, height) * maxradius * 0.2f * drawnMask.smoothness;
#ifdef _OPENMP
#       pragma omp parallel if (multithread)
#endif
        gaussianBlur(mask, mask, mask_w, mask_h, sigma);
    }

    if (drawnMask.feather > 0) {
        int radius = int(drawnMask.feather / 100.0 * std::min(width, height) * 0.1 + 0.5);
        if (radius > 0) {
            guidedFilter(guide, mask, mask, radius, 1e-7, multithread);
        }
    }

    if (!ccurve.isIdentity()) {
#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < mask_h; ++y) {
            for (int x = 0; x < mask_w; ++x) {
                mask[y][x] = ccurve.getVal(mask[y][x]);
            }
        }
    }

    if (needscale) {
#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < mask_h; ++y) {
            for (int x = 0; x < mask_w; ++x) {
                mask[y][x] = (mask[y][x] * 2.f - 1.f);
            }
        }
    }

#if 0
    if (mask_w > 400) {
        Imagefloat tmp(mask_w, mask_h);
        for (int y = 0; y < mask_h; ++y) {
            for (int x = 0; x < mask_w; ++x) {
                tmp.r(y, x) = tmp.g(y, x) = tmp.b(y, x) = (mask[y][x] + 1.f) / 2.f * 65535.f;
            }
        }
        tmp.saveTIFF("/tmp/drawnmask.tif", 16);
    }
#endif
    
    return true;
}


template <class T>
void rgb2lab(Imagefloat::Mode mode, float R, float G, float B, float &L, float &a, float &b, const T ws[3][3])
{
    switch (mode) {
    case Imagefloat::Mode::RGB:
        Color::rgb2lab(R, G, B, L, a, b, ws);
        return;
    case Imagefloat::Mode::YUV:
        Color::yuv2rgb(G, B, R, R, G, B, ws);
        Color::rgb2lab(R, G, B, L, a, b, ws);
        return;
    case Imagefloat::Mode::XYZ:
        Color::XYZ2Lab(R, G, B, L, a, b);
        return;
    case Imagefloat::Mode::LAB:
        L = G;
        a = R;
        b = B;
        return;
    }
}


class DeltaEEvaluator {
public:
    DeltaEEvaluator(const std::vector<Mask> &masks)
    {
        masks_.reserve(masks.size());
        for (auto &m : masks) {
            masks_.push_back(m.deltaEMask);
            refs_.push_back(cmsCIELab());
            auto &de = masks_.back();
            auto &v = refs_.back();
            double h = de.H * 2.f * RT_PI / 360.f;
            v.L = de.L * 100.f;
            v.a = de.C * std::cos(h) * 100.f;
            v.b = de.C * std::sin(h) * 100.f;
        }
    }
    
    float operator()(size_t idx, float L, float a, float b)
    {
        auto &m = masks_[idx];
        if (!m.enabled || m.weight_L + m.weight_C + m.weight_H == 0) {
            return 1.f;
        }
        cmsCIELab c = { L * 100.f, a * 100.f, b * 100.f };
        auto &r = refs_[idx];
        const auto W = [](int w) -> double { return w ? 100.0 / w : 1000.0; };
        auto d = cmsCIE2000DeltaE(&r, &c, W(m.weight_L), W(m.weight_C), W(m.weight_H));
        return getval(m.range, 1.0 + LIM01(m.decay/100.0), d);
    }

private:
    float getval(float range, float decay, float d)
    {
        if (d <= range * 0.4f) {
            return 1.f;
        } else if (d >= range * 5.f) {
            return 0.f;
        } else {
            // Generalised logistic function (see Wikipedia)
            return 1.f - 1.f / (1 + 100.f * xexpf(-decay * (d - range)));
        }
    }
    
    std::vector<DeltaEMask> masks_;
    std::vector<cmsCIELab> refs_;
};

} // namespace


bool generateLabMasks(Imagefloat *rgb, const std::vector<Mask> &masks, int offset_x, int offset_y, int full_width, int full_height, double scale, bool multithread, int show_mask_idx, std::vector<array2D<float>> *Lmask, std::vector<array2D<float>> *abmask)
{
    int n = masks.size();
    if (show_mask_idx >= n || !masks[show_mask_idx].enabled) {
        show_mask_idx = -1;
    }
    std::vector<std::unique_ptr<FlatCurve>> hmask(n);
    std::vector<std::unique_ptr<FlatCurve>> cmask(n);
    std::vector<std::unique_ptr<FlatCurve>> lmask(n);

    const int W = rgb->getWidth();
    const int H = rgb->getHeight();
    const auto mode = rgb->mode();

    const Mask dflt;

    const int begin_idx = max(show_mask_idx, 0);
    const int end_idx = (show_mask_idx < 0 ? n : show_mask_idx+1);
    bool has_mask = false;

    for (int i = begin_idx; i < end_idx; ++i) {
        auto &r = masks[i];
        if (r.deltaEMask.enabled) {
            has_mask = true;
        }
        if (r.parametricMask.enabled && !r.parametricMask.hue.empty() && r.parametricMask.hue[0] != FCT_Linear && r.parametricMask.hue != dflt.parametricMask.hue) {
            hmask[i].reset(new FlatCurve(r.parametricMask.hue, true));
            has_mask = true;
        }
        if (r.parametricMask.enabled && !r.parametricMask.chromaticity.empty() && r.parametricMask.chromaticity[0] != FCT_Linear && r.parametricMask.chromaticity != dflt.parametricMask.chromaticity) {
            cmask[i].reset(new FlatCurve(r.parametricMask.chromaticity, false));
            has_mask = true;
        }
        if (r.parametricMask.enabled && !r.parametricMask.lightness.empty() && r.parametricMask.lightness[0] != FCT_Linear && r.parametricMask.lightness != dflt.parametricMask.lightness) {
            lmask[i].reset(new FlatCurve(r.parametricMask.lightness, false));
            has_mask = true;
        }
    }

    assert(!abmask || abmask->size() == size_t(n));
    assert(!Lmask || Lmask->size() == size_t(n));
    
    for (int i = begin_idx; i < end_idx; ++i) {
        if (abmask) {
            (*abmask)[i](W, H);
        }
        if (Lmask) {
            (*Lmask)[i](W, H);
        }
    }

    array2D<float> guide(W, H);
    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(rgb->colorSpace());
    float wp[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            wp[i][j] = ws[i][j];
        }
    }
    
    // magic constant c_factor: normally chromaticity is in [0; 42000] (see color.h), but here we use the constant to match how the chromaticity pipette works (see improcfun.cc lines 4705-4706 and color.cc line 1930
    constexpr float c_factor = 327.68f / 48000.f;

    DeltaEEvaluator dE(masks);

#ifdef _OPENMP
#       pragma omp parallel if (multithread)
#endif
    {
#ifdef __SSE2__
        float cBuffer[W];
        float hBuffer[W];
        float lBuffer[W];
        float aBuffer[W];
        float bBuffer[W];
#endif
#ifdef _OPENMP
#           pragma omp for schedule(dynamic, 16)
#endif
        for (int y = 0; y < H; ++y) {
#ifdef __SSE2__
            for (int x = 0; x < W; ++x) {
                rgb2lab(mode, rgb->r(y, x), rgb->g(y, x), rgb->b(y, x), lBuffer[x], aBuffer[x], bBuffer[x], wp);
            }
            if (has_mask) {
                // vectorized precalculation
                Color::Lab2Lch(aBuffer, bBuffer, cBuffer, hBuffer, W);
                fastlin2log(cBuffer, c_factor, 10.f, W);
            }
#endif
            for (int x = 0; x < W; ++x) {
#ifdef __SSE2__
                const float l = lBuffer[x] / 32768.f;
                const float a = aBuffer[x] / 42000.f;
                const float b = bBuffer[x] / 42000.f;
#else
                float l, a, b;
                rgb2lab(mode, rgb->r(y, x), rgb->g(y, x), rgb->b(y, x), l, a, b, wp);
                l /= 32768.f;
                a /= 42000.f;
                b /= 42000.f;
#endif
                guide[y][x] = LIM01(l);

                if (has_mask) {
#ifdef __SSE2__
                    // use precalculated values
                    const float c = cBuffer[x];
                    float h = hBuffer[x];
#else
                    float c, h;
                    Color::Lab2Lch(a, b, c, h);
                    c = xlin2log(c * c_factor, 10.f);
#endif
                        
                    h = Color::huelab_to_huehsv2(h);
                    h += 1.f/6.f; // offset the hue because we start from purple instead of red
                    if (h > 1.f) {
                        h -= 1.f;
                    }
                    h = xlin2log(h, 3.f);

                    for (int i = begin_idx; i < end_idx; ++i) {
                        auto &hm = hmask[i];
                        auto &cm = cmask[i];
                        auto &lm = lmask[i];
                        float blend = LIM01(dE(i, l, a, b) * (hm ? hm->getVal(h) : 1.f) * (cm ? cm->getVal(c) : 1.f) * (lm ? lm->getVal(l) : 1.f));
                        if (Lmask) {
                            (*Lmask)[i][y][x] = blend;
                        }
                        if (abmask) {
                            (*abmask)[i][y][x] = blend;
                        }
                    }
                }
            }
        }
    }

    if (has_mask) {
        for (int i = begin_idx; i < end_idx; ++i) {
            float blur = masks[i].parametricMask.blur;
            blur = blur < 0.f ? -1.f/blur : 1.f + blur;
            int r1 = max(int(4 / scale * blur + 0.5), 1);
            int r2 = max(int(25 / scale * blur + 0.5), 1);
            if (abmask) {
                rtengine::guidedFilter(guide, (*abmask)[i], (*abmask)[i], r1, 0.001, multithread);
            }
            if (Lmask) {
                rtengine::guidedFilter(guide, (*Lmask)[i], (*Lmask)[i], r2, 0.0001, multithread);
            }
        }
    } else {
        for (int i = begin_idx; i < end_idx; ++i) {
#ifdef _OPENMP
#           pragma omp parallel for if (multithread)
#endif
            for (int y = 0; y < H; ++y) {
                if (Lmask) {
                    float *v = (*Lmask)[i][y];
                    std::fill(v, v + W, 1.f);
                }
                if (abmask) {
                    float *v = (*abmask)[i][y];
                    std::fill(v, v + W, 1.f);
                }
            }
        }
    }

    if (full_width < 0) {
        full_width = W;
    }
    if (full_height < 0) {
        full_height = H;
    }

    array2D<float> amask;

    for (int i = begin_idx; i < end_idx; ++i) {
        if (generate_area_mask(offset_x, offset_y, full_width, full_height, guide, masks[i].areaMask, scale, multithread, amask)) {
#ifdef _OPENMP
#           pragma omp parallel for if (multithread)
#endif
            for (int y = 0; y < H; ++y) {
                for (int x = 0; x < W; ++x) {
                    if (abmask) {
                        (*abmask)[i][y][x] *= amask[y][x];
                    }
                    if (Lmask) {
                        (*Lmask)[i][y][x] *= amask[y][x];
                    }
                }
            }
        }
    }

    
    float s_scale = std::sqrt(scale);
    for (int i = begin_idx; i < end_idx; ++i) {
        if (masks[i].parametricMask.enabled && masks[i].parametricMask.contrastThreshold != 0) {
            amask(W, H);
            float thresh = float(std::abs(masks[i].parametricMask.contrastThreshold))/100.f * s_scale;
            float blur = std::max(masks[i].parametricMask.blur, 2.0) / s_scale;
            bool neg = masks[i].parametricMask.contrastThreshold < 0;
            buildBlendMask(guide, amask, W, H, thresh, 1.f, false, blur, 32768.f);
#ifdef _OPENMP
#           pragma omp parallel for if (multithread)
#endif
            for (int y = 0; y < H; ++y) {
                for (int x = 0; x < W; ++x) {
                    float f = neg ? 1.f - amask[y][x] : amask[y][x];
                    if (abmask) {
                        (*abmask)[i][y][x] *= f;
                    }
                    if (Lmask) {
                        (*Lmask)[i][y][x] *= f;
                    }
                }
            }
        }
    }

    for (int i = begin_idx; i < end_idx; ++i) {
        if (generate_drawn_mask(offset_x, offset_y, full_width, full_height, masks[i].drawnMask, guide, multithread, amask)) {
            const bool add = masks[i].drawnMask.addmode;
            const float alpha = 1.f - LIM01(masks[i].drawnMask.transparency);
#ifdef _OPENMP
#           pragma omp parallel for if (multithread)
#endif
            for (int y = 0; y < H; ++y) {
                for (int x = 0; x < W; ++x) {
                    const float f = alpha * amask[y][x];
                    if (add) {
                        if (abmask) {
                            (*abmask)[i][y][x] = LIM01((*abmask)[i][y][x] + f);
                        }
                        if (Lmask) {
                            (*Lmask)[i][y][x] = LIM01((*Lmask)[i][y][x] + f);
                        }
                    } else {
                        if (abmask) {
                            (*abmask)[i][y][x] *= f;
                        }
                        if (Lmask) {
                            (*Lmask)[i][y][x] *= f;
                        }
                    }
                }
            }
        }
    }
    
    for (int i = begin_idx; i < end_idx; ++i) {
        if (masks[i].inverted) {
#ifdef _OPENMP
#           pragma omp parallel for if (multithread)
#endif
            for (int y = 0; y < H; ++y) {
                for (int x = 0; x < W; ++x) {
                    if (abmask) {
                        (*abmask)[i][y][x] = 1.f - (*abmask)[i][y][x];
                    }
                    if (Lmask) {
                        (*Lmask)[i][y][x] = 1.f - (*Lmask)[i][y][x];
                    }
                }
            }
        }
    }

    if (show_mask_idx >= 0) {
        TMatrix iws = ICCStore::getInstance()->workingSpaceInverseMatrix(rgb->colorSpace());
        float iwp[3][3];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                iwp[i][j] = iws[i][j];
            }
        }
        
        auto *smask = abmask ? abmask : Lmask;
        
#ifdef _OPENMP
        #pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                auto blend = smask ? (*smask)[show_mask_idx][y][x] : 0.f;
                float l, a, b;
                rgb2lab(mode, rgb->r(y, x), rgb->g(y, x), rgb->b(y, x), l, a, b, wp);
                a = 0.f;
                b = blend * 42000.f;
                l = LIM(l + 32768.f * blend, 0.f, 32768.f);
                Color::lab2rgb(l, a, b, rgb->r(y, x), rgb->g(y, x), rgb->b(y, x), iwp);
            }
        }
        rgb->assignMode(Imagefloat::Mode::RGB);

        return false;
    }

    return true;
}


void fillPipetteLabMasks(Imagefloat *rgb, PlanarWhateverData<float>* editWhatever, LabMasksEditID id, bool multithread)
{
    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(rgb->colorSpace());
    float wp[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            wp[i][j] = ws[i][j];
        }
    }
    const int W = rgb->getWidth();
    const int H = rgb->getHeight();
    const auto mode = rgb->mode();
    
#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            float v = 0.f;
            float l, a, b;
            rgb2lab(mode, rgb->r(y, x), rgb->g(y, x), rgb->b(y, x), l, a, b, wp);
            switch (id) {
            case LabMasksEditID::H:
                v = Color::huelab_to_huehsv2(xatan2f(b, a));
                break;
            case LabMasksEditID::C:
                v = LIM01<float>(std::sqrt(SQR(a) + SQR(b) + 0.001f) / 48000.f);
                break;
            case LabMasksEditID::L:
                v = LIM01<float>(l / 32768.f);
                break;
            }
            editWhatever->v(y, x) = v;
        }
    }
}


bool getDeltaEColor(Imagefloat *rgb, int x, int y, int offset_x, int offset_y, int full_width, int full_height, double scale, float &L, float &C, float &H)
{
    std::vector<float> med_L, med_a, med_b;
    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(rgb->colorSpace());
    const auto mode = rgb->mode();
    x = x / scale - offset_x;
    y = y / scale - offset_y;
    const int off = int(32.0 / scale + 0.5);
    if (x < 0 || x >= rgb->getWidth() || y < 0 || y >= rgb->getHeight()) {
        return false;
    }
    for (int i = max(y-off, 0), end = min(y+off, rgb->getHeight()); i < end; ++i) {
        for (int j = max(x-off, 0), end = min(x+off, rgb->getWidth()); j < end; ++j) {
            float L, a, b;
            rgb2lab(mode, rgb->r(i, j), rgb->g(i, j), rgb->b(i, j), L, a, b, ws);
            med_L.push_back(L);
            med_a.push_back(a);
            med_b.push_back(b);
        }
    }

    std::sort(med_L.begin(), med_L.end());
    std::sort(med_a.begin(), med_a.end());
    std::sort(med_b.begin(), med_b.end());

    auto idx = med_L.size()/2;
    L = med_L[idx] / 32768.f;
    float a = med_a[idx] / 42000.f;
    float b = med_b[idx] / 42000.f;
    C = sqrtf(SQR(a) + SQR(b));
    H = xatan2f(b, a) * 360.f / (2.f * RT_PI);
    if (H < 0.f) {
        H += 360.f;
    } else if (H > 360.f) {
        H -= 360.f;
    }
    return true;
}

} // namespace rtengine
