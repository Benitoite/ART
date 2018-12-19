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
#include "sleef.c"
#include "coord.h"
#include "gauss.h"
#include "color.h"
#include "curves.h"


namespace rtengine {

using procparams::AreaMask;
using procparams::LabCorrectionMask;

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


bool generate_area_mask(int ox, int oy, int width, int height, array2D<float> &mask, const AreaMask &areaMask, bool enabled, float blur, bool multithread)
{
    if (!enabled || areaMask.shapes.empty() || areaMask.isTrivial()) {
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
    array2D<float> guide(mask.width(), mask.height(), mask, 0);
    float *maskdata = mask;
    std::fill(maskdata, maskdata + (mask.width() * mask.height()), bgcolor);

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
                                mask[point.y+j][point.x+i] = fgcolor;
                            }
                        }
                    }
                }
            }
        }
    }

    // guided feathering and contrast
    int radius = std::max(int(areaMask.feather / 100.0 * min_feather), 1);
    guidedFilter(guide, mask, mask, radius, 1e-7, multithread);
    const float c = float(areaMask.contrast) / 4.f;
    const auto contrast =
        [c](float x) -> float
        {
            if (c <= 0) {
                return x;
            }
            constexpr float s = 1.f;
            constexpr float a = 0.5f;
            float y = 0.f;
            if (x <= 0.5f) {
                y = a * std::pow(x/a, c);
            } else {
                y = 1.f - (1-a) * std::pow((1 - x) / (1 - a), c);
            }
            return s*y + (1.f-s)*x;
        };
#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < mask.height(); ++y) {
        for (int x = 0; x < mask.width(); ++x) {
            float v = LIM01(mask[y][x]);
            if (!areaMask.inverted) {
                v = 1.f - v;
            }
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

} // namespace


bool generateLabMasks(LabImage *lab, const std::vector<LabCorrectionMask> &masks, int offset_x, int offset_y, int full_width, int full_height, double scale, bool multithread, int show_mask_idx, std::vector<array2D<float>> *Lmask, std::vector<array2D<float>> *abmask)
{
    int n = masks.size();
    if (show_mask_idx >= n) {
        show_mask_idx = -1;
    }
    std::vector<std::unique_ptr<FlatCurve>> hmask(n);
    std::vector<std::unique_ptr<FlatCurve>> cmask(n);
    std::vector<std::unique_ptr<FlatCurve>> lmask(n);

    const int begin_idx = max(show_mask_idx, 0);
    const int end_idx = (show_mask_idx < 0 ? n : show_mask_idx+1);

    for (int i = begin_idx; i < end_idx; ++i) {
        auto &r = masks[i];
        if (!r.hueMask.empty() && r.hueMask[0] != FCT_Linear) {
            hmask[i].reset(new FlatCurve(r.hueMask, true));
        }
        if (!r.chromaticityMask.empty() && r.chromaticityMask[0] != FCT_Linear) {
            cmask[i].reset(new FlatCurve(r.chromaticityMask, false));
        }
        if (!r.lightnessMask.empty() && r.lightnessMask[0] != FCT_Linear) {
            lmask[i].reset(new FlatCurve(r.lightnessMask, false));
        }
    }

    assert(!abmask || abmask->size() == size_t(n));
    assert(!Lmask || Lmask->size() == size_t(n));
    
    for (int i = begin_idx; i < end_idx; ++i) {
        if (abmask) {
            (*abmask)[i](lab->W, lab->H);
        }
        if (Lmask) {
            (*Lmask)[i](lab->W, lab->H);
        }
    }

    array2D<float> guide(lab->W, lab->H);

    // magic constant c_factor: normally chromaticity is in [0; 42000] (see color.h), but here we use the constant to match how the chromaticity pipette works (see improcfun.cc lines 4705-4706 and color.cc line 1930
    constexpr float c_factor = 327.68f / 48000.f;

#ifdef _OPENMP
    #pragma omp parallel if (multithread)
#endif
    {
#ifdef __SSE2__
        float cBuffer[lab->W];
        float hBuffer[lab->W];
#endif
#ifdef _OPENMP
        #pragma omp for schedule(dynamic, 16)
#endif
        for (int y = 0; y < lab->H; ++y) {
#ifdef __SSE2__
            // vectorized precalculation
            Color::Lab2Lch(lab->a[y], lab->b[y], cBuffer, hBuffer, lab->W);
            fastlin2log(cBuffer, c_factor, 10.f, lab->W);
#endif
            for (int x = 0; x < lab->W; ++x) {
                const float l = lab->L[y][x] / 32768.f;
                guide[y][x] = LIM01(l);
#ifdef __SSE2__
                // use precalculated values
                const float c = cBuffer[x];
                float h = hBuffer[x];
#else
                float c, h;
                Color::Lab2Lch(lab->a[y][x], lab->b[y][x], c, h);
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
                    float blend = LIM01((hm ? hm->getVal(h) : 1.f) * (cm ? cm->getVal(c) : 1.f) * (lm ? lm->getVal(l) : 1.f));
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

    for (int i = begin_idx; i < end_idx; ++i) {
        float blur = masks[i].maskBlur;
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

    if (full_width < 0) {
        full_width = lab->W;
    }
    if (full_height < 0) {
        full_height = lab->H;
    }

    for (int i = begin_idx; i < end_idx; ++i) {
        if (generate_area_mask(offset_x, offset_y, full_width, full_height, guide, masks[i].areaMask, masks[i].areaEnabled, masks[i].maskBlur, multithread)) {
#ifdef _OPENMP
#           pragma omp parallel for if (multithread)
#endif
            for (int y = 0; y < lab->H; ++y) {
                for (int x = 0; x < lab->W; ++x) {
                    if (abmask) {
                        (*abmask)[i][y][x] *= guide[y][x];
                    }
                    if (Lmask) {
                        (*Lmask)[i][y][x] *= guide[y][x];
                    }
                }
            }
        }
    }

    if (show_mask_idx >= 0) {
        auto *smask = abmask ? abmask : Lmask;
        
#ifdef _OPENMP
        #pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < lab->H; ++y) {
            for (int x = 0; x < lab->W; ++x) {
                auto blend = smask ? (*smask)[show_mask_idx][y][x] : 0.f;
                lab->a[y][x] = 0.f;
                lab->b[y][x] = blend * 42000.f;
                lab->L[y][x] = LIM(lab->L[y][x] + 32768.f * blend, 0.f, 32768.f);
            }
        }

        return false;
    }

    return true;
}


void fillPipetteLabMasks(LabImage *lab, PlanarWhateverData<float>* editWhatever, LabMasksEditID id, bool multithread)
{
#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < lab->H; ++y) {
        for (int x = 0; x < lab->W; ++x) {
            float v = 0.f;
            switch (id) {
            case LabMasksEditID::H:
                v = Color::huelab_to_huehsv2(xatan2f(lab->b[y][x], lab->a[y][x]));
                break;
            case LabMasksEditID::C:
                v = LIM01<float>(std::sqrt(SQR(lab->a[y][x]) + SQR(lab->b[y][x]) + 0.001f) / 48000.f);
                break;
            case LabMasksEditID::L:
                v = LIM01<float>(lab->L[y][x] / 32768.f);
                break;
            }
            editWhatever->v(y, x) = v;
        }
    }
}

} // namespace rtengine
