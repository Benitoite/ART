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

#include "improcfun.h"
#include "gauss.h"
#include "sleef.c"
#include "opthelper.h"
#include "guidedfilter.h"

namespace rtengine {

namespace {

void sh(Imagefloat *img, const ProcParams *params, double scale, bool multithread)
{
    const int width = img->getWidth();
    const int height = img->getHeight();

    array2D<float> mask(width, height);
    array2D<float> L(width, height);
    const float radius = float(params->sh.radius) * 10 / scale;
    LUTf f(32768);

    const auto apply =
        [&](int amount, int tonalwidth, bool hl) -> void
        {
            const float thresh = tonalwidth * 327.68f;
            const float scale = hl ? (thresh > 0.f ? 0.9f / thresh : 1.f) : thresh * 0.9f;

#ifdef _OPENMP
#               pragma omp parallel for if (multithread)
#endif
            for (int y = 0; y < height; ++y) {
                for (int x = 0; x < width; ++x) {
                    float l = img->g(y, x);
                    float l1 = l / 32768.f;
                    if (hl) {
                        mask[y][x] = (l > thresh) ? 1.f : pow4(l * scale);
                        L[y][x] = 1.f - l1;
                    } else {
                        mask[y][x] = l <= thresh ? 1.f : pow4(scale / l);
                        L[y][x] = l1;
                    }
                }
            }

            rtengine::guidedFilter(L, mask, mask, radius, 0.075, multithread, 4);

            const float base = std::pow(4.f, float(std::abs(amount))/100.f);
            const float gamma = 3.f * (hl == (amount >= 0)) ? base : 1.f / base;

            const float contrast = std::pow(3.f*(1.f-std::pow(float(amount)/100.f, 1.3f))+2.f, float(min(std::abs(amount), 50))/100.f);
            DiagonalCurve sh_contrast({
                    DCT_NURBS,
                    0, 0,
                    0.125, std::pow(0.125 / 0.25, contrast) * 0.25, 
                    1, 1
                });

#ifdef _OPENMP
#               pragma omp parallel for if (multithread)
#endif
            for (int l = 0; l < 32768; ++l) {
                auto base = pow_F(l / 32768.f, gamma);
                if (!hl) {
                    base = sh_contrast.getVal(base);
                }
                f[l] = base * 32768.f;
            }

#ifdef _OPENMP
#               pragma omp parallel for schedule(dynamic,16) if (multithread)
#endif
            for (int y = 0; y < height; ++y) {
                for (int x = 0; x < width; ++x) {
                    float blend = LIM01(mask[y][x]);
                    float orig = 1.f - blend;
                    float l = img->g(y, x);
                    if (l >= 0 && l < 32768.f) {
                        float ll = intp(blend, f[l], l);
                        if (!hl && l > 1.f) {
                            float &a = img->r(y, x);
                            float &b = img->b(y, x);
                            // when pushing shadows, scale also the chromaticity
                            float s = max(ll / l * 0.5f, 1.f) * blend;
                            a = a * s + a * orig;
                            b = b * s + b * orig;
                        }
                        img->g(y, x) = ll;
                    }
                }
            }
        };

    if (params->sh.highlights) {
        apply(params->sh.highlights * 0.7, params->sh.htonalwidth, true);
    }

    if (params->sh.shadows) {
        apply(params->sh.shadows * 0.6, params->sh.stonalwidth, false);
    }
}

} // namespace


void ImProcFunctions::shadowsHighlights(Imagefloat *rgb)
{
    if (!params->sh.enabled || (!params->sh.highlights && !params->sh.shadows)){
        return;
    }

    rgb->setMode(Imagefloat::Mode::LAB, multiThread);
    sh(rgb, params, scale, multiThread);
}

} // namespace rtengine
