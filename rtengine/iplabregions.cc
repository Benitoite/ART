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

namespace rtengine {

void ImProcFunctions::labColorCorrectionRegions(LabImage *lab)
{
    if (!params->colorToning.enabled || params->colorToning.method != "LabRegions") {
        return;
    }

    const float factor = ColorToningParams::LABGRID_CORR_MAX * 3.f;
    const float scaling = 1.f;

    size_t n = params->colorToning.labregions.size();
    std::vector<std::unique_ptr<FlatCurve>> hmask(n);
    std::vector<std::unique_ptr<FlatCurve>> cmask(n);
    std::vector<std::unique_ptr<FlatCurve>> lmask(n);

    for (size_t i = 0; i < n; ++i) {
        auto &r = params->colorToning.labregions[i];
        hmask[i].reset(new FlatCurve(r.hueMask, true));
        cmask[i].reset(new FlatCurve(r.chromaticityMask, false));
        lmask[i].reset(new FlatCurve(r.lightnessMask, false));
    }

    float min_c = RT_INFINITY_F, max_c = -RT_INFINITY_F;
    float min_l = RT_INFINITY_F, max_l = -RT_INFINITY_F;

    const auto logscale =
        [](float x, float base) -> float
        {
            return std::log(x * (base - 1.0f) + 1.0f) / std::log(base);
        };
    
#ifdef _OPENMP
    #pragma omp parallel for if (multiThread)
#endif
    for (int y = 0; y < lab->H; ++y) {
        for (int x = 0; x < lab->W; ++x) {
            float l = lab->L[y][x];
            float a = lab->a[y][x];
            float b = lab->b[y][x];
            float c, h;
            Color::Lab2Lch(a, b, c, h);
            float c1 = logscale(c / 42000.f, 1e6);
            float h1 = (RT_PI + h) * RT_1_PI;
            float l1 = l / 32768.f;

#ifdef _OPENMP
            #pragma omp critical
            {
                min_c = min(min_c, c1);
                max_c = max(max_c, c1);
                min_l = min(min_l, l1);
                max_l = max(max_l, l1);
            }
#endif
            
            for (size_t i = 0; i < n; ++i) {
                auto &r = params->colorToning.labregions[i];
                auto &hm = hmask[i];
                auto &cm = cmask[i];
                auto &lm = lmask[i];
                float blend = hm->getVal(h1) * cm->getVal(c1) * lm->getVal(l1);
                float s = 1.f + r.saturation / 100.f;
                float a_new = s * (a + 32768.f * r.a / factor / scaling);
                float b_new = s * (b + 32768.f * r.b / factor / scaling);
                float l_new = l * (1.f + r.lightness / 100.f);
                lab->L[y][x] = intp(blend, l_new, l);
                lab->a[y][x] = intp(blend, a_new, a);
                lab->b[y][x] = intp(blend, b_new, b);
            }
        }
    }

    std::cout << "LAB regions: min_c = " << min_c << ", max_c = " << max_c
              << ", min_l = " << min_l << ", max_l = " << max_l << std::endl;
}

} // namespace rtengine
