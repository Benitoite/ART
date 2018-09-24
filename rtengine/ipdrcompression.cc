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

namespace {

void gamma_compress(Imagefloat *rgb, float power, float slope, bool multithread)
{
    int w = rgb->getWidth();
    int h = rgb->getHeight();

    const auto apply =
        [=](float x) -> float
        {
            x = (x / 65535.0) * slope;
            if (x >= 0) {
                return powf(x, power) * 65535.0;
            } else {
                return x;
            }
        };

#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < h; y++) {
        for (int x = 0; x < w; x++) {
            rgb->r(y, x) = apply(rgb->r(y, x));
            rgb->g(y, x) = apply(rgb->g(y, x));
            rgb->b(y, x) = apply(rgb->b(y, x));
        }
    }
}

} // namespace


void ImProcFunctions::dynamicRangeCompression(Imagefloat *rgb)
{
    if (params->drcomp.enabled) {
        if (params->drcomp.method == procparams::DRCompressionParams::DR_COMP_GAMMA) {
            gamma_compress(rgb, params->drcomp.power, params->drcomp.slope, multiThread);
        } else {
            ToneMapFattal02(rgb);
        }
    }
}

} // namespace rtengine
