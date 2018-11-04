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

namespace rtengine {

namespace {

// taken from darktable (src/iop/profile_gamma.c)
/*
   copyright (c) 2009--2010 johannes hanika.
   copyright (c) 2014 LebedevRI.
   copyright (c) 2018 Aurélien Pierre.

   darktable is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   darktable is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/ 
void log_compress(Imagefloat *rgb, float dynamic_range, float gray_point, float shadows_range, bool multithread)
{
    int w = rgb->getWidth();
    int h = rgb->getHeight();

    const float gray = (gray_point / 100.f) * 65535.f;
    const float noise = pow_F(2.f, -16.f);
    const float log2 = xlogf(2.f);

    const auto apply =
        [=](float x) -> float
        {
            x = max(x / gray, noise);
            x = max((xlogf(x)/log2 - shadows_range) / dynamic_range, noise);
            return x * 65535.f;
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
        if (params->drcomp.method == procparams::DRCompressionParams::DR_COMP_LOG) {
            log_compress(rgb, params->drcomp.dynamicRange, params->drcomp.grayPoint, params->drcomp.shadowsRange, multiThread);
        } else {
            ToneMapFattal02(rgb);
        }
    }
}

} // namespace rtengine
