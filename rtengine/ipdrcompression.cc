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


void ImProcFunctions::dynamicRangeCompression(Imagefloat *rgb)
{
    ToneMapFattal02(rgb);
}


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
void ImProcFunctions::logEncoding(float *r, float *g, float *b, int istart, int jstart, int tW, int tH, int TS)
{
    if (!params->logenc.enabled) {
        return;
    }
    
    const float gray = (params->logenc.grayPoint / 100.f) * 65535.f;
    const float shadows_range = params->logenc.shadowsRange;
    const float dynamic_range = params->logenc.dynamicRange;
    const float noise = pow_F(2.f, -16.f);
    const float log2 = xlogf(2.f);
    const float base = pow_F(2.f, -params->logenc.brightness);

    const auto apply =
        [=](float x) -> float
        {
            x = max(x / gray, noise);
            x = max((xlogf(x)/log2 - shadows_range) / dynamic_range, noise);
            return (base != 1.f ? xlog2lin(x, base) : x) * 65535.f;
        };

    for (int i = istart, ti = 0; i < tH; i++, ti++) {
        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
            int idx = ti * TS + tj;
            r[idx] = apply(r[idx]);
            g[idx] = apply(g[idx]);
            b[idx] = apply(b[idx]);
        }
    }
}

} // namespace rtengine
