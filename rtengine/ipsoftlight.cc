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

inline float sl(float blend, float x)
{
    if (!OOG(x)) {
        const float orig = 1.f - blend;
        float v = Color::gamma_srgb(x) / MAXVALF;
        // Pegtop's formula from
        // https://en.wikipedia.org/wiki/Blend_modes#Soft_Light
        float v2 = v * v;
        float v22 = v2 * 2.f;
        v = v2 + v22 - v22 * v;
        x = blend * Color::igamma_srgb(v * MAXVALF) + orig * x;
    }
    return x;
}

} // namespace


void ImProcFunctions::softLight(Imagefloat *rgb)
{
    const bool sl_enabled = params->softlight.enabled && params->softlight.strength > 0;
    if (!sl_enabled) {
        return;
    }

    const float blend = params->softlight.strength / 100.f;

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int y = 0; y < rgb->getHeight(); ++y) {
        for (int x = 0; x < rgb->getWidth(); ++x) {
            rgb->r(y, x) = sl(blend, rgb->r(y, x));
            rgb->g(y, x) = sl(blend, rgb->g(y, x));
            rgb->b(y, x) = sl(blend, rgb->b(y, x));
        }
    }
}

} // namespace rtengine
