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

namespace rtengine {

void ImProcFunctions::blackAndWhite(Imagefloat *img)
{
    if (!params->blackwhite.enabled) {
        return;
    }

    img->setMode(Imagefloat::Mode::RGB, multiThread);

    float bwr = float(params->blackwhite.mixerRed);
    float bwg = float(params->blackwhite.mixerGreen);
    float bwb = float(params->blackwhite.mixerBlue);
    float bwrgam = float(params->blackwhite.gammaRed);
    float bwggam = float(params->blackwhite.gammaGreen);
    float bwbgam = float(params->blackwhite.gammaBlue);
    bool chmix = params->blackwhite.method == "ChannelMixer";

    float gamvalr = 125.f;
    float gamvalg = 125.f;
    float gamvalb = 125.f;

    if (bwrgam < 0) {
        gamvalr = 100.f;
    }

    if (bwggam < 0) {
        gamvalg = 100.f;
    }

    if (bwbgam < 0) {
        gamvalb = 100.f;
    }

    float gammabwr = 1.f;
    float gammabwg = 1.f;
    float gammabwb = 1.f;
    {
        gammabwr = 1.f - bwrgam / gamvalr;
        gammabwg = 1.f - bwggam / gamvalg;
        gammabwb = 1.f - bwbgam / gamvalb;
    }
    bool hasgammabw = gammabwr != 1.f || gammabwg != 1.f || gammabwb != 1.f;
    
    const int W = img->getWidth();
    const int H = img->getHeight();

    if (!chmix) {
#ifdef _OPENMP
#       pragma omp parallel for if (multiThread)
#endif
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                float r = img->r(y, x);
                float g = img->g(y, x);
                float b = img->b(y, x);
                float Y = Color::rgbLuminance(r, g, b);
                r = g = b = Y;
                if (hasgammabw) {
                    Color::trcGammaBW(r, g, b, gammabwr, gammabwg, gammabwb);
                }
                img->r(y, x) = r;
                img->g(y, x) = g;
                img->b(y, x) = b;
            }
        }
    } else {
        float kcorec = 1.f;
        float filcor;
        double rrm, ggm, bbm;
        Color::computeBWMixerConstants(params->blackwhite.setting, params->blackwhite.filter, "", filcor, bwr, bwg, bwb, kcorec, rrm, ggm, bbm);

#ifdef _OPENMP
#       pragma omp parallel for if (multiThread)
#endif
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                img->r(y, x) = img->g(y, x) = img->b(y, x) = ((bwr * img->r(y, x) + bwg * img->g(y, x) + bwb * img->b(y, x)) * kcorec);

#ifndef __SSE2__
                if (hasgammabw) {
                    Color::trcGammaBW(img->r(y, x), img->g(y, x), img->b(y, x), gammabwr, gammabwg, gammabwb);
                }
#endif
            }

#ifdef __SSE2__
            if (hasgammabw) {
                Color::trcGammaBWRow(img->r(y), img->g(y), img->b(y), W, gammabwr, gammabwg, gammabwb);
            }
#endif
        }
    }
}

} // namespace rtengine
