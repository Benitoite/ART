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

namespace rtengine {

void ImProcFunctions::channelMixer(Imagefloat *img)
{
    const bool mixchannels = params->chmixer.enabled &&
        (params->chmixer.red[0] != 100 ||
         params->chmixer.red[1] != 0 ||
         params->chmixer.red[2] != 0 ||
         params->chmixer.green[0] != 0 ||
         params->chmixer.green[1] != 100 ||
         params->chmixer.green[2] != 0 ||
         params->chmixer.blue[0] != 0 ||
         params->chmixer.blue[1] != 0 ||
         params->chmixer.blue[2] != 100);

    if (mixchannels) {
        img->setMode(Imagefloat::Mode::RGB, multiThread);
        
        const float chMixRR = float(params->chmixer.red[0])/10.f;
        const float chMixRG = float(params->chmixer.red[1])/10.f;
        const float chMixRB = float(params->chmixer.red[2])/10.f;
        const float chMixGR = float(params->chmixer.green[0])/10.f;
        const float chMixGG = float(params->chmixer.green[1])/10.f;
        const float chMixGB = float(params->chmixer.green[2])/10.f;
        const float chMixBR = float(params->chmixer.blue[0])/10.f;
        const float chMixBG = float(params->chmixer.blue[1])/10.f;
        const float chMixBB = float(params->chmixer.blue[2])/10.f;
        
#ifdef _OPENMP
#       pragma omp parallel for if (multiThread)
#endif
        for (int y = 0; y < img->getHeight(); ++y) {
            for (int x = 0; x < img->getWidth(); ++x) {
                float r = img->r(y, x);
                float g = img->g(y, x);
                float b = img->b(y, x);

                float rmix = (r * chMixRR + g * chMixRG + b * chMixRB) / 100.f;
                float gmix = (r * chMixGR + g * chMixGG + b * chMixGB) / 100.f;
                float bmix = (r * chMixBR + g * chMixBG + b * chMixBB) / 100.f;

                img->r(y, x) = rmix;
                img->g(y, x) = gmix;
                img->b(y, x) = bmix;
            }
        }
    }
}

} // namespace rtengine
