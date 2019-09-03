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

void ImProcFunctions::rgbCurves(Imagefloat *img)
{
    PlanarWhateverData<float> *editWhatever = nullptr;
    EditUniqueID eid = pipetteBuffer ? pipetteBuffer->getEditID() : EUID_None;
    if ((eid == EUID_RGB_R || eid == EUID_RGB_G || eid == EUID_RGB_B) && pipetteBuffer->getDataProvider()->getCurrSubscriber()->getPipetteBufferType() == BT_SINGLEPLANE_FLOAT) {
        editWhatever = pipetteBuffer->getSinglePlaneBuffer();
    }
    
    if (!params->rgbCurves.enabled) {
        if (editWhatever) {
            editWhatever->fill(0.f);
        }
        return;
    }
    
    img->setMode(Imagefloat::Mode::RGB, multiThread);

    LUTf rCurve, gCurve, bCurve;
    CurveFactory::RGBCurve(params->rgbCurves.rcurve, rCurve, scale);
    CurveFactory::RGBCurve(params->rgbCurves.gcurve, gCurve, scale);
    CurveFactory::RGBCurve(params->rgbCurves.bcurve, bCurve, scale);

    const int W = img->getWidth();
    const int H = img->getHeight();

    if (editWhatever) {
        float **chan = nullptr;
        switch (eid) {
        case EUID_RGB_R:
            chan = img->r.ptrs;
            break;
        case EUID_RGB_G:
            chan = img->g.ptrs;
            break;
        case EUID_RGB_B:
            chan = img->b.ptrs;
            break;
        default:
            assert(false);
        }

#ifdef _OPENMP
#       pragma omp parallel for if (multiThread)
#endif
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                editWhatever->v(y, x) = Color::gamma2curve[chan[y][x]] / 65536.f;
            }
        }
    }

    if (rCurve || gCurve || bCurve) { // if any of the RGB curves is engaged
#ifdef _OPENMP
#       pragma omp parallel for if (multiThread)
#endif
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                if (rCurve) {
                    img->r(y, x) = rCurve[img->r(y, x)];
                }
                if (gCurve) {
                    img->g(y, x) = gCurve[img->g(y, x)];
                }
                if (bCurve) {
                    img->b(y, x) = bCurve[img->b(y, x)];
                }
            }
        }
    }
}

} // namespace rtengine
