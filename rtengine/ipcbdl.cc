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
#include "labmasks.h"
#include "array2D.h"

namespace rtengine {

bool ImProcFunctions::contrastByDetailLevels(Imagefloat *rgb, int offset_x, int offset_y, int full_width, int full_height)
{
    PlanarWhateverData<float> *editWhatever = nullptr;
    EditUniqueID eid = pipetteBuffer ? pipetteBuffer->getEditID() : EUID_None;

    if ((eid == EUID_LabMasks_H2 || eid == EUID_LabMasks_C2 || eid == EUID_LabMasks_L2) && pipetteBuffer->getDataProvider()->getCurrSubscriber()->getPipetteBufferType() == BT_SINGLEPLANE_FLOAT) {
        editWhatever = pipetteBuffer->getSinglePlaneBuffer();
    }
    
    if (params->dirpyrequalizer.enabled && rgb->getWidth() >= 8 && rgb->getHeight() >= 8) {
        LabImage tmplab(rgb->getWidth(), rgb->getHeight());
        rgb2lab(*rgb, tmplab);
        LabImage *lab = &tmplab;
        
        if (editWhatever) {
            LabMasksEditID id = static_cast<LabMasksEditID>(int(eid) - EUID_LabMasks_H2);
            fillPipetteLabMasks(lab, editWhatever, id, multiThread);
        }
        
        int n = params->dirpyrequalizer.levels.size();
        int show_mask_idx = params->dirpyrequalizer.showMask;
        if (show_mask_idx >= n) {
            show_mask_idx = -1;
        }
        std::vector<array2D<float>> mask(n);
        if (!generateLabMasks(lab, params->dirpyrequalizer.labmasks, offset_x, offset_y, full_width, full_height, scale, multiThread, show_mask_idx, &mask, nullptr)) {
            lab2rgb(*lab, *rgb);
            return true; // show mask is active, nothing more to do
        }

        array2D<float> L(lab->W, lab->H, lab->L, 0);

        // double mult[6];
        // const double scale_factor = 1.0;//min(1.5 / scale, 1.0);

        for (int i = 0; i < n; ++i) {
            auto &l = params->dirpyrequalizer.levels[i];
            // for (int k = 0; k < 6; ++k) {
            //     mult[k] = 1.0 + (l.mult[k] - 1.0) * scale_factor;
            // }
            const double threshold = l.threshold / scale;
            dirpyr_equalizer(lab->L, L, lab->W, lab->H, lab->a, lab->b, l.mult, /*l.*/threshold, 0.0, 0.f, 0.f, 0.f, std::max(scale, 1.0));
            const auto &blend = mask[i];

#ifdef _OPENMP
#           pragma omp parallel for if (multiThread)
#endif
            for (int y = 0; y < lab->H; ++y) {
                for (int x = 0; x < lab->W; ++x) {
                    // float l = lab->L[y][x];
                    lab->L[y][x] = intp(blend[y][x], L[y][x], lab->L[y][x]);
                    // L[y][x] = l;
                }
            }
        }

        lab2rgb(*lab, *rgb);
    } else if (editWhatever) {
        editWhatever->fill(0.f);
    }

    return false;
}

} // namespace rtengine
