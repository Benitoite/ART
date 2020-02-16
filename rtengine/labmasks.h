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

#pragma once

#include "procparams.h"
#include "array2D.h"
#include "labimage.h"
#include "imagefloat.h"


namespace rtengine {

bool generateLabMasks(Imagefloat *rgb, const std::vector<procparams::Mask> &masks, int offset_x, int offset_y, int full_width, int full_height, double scale, bool multithread, int show_mask_idx, std::vector<array2D<float>> *Lmask, std::vector<array2D<float>> *abmask);

enum class LabMasksEditID { H = 0, C, L };
void fillPipetteLabMasks(Imagefloat *rgb, PlanarWhateverData<float> *editWhatever, LabMasksEditID id, bool multithread);

bool getDeltaEColor(Imagefloat *rgb, int x, int y, int offset_x, int offset_y, int full_width, int full_height, double scale, float &L, float &C, float &H);

} // namespace rtengine
