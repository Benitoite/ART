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
// extracted and datapted from ImProcFunctions (improcfun.cc,
// dirpyr_equalizer.cc) of RawTherapee

#pragma once

#include "improcfun.h"

namespace rtengine { namespace cbdl {

// pyramid wavelet
void dirpyr_equalizer(float ** src, float ** dst, int srcwidth, int srcheight, float ** l_a, float ** l_b, const double * mult, const double dirpyrThreshold, const double skinprot, float b_l, float t_l, float t_r, float scale, bool multiThread);    //Emil's directional pyramid wavelet


}} // namespace rtengine::cbdl
