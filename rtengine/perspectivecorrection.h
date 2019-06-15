/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 Alberto Griggio <alberto.griggio@gmail.com>
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

#include "coord2d.h"
#include "procparams.h"

namespace rtengine {

class PerspectiveCorrection {
public:
    PerspectiveCorrection();
    void init(int width, int height, const procparams::PerspectiveParams &params, bool fill);
    void operator()(double &x, double &y);

private:
    void correct(double &x, double &y, double scale, double offx, double offy);
    void calc_scale(int w, int h, const procparams::PerspectiveParams &params, bool fill);
    bool test_scale(int w, int h, double scale);
    
    bool ok_;
    double scale_;
    double offx_;
    double offy_;
    double scalein_;
    double offxin_;
    double offyin_;
    float ihomograph_[3][3];
};

} // namespace rtengine
