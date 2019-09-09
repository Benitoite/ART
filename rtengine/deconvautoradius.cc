/* -*- C++ -*-
 *
 *  This file is part of ART.
 *
 *  Copyright (c) 2019 Alberto Griggio <alberto.griggio@gmail.com>
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

#include "rawimagesource.h"
//#define BENCHMARK
#include "StopWatch.h"
#include <iostream>
#ifdef _OPENMP
#include <omp.h>
#endif


namespace rtengine {

extern const Settings *settings;

namespace {

//-----------------------------------------------------------------------------
// Copyright 2019 by Ingo Weyrich -- released under GPL
//-----------------------------------------------------------------------------
float calcRadiusBayer(const float * const *rawData, int W, int H, float lowerLimit, float upperLimit, const unsigned int fc[2])
{
    BENCHFUN
    float maxRatio = 1.f;
#ifdef _OPENMP
    #pragma omp parallel for reduction(max:maxRatio) schedule(dynamic, 16)
#endif
    for (int row = 4; row < H - 4; ++row) {
        for (int col = 5 + (fc[row & 1] & 1); col < W - 4; col += 2) {
            const float val00 = rawData[row][col];
            if (val00 > 0.f) {
                const float val1m1 = rawData[row + 1][col - 1];
                const float val1p1 = rawData[row + 1][col + 1];
                const float maxVal0 = std::max(val00, val1m1);
                if (val1m1 > 0.f && maxVal0 > lowerLimit) {
                    const float minVal = std::min(val00, val1m1);
                    if (UNLIKELY(maxVal0 > maxRatio * minVal)) {
                        bool clipped = false;
                        if (maxVal0 == val00) { // check for influence by clipped green in neighborhood
                            if (rtengine::max(rawData[row - 1][col - 1], rawData[row - 1][col + 1], val1p1) >= upperLimit) {
                                clipped = true;
                            }
                        } else { // check for influence by clipped green in neighborhood
                            if (rtengine::max(rawData[row][col - 2], val00, rawData[row + 2][col - 2], rawData[row + 2][col]) >= upperLimit) {
                                clipped = true;
                            }
                        }
                        if (!clipped) {
                            maxRatio = maxVal0 / minVal;
                        }
                    }
                }
                const float maxVal1 = std::max(val00, val1p1);
                if (val1p1 > 0.f && maxVal1 > lowerLimit) {
                    const float minVal = std::min(val00, val1p1);
                    if (UNLIKELY(maxVal1 > maxRatio * minVal)) {
                        if (maxVal1 == val00) { // check for influence by clipped green in neighborhood
                            if (rtengine::max(rawData[row - 1][col - 1], rawData[row - 1][col + 1], val1p1) >= upperLimit) {
                                continue;
                            }
                        } else { // check for influence by clipped green in neighborhood
                            if (rtengine::max(val00, rawData[row][col + 2], rawData[row + 2][col], rawData[row + 2][col + 2]) >= upperLimit) {
                                continue;
                            }
                        }
                        maxRatio = maxVal1 / minVal;
                    }
                }
            }
        }
    }
    if (settings->verbose) {
        std::cout << "Bayer auto deconv radius - maxRatio : " << maxRatio << std::endl;
    }
    float radius = std::sqrt((1.f / (std::log(1.f / maxRatio) /  2.f)) / -2.f);
    if (settings->verbose) {
        std::cout << "                           radius : " << radius << std::endl;
    }
    return radius;
}

} // namespace

bool RawImageSource::getDeconvAutoRadius(float *out)
{
    if (ri->getSensorType() == ST_BAYER) {
        if (!out) {
            return true; // only check whether this is supported
        }
        const float clipVal = (ri->get_white(1) - ri->get_cblack(1)) * scale_mul[1];
        const unsigned int fc[2] = {FC(0,0), FC(1,0)};
        *out = calcRadiusBayer(rawData, W, H, 1000.f, clipVal, fc);        
        return true;
    }
    return false;
}

} // namespace rtengine
