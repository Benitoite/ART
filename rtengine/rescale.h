/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2018 Alberto Griggio <alberto.griggio@gmail.com>
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

namespace rtengine {

template <class T>
float getBilinearValue(T **src, int W, int H, float x, float y)
{
    // Get integer and fractional parts of numbers
    int xi = x;
    int yi = y;
    float xf = x - xi;
    float yf = y - yi;
    int xi1 = std::min(xi + 1, W - 1);
    int yi1 = std::min(yi + 1, H - 1);

    float bl = src[yi][xi];
    float br = src[yi][xi1];
    float tl = src[yi1][xi];
    float tr = src[yi1][xi1];

    // interpolate
    float b = xf * br + (1.f - xf) * bl;
    float t = xf * tr + (1.f - xf) * tl;
    float pxf = yf * t + (1.f - yf) * b;
    return pxf;
}


template <class T1, class T2>
void rescaleBilinear(T1 **src, int Ws, int Hs, T2 **dst, int Wd, int Hd, bool multithread)
{
    float col_scale = float (Ws) / float (Wd);
    float row_scale = float (Hs) / float (Hd);

#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif

    for (int y = 0; y < Hd; ++y) {
        float ymrs = y * row_scale;

        for (int x = 0; x < Wd; ++x) {
            dst[y][x] = getBilinearValue(src, Ws, Hs, x * col_scale, ymrs);
        }
    }
}


template <class T1, class T2>
void rescaleNearest(T1 **src, int Ws, int Hs, T2 **dst, int Wd, int Hd, bool multithread)
{
    const int width = Ws;
    const int height = Hs;
    const int nw = Wd;
    const int nh = Hd;

#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif

    for (int y = 0; y < nh; ++y) {
        int sy = y * height / nh;

        for (int x = 0; x < nw; ++x) {
            int sx = x * width / nw;
            dst[y][x] = src[sy][sx];
        }
    }
}


} // namespace rtengine
