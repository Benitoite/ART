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

/*
 * Haze removal using the algorithm described in the paper:
 *
 * Single Image Haze Removal Using Dark Channel Prior
 * by He, Sun and Tang
 *
 * using a guided filter for the "soft matting" of the transmission map
 *
 */  

#include "improcfun.h"
#include "guidedfilter.h"
#include "rt_math.h"

namespace rtengine {

namespace {

std::pair<int, int> get_dark_channel(const Imagefloat &src, array2D<float> &dst,
                                     int patchsize, float *scale=nullptr)
{
    const int w = src.getWidth();
    const int h = src.getHeight();

    float maxval = -RT_INFINITY_F;
    int maxx = 0, maxy = 0;

    for (int y = 0; y < src.getHeight(); y += patchsize) {
        for (int x = 0; x < src.getWidth(); x += patchsize) {
            float val = RT_INFINITY_F;
            int ysz = std::min(y+patchsize, h) - y;
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for (int yy = 0; yy < ysz; ++yy) {
                for (int xx = 0; x+xx < std::min(x+patchsize, w); ++xx) {
                    float r = src.r(y+yy, x+xx) / 65535.f;
                    float g = src.g(y+yy, x+xx) / 65535.f;
                    float b = src.b(y+yy, x+xx) / 65535.f;
                    if (scale) {
                        r /= scale[0];
                        g /= scale[1];
                        b /= scale[2];
                    }
                    val = min(val, r, g, b);
                }
            }
            if (val > maxval) {
                maxval = val;
                maxx = x;
                maxy = y;
            }
#ifdef _OPENMP
            #pragma omp parallel for
#endif
            for (int yy = 0; yy < ysz; ++yy) {
                for (int xx = 0; x+xx < std::min(x+patchsize, w); ++xx) {
                    dst[y+yy][x+xx] = val;
                }
            }
        }
    }

    return std::make_pair(maxx, maxy);
}

} // namespace


void ImProcFunctions::dehaze(Imagefloat *img)
{
    if (!params->dehaze.enabled) {
        return;
    }
    
    const int W = img->getWidth();
    const int H = img->getHeight();
    const float strength = LIM01(float(params->dehaze.strength) / 100.f);
    
    array2D<float> dark(W, H);
    const int patchsize = W / 50;
    auto p = get_dark_channel(*img, dark, patchsize);
    float ambient[3];
    float maxl = -RT_INFINITY_F;

    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
    

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int y = p.second; y < std::min(p.second+patchsize, H); ++y) {
        for (int x = p.first; x < std::min(p.first+patchsize, H); ++x) {
            float r = img->r(y, x) / 65535.f;
            float g = img->g(y, x) / 65535.f;
            float b = img->b(y, x) / 65535.f;
            float l = Color::rgbLuminance(r, g, b, ws);
            if (l > maxl) {
                ambient[0] = r;
                ambient[1] = g;
                ambient[2] = b;
                maxl = l;
            }
        }
    }

    get_dark_channel(*img, dark, patchsize, ambient);
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            dark[y][x] = 1.f - strength * dark[y][x];
        }
    }

    array2D<float> Y(W, H);
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            Y[y][x] = Color::rgbLuminance(img->r(y, x), img->g(y, x), img->b(y, x), ws) / 65535.f;
        }
    }

    const int radius = patchsize * 4;
    const float epsilon = 0.05;
    array2D<float> &t = dark;
    guidedFilter(Y, t, t, radius, epsilon, true);

    const float t0 = 0.1;
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            float mt = std::max(t[y][x], t0);
            float r = (img->r(y, x) / 65535.f - ambient[0]) / mt + ambient[0];
            float g = (img->g(y, x) / 65535.f - ambient[1]) / mt + ambient[1];
            float b = (img->b(y, x) / 65535.f - ambient[2]) / mt + ambient[2];
            img->r(y, x) = r * 65535.f;
            img->g(y, x) = g * 65535.f;
            img->b(y, x) = b * 65535.f;
        }
    }
}


} // namespace rtengine
