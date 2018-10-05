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
#include <iostream>

extern Options options;

namespace rtengine {

namespace {

#if 0
#  define DEBUG_DUMP(arr)                                                 \
    do {                                                                \
        Imagefloat im(arr.width(), arr.height());                      \
        const char *out = "/tmp/" #arr ".tif";                     \
        for (int y = 0; y < im.getHeight(); ++y) {                      \
            for (int x = 0; x < im.getWidth(); ++x) {                   \
                im.r(y, x) = im.g(y, x) = im.b(y, x) = arr[y][x] * 65535.f; \
            }                                                           \
        }                                                               \
        im.saveTIFF(out, 16);                                           \
    } while (false)
#else
#  define DEBUG_DUMP(arr)
#endif


std::pair<int, int> get_dark_channel(const Imagefloat &src, array2D<float> &dst,
                                     int patchsize, float *ambient=nullptr)
{
    const int w = src.getWidth();
    const int h = src.getHeight();

    float maxval = -RT_INFINITY_F;
    int maxx = 0, maxy = 0;

    for (int y = 0; y < src.getHeight(); y += patchsize) {
        int pH = std::min(y+patchsize, h);
        for (int x = 0; x < src.getWidth(); x += patchsize) {
            float val = RT_INFINITY_F;
            int pW = std::min(x+patchsize, w);
#ifdef _OPENMP
            #pragma omp parallel for shared(val)
#endif
            for (int yy = y; yy < pH; ++yy) {
                for (int xx = x; xx < pW; ++xx) {
                    float r = src.r(yy, xx);
                    float g = src.g(yy, xx);
                    float b = src.b(yy, xx);
                    if (ambient) {
                        r /= ambient[0];
                        g /= ambient[1];
                        b /= ambient[2];
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
            for (int yy = y; yy < pH; ++yy) {
                std::fill(dst[yy]+x, dst[yy]+pW, val);
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

    img->normalizeFloatTo1();
    
    const int W = img->getWidth();
    const int H = img->getHeight();
    const float strength = LIM01(float(params->dehaze.strength) / 100.f);

    if (options.rtSettings.verbose) {
        std::cout << "dehaze: strength = " << strength << std::endl;
    }
    
    array2D<float> dark(W, H);
    const int patchsize = std::max(W / 50, 2);
    auto p = get_dark_channel(*img, dark, patchsize);
    float ambient[3];
    float maxl = -RT_INFINITY_F;

    DEBUG_DUMP(dark);

    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
    
    const int pW = std::min(p.first+patchsize, W);
    const int pH = std::min(p.second+patchsize, H);

    if (options.rtSettings.verbose) {
        std::cout << "dehaze: computing ambient light from patch at "
                  << p.first << ", " << p.second << std::endl;
    }
    
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int y = p.second; y < pH; ++y) {
        for (int x = p.first; x < pW; ++x) {
            float r = img->r(y, x);
            float g = img->g(y, x);
            float b = img->b(y, x);
            float l = Color::rgbLuminance(r, g, b, ws);
            if (l > maxl) {
                ambient[0] = r;
                ambient[1] = g;
                ambient[2] = b;
                maxl = l;
            }
        }
    }

    if (options.rtSettings.verbose) {
        std::cout << "dehaze: ambient light is "
                  << ambient[0] << ", " << ambient[1] << ", " << ambient[2]
                  << std::endl;
    }
    
    if (min(ambient[0], ambient[1], ambient[2]) < 0.1f) {
        if (options.rtSettings.verbose) {
            std::cout << "dehaze: no haze detected" << std::endl;
        }
        img->normalizeFloatTo65535();
        return; // probably no haze at all
    }

    array2D<float> &t_tilde = dark;
    get_dark_channel(*img, dark, patchsize, ambient);
    DEBUG_DUMP(t_tilde);
    
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
            Y[y][x] = Color::rgbLuminance(img->r(y, x), img->g(y, x), img->b(y, x), ws);
        }
    }

    const int radius = patchsize * 4;
    const float epsilon = 0.05;
    array2D<float> &t = t_tilde;
    
    guidedFilter(Y, t_tilde, t, radius, epsilon, true);

    DEBUG_DUMP(t);

    const float t0 = 0.1;
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            float mt = std::max(t[y][x], t0);
            float r = (img->r(y, x) - ambient[0]) / mt + ambient[0];
            float g = (img->g(y, x) - ambient[1]) / mt + ambient[1];
            float b = (img->b(y, x) - ambient[2]) / mt + ambient[2];
            img->r(y, x) = r;
            img->g(y, x) = g;
            img->b(y, x) = b;
        }
    }

    img->normalizeFloatTo65535();
}


} // namespace rtengine
