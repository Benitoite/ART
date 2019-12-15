/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
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

#include "improcfun.h"
#include "gauss.h"
#include "bilateral2.h"
#include "jaggedarray.h"
#include "rt_math.h"
#include "sleef.c"
#include "opthelper.h"
//#define BENCHMARK
#include "StopWatch.h"
#include "rt_algo.h"
using namespace std;


namespace rtengine { extern const Settings *settings; }


namespace {

template <bool reverse, bool usm>
void apply_gamma(float **Y, int W, int H, float pivot, float gamma, bool multiThread)
{
    BENCHFUN

    if (!reverse) {
        gamma = 1.f/gamma;
    }

    LUTf glut(65536);
    glut[0] = 0;
    const float d = 65535.f * pivot;
    for (int i = 1; i < 65536; ++i) {
        glut[i] = pow_F(float(i)/d, gamma) * pivot;
        if (reverse || usm) {
            glut[i] *= 65535.f;
        }
    }
        
#ifdef _OPENMP
#    pragma omp parallel for if (multiThread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            float l = Y[y][x];
            if (LIKELY(l >= 0.f && l < 65536.f)) {
                if (reverse && !usm) {
                    l *= 65535.f;
                }
                l = glut[l];
            } else {
                l = pow_F(std::max(l / d, 1e-18f), gamma) * pivot;
                if (reverse || usm) {
                    l *= 65535.f;
                }
            }
            Y[y][x] = l;
        }
    }
}


void sharpenHaloCtrl(float** luminance, float** blurmap, float** base, float** blend, int W, int H, const SharpeningParams &sharpenParam, bool multiThread)
{

    const float scale = (100.f - sharpenParam.halocontrol_amount) * 0.01f;
    const float sharpFac = sharpenParam.amount * 0.01f;
    float** nL = base;

#ifdef _OPENMP
#   pragma omp parallel for if (multiThread)
#endif

    for (int i = 2; i < H - 2; i++) {
        float max1 = 0, max2 = 0, min1 = 0, min2 = 0;

        for (int j = 2; j < W - 2; j++) {
            // compute 3 iterations, only forward
            float np1 = 2.f * (nL[i - 2][j] + nL[i - 2][j + 1] + nL[i - 2][j + 2] + nL[i - 1][j] + nL[i - 1][j + 1] + nL[i - 1][j + 2] + nL[i]  [j] + nL[i]  [j + 1] + nL[i]  [j + 2]) / 27.f + nL[i - 1][j + 1] / 3.f;
            float np2 = 2.f * (nL[i - 1][j] + nL[i - 1][j + 1] + nL[i - 1][j + 2] + nL[i]  [j] + nL[i]  [j + 1] + nL[i]  [j + 2] + nL[i + 1][j] + nL[i + 1][j + 1] + nL[i + 1][j + 2]) / 27.f + nL[i]  [j + 1] / 3.f;
            float np3 = 2.f * (nL[i]  [j] + nL[i]  [j + 1] + nL[i]  [j + 2] + nL[i + 1][j] + nL[i + 1][j + 1] + nL[i + 1][j + 2] + nL[i + 2][j] + nL[i + 2][j + 1] + nL[i + 2][j + 2]) / 27.f + nL[i + 1][j + 1] / 3.f;

            // Max/Min of all these deltas and the last two max/min
            float maxn = rtengine::max(np1, np2, np3);
            float minn = rtengine::min(np1, np2, np3);
            float max_ = rtengine::max(max1, max2, maxn);
            float min_ = rtengine::min(min1, min2, minn);

            // Shift the queue
            max1 = max2;
            max2 = maxn;
            min1 = min2;
            min2 = minn;
            float labL = luminance[i][j];

            if (max_ < labL) {
                max_ = labL;
            }

            if (min_ > labL) {
                min_ = labL;
            }

            // deviation from the environment as measurement
            float diff = nL[i][j] - blurmap[i][j];

            constexpr float upperBound = 2000.f;  // WARNING: Duplicated value, it's baaaaaad !
            float delta = sharpenParam.threshold.multiply<float, float, float>(
                              rtengine::min(fabsf(diff), upperBound),   // X axis value = absolute value of the difference
                              sharpFac * diff               // Y axis max value = sharpening.amount * signed difference
                          );
            float newL = labL + delta;

            // applying halo control
            if (newL > max_) {
                newL = max_ + (newL - max_) * scale;
            } else if (newL < min_) {
                newL = min_ - (min_ - newL) * scale;
            }

            luminance[i][j] = intp(blend[i][j], newL, luminance[i][j]);
        }
    }
}


void deconvsharpening(float **luminance, float **blend, char **impulse, int W, int H, const SharpeningParams &sharpenParam, double scale, bool multiThread)
{
    if (sharpenParam.deconvamount == 0) {
        return;
    }
BENCHFUN
    //apply_gamma<false, false>(luminance, W, H, 0.18f, 3.f, multiThread);

    JaggedArray<float> tmp(W, H);
    JaggedArray<float> tmpI(W, H);
    JaggedArray<float> out(W, H);

#ifdef _OPENMP
#   pragma omp parallel for if (multiThread)
#endif
    for (int i = 0; i < H; i++) {
        for(int j = 0; j < W; j++) {
            tmpI[i][j] = std::max(luminance[i][j], 0.f);
            out[i][j] = RT_NAN;
        }
    }

    const double sigma = sharpenParam.deconvradius / scale;
    const float amount = sharpenParam.deconvamount / 100.f;
    const int maxiter = 20;
    const float delta_factor = 0.2f;
    
    const auto get_output =
        [&](int i, int j) -> float
        {
            float b = impulse[i][j] ? 0.f : blend[i][j] * amount;
            return intp(b, std::max(tmpI[i][j], 0.0f), luminance[i][j]);
        };

    const auto check_stop =
        [&](int y, int x) -> void
        {
            if (LIKELY(std::isnan(out[y][x]))) {
                float l = luminance[y][x];
                float delta = l * delta_factor;
                if (UNLIKELY(std::abs(tmpI[y][x] - l) > delta)) {
                    out[y][x] = get_output(y, x);
                }
            }
        };

#ifdef _OPENMP
#   pragma omp parallel if (multiThread)
#endif
    {
        for (int k = 0; k < maxiter; k++) {
            gaussianBlur(tmpI, tmp, W, H, sigma, nullptr, GAUSS_DIV, luminance);
            gaussianBlur(tmp, tmpI, W, H, sigma, nullptr, GAUSS_MULT);
#ifdef _OPENMP
#           pragma omp for
#endif
            for (int y = 0; y < H; ++y) {
                for (int x = 0; x < W; ++x) {
                    check_stop(y, x);
                }
            }
        }

#ifdef _OPENMP
#       pragma omp for
#endif
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                float l = out[i][j];
                if (std::isnan(l)) {
                    l = get_output(i, j);
                }
                luminance[i][j] = l;
            }
        }
    }

    //apply_gamma<true, false>(luminance, W, H, 0.18f, 3.f, multiThread);    
}


void unsharp_mask(float **Y, float **blend, int W, int H, const SharpeningParams &sharpenParam, double scale, bool multiThread)
{
BENCHFUN
    apply_gamma<false, true>(Y, W, H, 1.f, 3.f, multiThread);

    float** b3 = nullptr;

    if (sharpenParam.edgesonly) {
        b3 = new float*[H];

        for (int i = 0; i < H; i++) {
            b3[i] = new float[W];
        }
    }

    JaggedArray<float> b2(W, H);

#ifdef _OPENMP
#   pragma omp parallel if (multiThread)
#endif
    {

        if (!sharpenParam.edgesonly) {
            gaussianBlur(Y, b2, W, H, sharpenParam.radius / scale);
        } else {
            bilateral<float, float>(Y, (float**)b3, b2, W, H, sharpenParam.edges_radius / scale, sharpenParam.edges_tolerance, multiThread);
            gaussianBlur (b3, b2, W, H, sharpenParam.radius / scale);
        }
    }
    float** base = Y;

    if (sharpenParam.edgesonly) {
        base = b3;
    }

    if (!sharpenParam.halocontrol) {
#ifdef _OPENMP
#       pragma omp parallel for if (multiThread)
#endif

        for (int i = 0; i < H; i++)
            for (int j = 0; j < W; j++) {
                constexpr float upperBound = 2000.f;  // WARNING: Duplicated value, it's baaaaaad !
                float diff = base[i][j] - b2[i][j];
                float delta = sharpenParam.threshold.multiply<float, float, float>(
                                  std::min(fabsf(diff), upperBound),                   // X axis value = absolute value of the difference, truncated to the max value of this field
                                  sharpenParam.amount * diff * 0.01f        // Y axis max value
                              );
                Y[i][j] = intp(blend[i][j], Y[i][j] + delta, Y[i][j]);
            }
    } else {
        if (!sharpenParam.edgesonly) {
            // make a deep copy of lab->L
            JaggedArray<float> labCopy(W, H);

#ifdef _OPENMP
#           pragma omp parallel for if (multiThread)
#endif

            for( int i = 0; i < H; i++ )
                for( int j = 0; j < W; j++ ) {
                    labCopy[i][j] = Y[i][j];
                }

            sharpenHaloCtrl(Y, b2, labCopy, blend, W, H, sharpenParam, multiThread);
        } else {
            sharpenHaloCtrl(Y, b2, base, blend, W, H, sharpenParam, multiThread);
        }

    }

    if (sharpenParam.edgesonly) {
        for (int i = 0; i < H; i++) {
            delete [] b3[i];
        }

        delete [] b3;
    }

    apply_gamma<true, true>(Y, W, H, 1.f, 3.f, multiThread);
}


} // namespace

namespace rtengine {

bool ImProcFunctions::sharpening(Imagefloat *rgb, const SharpeningParams &sharpenParam, bool showMask)
{
    const int W = rgb->getWidth();
    const int H = rgb->getHeight();

    if ((!sharpenParam.enabled) || sharpenParam.amount < 1 || W < 8 || H < 8) {
        return false;
    }

    rgb->setMode(Imagefloat::Mode::YUV, multiThread);
    float **Y = rgb->g.ptrs;
    float s_scale = std::sqrt(scale);
    float contrast = pow_F(sharpenParam.contrast / 100.f, 1.2f) * s_scale;
    JaggedArray<float> blend(W, H);
    buildBlendMask(Y, blend, W, H, contrast, 1.f, false, 2.f / s_scale);
    
    if (showMask) {
#ifdef _OPENMP
#       pragma omp parallel for if (multiThread)
#endif
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                Y[i][j] = blend[i][j] * 65536.f;
            }
        }
        return true;
    }

    std::unique_ptr<JaggedArray<char>> impulse;
    if (sharpenParam.method == "rld") {
        impulse.reset(new JaggedArray<char>(W, H));
        markImpulse(W, H, Y, *impulse, 2.f);
    }
    
    
    if (sharpenParam.method == "rld") {
        deconvsharpening(Y, blend, *impulse, W, H, sharpenParam, scale, multiThread);
    } else {
        unsharp_mask(Y, blend, W, H, sharpenParam, scale, multiThread);
    }


    return false;
}

} // namespace rtengine
