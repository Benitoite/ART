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


namespace {

void sharpenHaloCtrl (float** luminance, float** blurmap, float** base, float** blend, int W, int H, const SharpeningParams &sharpenParam)
{

    const float scale = (100.f - sharpenParam.halocontrol_amount) * 0.01f;
    const float sharpFac = sharpenParam.amount * 0.01f;
    float** nL = base;

#ifdef _OPENMP
    #pragma omp parallel for
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

void dcdamping (float** aI, float** aO, float damping, int W, int H)
{

    const float dampingFac = -2.0 / (damping * damping);

#ifdef __SSE2__
    vfloat Iv, Ov, Uv, zerov, onev, fourv, fivev, dampingFacv, Tv, Wv, Lv;
    zerov = _mm_setzero_ps();
    onev = F2V(1.f);
    fourv = F2V(4.f);
    fivev = F2V(5.f);
    dampingFacv = F2V(dampingFac);
#endif
#ifdef _OPENMP
    #pragma omp for
#endif

    for (int i = 0; i < H; i++) {
        int j = 0;
#ifdef __SSE2__

        for (; j < W - 3; j += 4) {
            Iv = LVFU(aI[i][j]);
            Ov = LVFU(aO[i][j]);
            Lv = xlogf(Iv / Ov);
            Wv = Ov - Iv;
            Uv = (Ov * Lv + Wv) * dampingFacv;
            Uv = vminf(Uv, onev);
            Tv = Uv * Uv;
            Tv = Tv * Tv;
            Uv = Tv * (fivev - Uv * fourv);
            Uv = (Wv / Iv) * Uv + onev;
            Uv = vselfzero(vmaskf_gt(Iv, zerov), Uv);
            Uv = vselfzero(vmaskf_gt(Ov, zerov), Uv);
            STVFU(aI[i][j], Uv);
        }

#endif

        for(; j < W; j++) {
            float I = aI[i][j];
            float O = aO[i][j];

            if (O <= 0.f || I <= 0.f) {
                aI[i][j] = 0.f;
                continue;
            }

            float U = (O * xlogf(I / O) - I + O) * dampingFac;
            U = rtengine::min(U, 1.0f);
            U = U * U * U * U * (5.f - U * 4.f);
            aI[i][j] = (O - I) / I * U + 1.f;
        }
    }
}

void deconvsharpening(float** luminance, float ** L, float** tmp, int W, int H, const SharpeningParams &sharpenParam, float contrast, double scale)
{
    const auto blurradius = sharpenParam.blurradius / scale;
    if (sharpenParam.deconvamount == 0 && blurradius < 0.25f) {
        return;
    }
BENCHFUN
    JaggedArray<float> tmpI(W, H);

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int i = 0; i < H; i++) {
        for(int j = 0; j < W; j++) {
            tmpI[i][j] = std::max(luminance[i][j], 0.f);
        }
    }

    // calculate contrast based blend factors to reduce sharpening in regions with low contrast
    JaggedArray<float> blend(W, H);
    //float contrast = sharpenParam.contrast / 100.f;
    buildBlendMask(L, blend, W, H, contrast, 1.f);
    JaggedArray<float>* blurbuffer = nullptr;

    if (blurradius >= 0.25f) {
        blurbuffer = new JaggedArray<float>(W, H);
        JaggedArray<float> &blur = *blurbuffer;
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            gaussianBlur(tmpI, blur, W, H, blurradius);
#ifdef _OPENMP
            #pragma omp for
#endif
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    blur[i][j] = intp(blend[i][j], luminance[i][j], std::max(blur[i][j], 0.0f));
                }
            }
        }
    }
    const float damping = sharpenParam.deconvdamping / 5.0;
    const bool needdamp = sharpenParam.deconvdamping > 0;
    const double sigma = sharpenParam.deconvradius / scale;
    const float amount = sharpenParam.deconvamount / 100.f;

    JaggedArray<char> impulse(W, H);
    markImpulse(W, H, luminance, impulse, 2.f);

#ifdef _OPENMP
    #pragma omp parallel
#endif
    {
        for (int k = 0; k < sharpenParam.deconviter; k++) {
            if (!needdamp) {
                // apply gaussian blur and divide luminance by result of gaussian blur
                gaussianBlur(tmpI, tmp, W, H, sigma, nullptr, GAUSS_DIV, luminance);
            } else {
                // apply gaussian blur + damping
                gaussianBlur(tmpI, tmp, W, H, sigma);
                dcdamping(tmp, luminance, damping, W, H);
            }
            gaussianBlur(tmp, tmpI, W, H, sigma, nullptr, GAUSS_MULT);
        } // end for

#ifdef _OPENMP
        #pragma omp for
#endif

        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                float b = impulse[i][j] ? 0.f : blend[i][j] * amount;
                luminance[i][j] = intp(b, std::max(tmpI[i][j], 0.0f), luminance[i][j]);
            }
        }

        if (blurradius >= 0.25f) {
            JaggedArray<float> &blur = *blurbuffer;
#ifdef _OPENMP
        #pragma omp for
#endif
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    luminance[i][j] = intp(blend[i][j], luminance[i][j], std::max(blur[i][j], 0.0f));
                }
            }
        }
    } // end parallel
    delete blurbuffer;
}

} // namespace

namespace rtengine {

extern const Settings* settings;

void ImProcFunctions::sharpening(Imagefloat *rgb, const SharpeningParams &sharpenParam, bool showMask)
{
    LabImage tmplab(rgb->getWidth(), rgb->getHeight());
    rgb2lab(*rgb, tmplab);
    LabImage *lab = &tmplab;

    if ((!sharpenParam.enabled) || sharpenParam.amount < 1 || lab->W < 8 || lab->H < 8) {
        return;
    }

    int W = lab->W, H = lab->H;
    float contrast = pow_F(sharpenParam.contrast / 100.f, 1.2f);

    if(showMask) {
        // calculate contrast based blend factors to reduce sharpening in regions with low contrast
        JaggedArray<float> blend(W, H);
        buildBlendMask(lab->L, blend, W, H, contrast, 1.f);
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                lab->L[i][j] = blend[i][j] * 32768.f;
            }
        }

        lab2rgb(*lab, *rgb);
        return;
    }

    JaggedArray<float> b2(W, H);

    if (sharpenParam.method == "rld") {
        // Imagefloat rgb(W, H);
        JaggedArray<float> Y(W, H);
        auto ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
        // lab2rgb(*lab, rgb, params->icm.workingProfile);
#ifdef _OPENMP
#       pragma omp parallel for
#endif
        for (int y = 0; y < lab->H; ++y) {
            for (int x = 0; x < lab->W; ++x) {
                Y[y][x] = Color::rgbLuminance(rgb->r(y, x), rgb->g(y, x), rgb->b(y, x), ws);
            }
        }
        deconvsharpening(Y, lab->L, b2, W, H, sharpenParam, contrast, scale);
#ifdef _OPENMP
#       pragma omp parallel for
#endif
        for (int y = 0; y < lab->H; ++y) {
            for (int x = 0; x < lab->W; ++x) {
                float oY = Color::rgbLuminance(rgb->r(y, x), rgb->g(y, x), rgb->b(y, x), ws);
                if (oY > 1e-5f) {
                    float f = Y[y][x] / oY;
                    rgb->r(y, x) *= f;
                    rgb->g(y, x) *= f;
                    rgb->b(y, x) *= f;
                }
            }
        }
        //rgb2lab(rgb, *lab, params->icm.workingProfile);
        return;
    }
BENCHFUN

    // Rest is UNSHARP MASK
    float** b3 = nullptr;

    if (sharpenParam.edgesonly) {
        b3 = new float*[H];

        for (int i = 0; i < H; i++) {
            b3[i] = new float[W];
        }
    }

    // calculate contrast based blend factors to reduce sharpening in regions with low contrast
    JaggedArray<float> blend(W, H);
    //float contrast = sharpenParam.contrast / 100.f;
    buildBlendMask(lab->L, blend, W, H, contrast);

    JaggedArray<float> blur(W, H);

    const auto blurradius = sharpenParam.blurradius / scale;
    if (blurradius >= 0.25f) {
#ifdef _OPENMP
        #pragma omp parallel
#endif
        {
            gaussianBlur(lab->L, blur, W, H, blurradius);
#ifdef _OPENMP
            #pragma omp for
#endif
            for (int i = 0; i < H; ++i) {
                for (int j = 0; j < W; ++j) {
                    blur[i][j] = intp(blend[i][j], lab->L[i][j], std::max(blur[i][j], 0.0f));
                }
            }
        }
    }


#ifdef _OPENMP
    #pragma omp parallel
#endif
    {

        if (!sharpenParam.edgesonly) {
            gaussianBlur (lab->L, b2, W, H, sharpenParam.radius / scale);
        } else {
            bilateral<float, float> (lab->L, (float**)b3, b2, W, H, sharpenParam.edges_radius / scale, sharpenParam.edges_tolerance, multiThread);
            gaussianBlur (b3, b2, W, H, sharpenParam.radius / scale);
        }
    }
    float** base = lab->L;

    if (sharpenParam.edgesonly) {
        base = b3;
    }

    if (!sharpenParam.halocontrol) {
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = 0; i < H; i++)
            for (int j = 0; j < W; j++) {
                constexpr float upperBound = 2000.f;  // WARNING: Duplicated value, it's baaaaaad !
                float diff = base[i][j] - b2[i][j];
                float delta = sharpenParam.threshold.multiply<float, float, float>(
                                  min(fabsf(diff), upperBound),                   // X axis value = absolute value of the difference, truncated to the max value of this field
                                  sharpenParam.amount * diff * 0.01f        // Y axis max value
                              );
                lab->L[i][j] = intp(blend[i][j], lab->L[i][j] + delta, lab->L[i][j]);
            }
    } else {
        if (!sharpenParam.edgesonly) {
            // make a deep copy of lab->L
            JaggedArray<float> labCopy(W, H);

#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for( int i = 0; i < H; i++ )
                for( int j = 0; j < W; j++ ) {
                    labCopy[i][j] = lab->L[i][j];
                }

            sharpenHaloCtrl (lab->L, b2, labCopy, blend, W, H, sharpenParam);
        } else {
            sharpenHaloCtrl (lab->L, b2, base, blend, W, H, sharpenParam);
        }

    }

    if (sharpenParam.edgesonly) {
        for (int i = 0; i < H; i++) {
            delete [] b3[i];
        }

        delete [] b3;
    }

    if (blurradius >= 0.25f) {
#ifdef _OPENMP
    #pragma omp parallel for
#endif
        for (int i = 0; i < H; ++i) {
            for (int j = 0; j < W; ++j) {
                lab->L[i][j] = intp(blend[i][j], lab->L[i][j], std::max(blur[i][j], 0.0f));
            }
        }
    }

    lab2rgb(*lab, *rgb);
}

} // namespace rtengine
