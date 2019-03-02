/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright 2018 Alberto Griggio <alberto.griggio@gmail.com>
 *
 *  
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

// adapted from the tone equalizer of darktable
/*
    This file is part of darktable,
    copyright (c) 2018 Aurelien Pierre.

    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
*/    

#include <array>
#include "improcfun.h"
#include "gauss.h"
#include "sleef.c"
#include "opthelper.h"
#include "guidedfilter.h"
#define BENCHMARK
#include "StopWatch.h"

namespace rtengine {

namespace {

typedef std::array<float, 5> KernelRow;
typedef std::array<KernelRow, 5> Kernel;

void laplacian_filter(const array2D<float> &Yorig, const array2D<float> &Ytoned,
                      array2D<float> &out, const Kernel &kernel, bool multithread)
{
    const int padding = floor(kernel.size() / 2);
    const int H = Yorig.height();
    const int W = Yorig.width();

    const auto padidx =
        [padding](int i, int o, int l) -> int
        {
            int ret = i + o - padding;
            return LIM(ret, 0, l-1);
        };

#ifdef _OPENMP
#    pragma omp parallel for if (multithread)
#endif
    for(int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            float orig = LIM01(Yorig[y][x]);
            float toned = LIM01(Ytoned[y][x]);

            float weight = 0.0f;

            // Convolution filter
            for (size_t m = 0; m < kernel.size(); ++m) {
                int yy = padidx(y, m, H);
                for (size_t n = 0; n < kernel.size(); ++n) {
                    int xx = padidx(x, n, W);
                    float neighbour_orig = LIM01(Yorig[yy][xx]);
                    float neighbour_toned = LIM01(Ytoned[yy][xx]);
                    weight += (neighbour_toned - toned) * (neighbour_orig - orig) * kernel[m][n];
                }
            }

            // Corrective ratio to apply on the toned image
            out[y][x] =  toned / (toned + weight);
        }
    }
}


void tone_eq(array2D<float> &R, array2D<float> &G, array2D<float> &B, array2D<float> &Y, array2D<float> &Yout, const ToneEqualizerParams &pp, const Glib::ustring &workingProfile, double scale, bool multithread)
{
    const int W = R.width();
    const int H = R.height();
    // array2D<float> Y(W, H);

    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(workingProfile);

#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            Y[y][x] = Color::rgbLuminance(R[y][x], G[y][x], B[y][x], ws);
        }
    }

    // const float detail = float(pp.detail) / 20.f;
    // const int r = 15 / scale * (detail >= 0.f ? 1.f + detail : 1.f / (1.f - detail));
    // const float epsilon = 0.035f * std::pow(2.f, detail);
    // if (r > 0) {
    //     rtengine::guidedFilter(Y, Y, Y, r, epsilon, multithread);
    // }

    const auto log2 =
        [](float x) -> float
        {
            static const float l2 = xlogf(2);
            return xlogf(x) / l2;
        };

    const auto exp2 =
        [](float x) -> float
        {
            return pow_F(2.f, x);
        };

     // Build the luma channels: band-pass filters with gaussian windows of
     // std 2 EV, spaced by 2 EV
    const float centers[12] = {
        -18.0f, -16.0f, -14.0f, -12.0f, -10.0f, -8.0f, -6.0f,
        -4.0f, -2.0f, 0.0f, 2.0f, 4.0f
    };

    const auto conv = [&](int v) -> float
                      {
                          return exp2(float(v) / 100.f * 2.f);
                      };

    const float factors[12] = {
        conv(pp.bands[0]), // -18 EV
        conv(pp.bands[0]), // -16 EV
        conv(pp.bands[0]), // -14 EV
        conv(pp.bands[0]), // -12 EV
        conv(pp.bands[0]), // -10 EV
        conv(pp.bands[0]), //  -8 EV
        conv(pp.bands[1]), //  -6 EV
        conv(pp.bands[2]), //  -4 EV
        conv(pp.bands[3]), //  -2 EV
        conv(pp.bands[4]), //   0 EV
        conv(pp.bands[4]), //   2 EV
        conv(pp.bands[4])  //   4 EV
    };

    const auto gauss =
        [](float b, float x) -> float
        {
            return xexpf((-SQR(x - b) / 4.0f));
        };

    // For every pixel luminance, the sum of the gaussian masks
    // evenly spaced by 2 EV with 2 EV std should be this constant
    float w_sum = 0.f;
    for (int i = 0; i < 12; ++i) {
        w_sum += gauss(centers[i], 0.f);
    }

    const auto process_pixel =
        [&](float y) -> float
        {
            // get the luminance of the pixel - just average channels
            const float luma = max(log2(max(y, 0.f)), -18.0f);

            // build the correction as the sum of the contribution of each
            // luminance channel to current pixel
            float correction = 0.0f;
            for (int c = 0; c < 12; ++c) {
                correction += gauss(centers[c], luma) * factors[c];
            }
            correction /= w_sum;

            return correction;
        };

#ifdef __SSE2__
    vfloat vfactors[12];
    vfloat vcenters[12];
    
    for (int i = 0; i < 12; ++i) {
        vfactors[i] = F2V(factors[i]);
        vcenters[i] = F2V(centers[i]);
    }

    const auto vgauss =
        [](vfloat b, vfloat x) -> vfloat
        {
            static const vfloat fourv = F2V(4.f);
            return xexpf((-SQR(x - b) / fourv));
        };

    vfloat zerov = F2V(0.f);
    
    vfloat vw_sum = zerov;
    for (int i = 0; i < 12; ++i) {
        vw_sum += vgauss(vcenters[i], zerov);
    }

    const vfloat noisev = F2V(-18.f);
    const vfloat xlog2v = F2V(xlogf(2.f));
    
    const auto vprocess_pixel =
        [&](vfloat y) -> vfloat
        {
            // get the luminance of the pixel - just average channels
            const vfloat luma = vmaxf(xlogf(vmaxf(y, zerov))/xlog2v, noisev);

            // build the correction as the sum of the contribution of each
            // luminance channel to current pixel
            vfloat correction = zerov;
            for (int c = 0; c < 12; ++c) {
                correction += vgauss(vcenters[c], luma) * vfactors[c];
            }
            correction /= vw_sum;

            return correction;
        };
#endif // __SSE2__
    
        
#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        int x = 0;
#ifdef __SSE2__
        for (; x < W - 3; x += 4) {
            vfloat cY = LVFU(Y[y][x]);
            vfloat corr = vprocess_pixel(cY);
            STVFU(R[y][x], LVFU(R[y][x]) * corr);
            STVFU(G[y][x], LVFU(G[y][x]) * corr);
            STVFU(B[y][x], LVFU(B[y][x]) * corr);
            STVFU(Yout[y][x], cY * corr);
        }
#endif // __SSE2__
        for (; x < W; ++x) {
            float cY = Y[y][x];
            float corr = process_pixel(cY);
            R[y][x] *= corr;
            G[y][x] *= corr;
            B[y][x] *= corr;
            Yout[y][x] = cY * corr;
        }
        // for (x = 0; x < W; ++x) {
        //     Yout[y][x] = Color::rgbLuminance(R[y][x], G[y][x], B[y][x], ws);
        // }
    }
}


#if 1
void tone_eq_gf(array2D<float> &R, array2D<float> &G, array2D<float> &B, const ToneEqualizerParams &pp, const Glib::ustring &workingProfile, double scale, bool multithread)
{
    const int W = R.width();
    const int H = R.height();
    array2D<float> Y(W, H);

    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(workingProfile);

#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            Y[y][x] = Color::rgbLuminance(R[y][x], G[y][x], B[y][x], ws);
        }
    }

    const int radius = float(pp.detail) / scale + 0.5f;
    const float epsilon = 0.005f + 0.001f * max(pp.detail - 3, 0);
    if (radius > 0) {
        rtengine::guidedFilter(Y, Y, Y, radius, epsilon, multithread);
    }
    
    const auto log2 =
        [](float x) -> float
        {
            static const float l2 = xlogf(2);
            return xlogf(x) / l2;
        };

    const auto exp2 =
        [](float x) -> float
        {
            return pow_F(2.f, x);
        };

     // Build the luma channels: band-pass filters with gaussian windows of
     // std 2 EV, spaced by 2 EV
    const float centers[12] = {
        -18.0f, -16.0f, -14.0f, -12.0f, -10.0f, -8.0f, -6.0f,
        -4.0f, -2.0f, 0.0f, 2.0f, 4.0f
    };

    const auto conv = [&](int v) -> float
                      {
                          return exp2(float(v) / 100.f * 2.f);
                      };

    const float factors[12] = {
        conv(pp.bands[0]), // -18 EV
        conv(pp.bands[0]), // -16 EV
        conv(pp.bands[0]), // -14 EV
        conv(pp.bands[0]), // -12 EV
        conv(pp.bands[0]), // -10 EV
        conv(pp.bands[0]), //  -8 EV
        conv(pp.bands[1]), //  -6 EV
        conv(pp.bands[2]), //  -4 EV
        conv(pp.bands[3]), //  -2 EV
        conv(pp.bands[4]), //   0 EV
        conv(pp.bands[4]), //   2 EV
        conv(pp.bands[4])  //   4 EV
    };

    const auto gauss =
        [](float b, float x) -> float
        {
            return xexpf((-SQR(x - b) / 4.0f));
        };

    // For every pixel luminance, the sum of the gaussian masks
    // evenly spaced by 2 EV with 2 EV std should be this constant
    float w_sum = 0.f;
    for (int i = 0; i < 12; ++i) {
        w_sum += gauss(centers[i], 0.f);
    }

    const auto process_pixel =
        [&](float y) -> float
        {
            // get the luminance of the pixel - just average channels
            const float luma = max(log2(max(y, 0.f)), -18.0f);

            // build the correction as the sum of the contribution of each
            // luminance channel to current pixel
            float correction = 0.0f;
            for (int c = 0; c < 12; ++c) {
                correction += gauss(centers[c], luma) * factors[c];
            }
            correction /= w_sum;

            return correction;
        };

#ifdef __SSE2__
    vfloat vfactors[12];
    vfloat vcenters[12];
    
    for (int i = 0; i < 12; ++i) {
        vfactors[i] = F2V(factors[i]);
        vcenters[i] = F2V(centers[i]);
    }

    const auto vgauss =
        [](vfloat b, vfloat x) -> vfloat
        {
            static const vfloat fourv = F2V(4.f);
            return xexpf((-SQR(x - b) / fourv));
        };

    vfloat zerov = F2V(0.f);
    vfloat vw_sum = F2V(w_sum);

    const vfloat noisev = F2V(-18.f);
    const vfloat xlog2v = F2V(xlogf(2.f));
    
    const auto vprocess_pixel =
        [&](vfloat y) -> vfloat
        {
            // get the luminance of the pixel - just average channels
            const vfloat luma = vmaxf(xlogf(vmaxf(y, zerov))/xlog2v, noisev);

            // build the correction as the sum of the contribution of each
            // luminance channel to current pixel
            vfloat correction = zerov;
            for (int c = 0; c < 12; ++c) {
                correction += vgauss(vcenters[c], luma) * vfactors[c];
            }
            correction /= vw_sum;

            return correction;
        };
#endif // __SSE2__
    
        
#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        int x = 0;
#ifdef __SSE2__
        for (; x < W - 3; x += 4) {
            vfloat cY = LVFU(Y[y][x]);
            vfloat corr = vprocess_pixel(cY);
            STVFU(R[y][x], LVFU(R[y][x]) * corr);
            STVFU(G[y][x], LVFU(G[y][x]) * corr);
            STVFU(B[y][x], LVFU(B[y][x]) * corr);
        }
#endif // __SSE2__
        for (; x < W; ++x) {
            float cY = Y[y][x];
            float corr = process_pixel(cY);
            R[y][x] *= corr;
            G[y][x] *= corr;
            B[y][x] *= corr;
        }
    }
}

#else
void tone_eq_gf(array2D<float> &R, array2D<float> &G, array2D<float> &B, const ToneEqualizerParams &pp, const Glib::ustring &workingProfile, double scale, bool multithread)
{
    const int W = R.width();
    const int H = R.height();
    array2D<float> Y(W, H);

    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(workingProfile);

#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            Y[y][x] = Color::rgbLuminance(R[y][x], G[y][x], B[y][x], ws);
        }
    }

    int radius = float(pp.detail) / scale + 0.5f;
    const float epsilon = 0.005f + 0.001f * max(pp.detail - 3, 0);

    const auto log2 =
        [](float x) -> float
        {
            static const float l2 = xlogf(2);
            return xlogf(x) / l2;
        };

    const auto exp2 =
        [](float x) -> float
        {
            return pow_F(2.f, x);
        };

     // Build the luma channels: band-pass filters with gaussian windows of
     // std 2 EV, spaced by 2 EV
    const float centers[12] = {
        -18.0f, -16.0f, -14.0f, -12.0f, -10.0f, -8.0f, -6.0f,
        -4.0f, -2.0f, 0.0f, 2.0f, 4.0f
    };

    const auto conv = [&](int v) -> float
                      {
                          return exp2(float(v) / 100.f * 2.f);
                      };

    const float factors[12] = {
        conv(pp.bands[0]), // -18 EV
        conv(pp.bands[0]), // -16 EV
        conv(pp.bands[0]), // -14 EV
        conv(pp.bands[0]), // -12 EV
        conv(pp.bands[0]), // -10 EV
        conv(pp.bands[0]), //  -8 EV (Blacks)
        conv(pp.bands[1]), //  -6 EV (Shadows)
        conv(pp.bands[2]), //  -4 EV (Midtones)
        conv(pp.bands[3]), //  -2 EV (Highlights)
        conv(pp.bands[4]), //   0 EV (Whites)
        conv(pp.bands[4]), //   2 EV
        conv(pp.bands[4])  //   4 EV
    };

    const auto gauss =
        [](float b, float x) -> float
        {
            return xexpf((-SQR(x - b) / 4.0f));
        };

    // For every pixel luminance, the sum of the gaussian masks
    // evenly spaced by 2 EV with 2 EV std should be this constant
    float w_sum = 0.f;
    for (int i = 0; i < 12; ++i) {
        w_sum += gauss(centers[i], 0.f);
    }

    const auto luma =
        [&](float y) -> float
        {
            return max(log2(max(y, 0.f)), -18.0f);
        };

#ifdef __SSE2__
    vfloat vfactors[12];
    vfloat vcenters[12];
    
    for (int i = 0; i < 12; ++i) {
        vfactors[i] = F2V(factors[i]);
        vcenters[i] = F2V(centers[i]);
    }

    const auto vgauss =
        [](vfloat b, vfloat x) -> vfloat
        {
            static const vfloat fourv = F2V(4.f);
            return xexpf((-SQR(x - b) / fourv));
        };

    vfloat zerov = F2V(0.f);
    vfloat vw_sum = F2V(w_sum);

    const vfloat noisev = F2V(-18.f);
    const vfloat xlog2v = F2V(xlogf(2.f));
    
    const auto vluma =
        [&](vfloat y) -> vfloat
        {
            return vmaxf(xlogf(vmaxf(y, zerov))/xlog2v, noisev);
        };
#endif // __SSE2__
    

    array2D<float> mask(W, H);
    array2D<float> corr(W, H);
    array2D<float> L(W, H);

    if (radius > 0) {
        rtengine::guidedFilter(Y, Y, Y, radius, epsilon, multithread);
        radius = 0;
    }

#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        int x = 0;
#ifdef __SSE2__
        for (; x < W - 3; x += 4) {
            STVFU(corr[y][x], zerov);
            STVFU(L[y][x], vluma(LVFU(Y[y][x])));
        }
#endif
        for (; x < W; ++x) {
            corr[y][x] = 0.f;
            L[y][x] = luma(Y[y][x]);
        }
    }

    for (int c = 0; c < 12; ++c) {
#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < H; ++y) {
            int x = 0;
#ifdef __SSE2__
            for (; x < W - 3; x += 4) {
                const vfloat l = LVFU(L[y][x]);
                STVFU(mask[y][x], vgauss(vcenters[c], l));
            }
#endif
            for (; x < W; ++x) {
                const float l = L[y][x];
                mask[y][x] = gauss(centers[c], l);
            }
        }

        if (radius > 0) {
            rtengine::guidedFilter(Y, mask, mask, radius, epsilon, multithread);
        }

#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < H; ++y) {
            int x = 0;
#ifdef __SSE2__
            for (; x < W - 3; x += 4) {
                vfloat m = LVFU(mask[y][x]);
                STVFU(corr[y][x], LVFU(corr[y][x]) + m * vfactors[c]);
            }
#endif
            for (; x < W; ++x) {
                corr[y][x] += mask[y][x] * factors[c];
            }
        }
    }
#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        int x = 0;
#ifdef __SSE2__
        for (; x < W - 3; x += 4) {
            vfloat f = LVFU(corr[y][x]) / vw_sum;
            STVFU(R[y][x], LVFU(R[y][x]) * f);
            STVFU(G[y][x], LVFU(G[y][x]) * f);
            STVFU(B[y][x], LVFU(B[y][x]) * f);
        }
#endif
        for (; x < W; ++x) {
            float f = corr[y][x] / w_sum;
            R[y][x] *= f;
            G[y][x] *= f;
            B[y][x] *= f;
        }
    }
}
#endif

} // namespace


void ImProcFunctions::toneEqualizer(Imagefloat *rgb)
{
    if (!params->toneEqualizer.enabled) {
        return;
    }

    BENCHFUN

    rgb->normalizeFloatTo1();

    const int W = rgb->getWidth();
    const int H = rgb->getHeight();
    
    array2D<float> R(W, H, rgb->r.ptrs, ARRAY2D_BYREFERENCE);
    array2D<float> G(W, H, rgb->g.ptrs, ARRAY2D_BYREFERENCE);
    array2D<float> B(W, H, rgb->b.ptrs, ARRAY2D_BYREFERENCE);

    if (true) {
        tone_eq_gf(R, G, B, params->toneEqualizer, params->icm.workingProfile, scale, multiThread);
    } else {
    array2D<float> Yorig(W, H);
    array2D<float> Ytoned(W, H);

    tone_eq(R, G, B, Yorig, Ytoned, params->toneEqualizer, params->icm.workingProfile, scale, multiThread);

    const auto gaussian_coef =
        [](float x, float mu, float sigma) -> float
        {
            return (std::exp(-(x - mu) * (x - mu)/(sigma * sigma) / 2.0f) / (std::sqrt(2.0f * RT_PI) * sigma));
        };

    const auto normalize =
        [](KernelRow &kernel) -> void
        {
            float sum = 0.0f;
            for (size_t i = 0; i < kernel.size(); ++i) {
                sum += kernel[i];
            }
            for (size_t i = 0; i < kernel.size(); ++i) {
                kernel[i] /= sum;
            }
        };


    if (params->toneEqualizer.detail > 0) {
        const float sigma = max(params->toneEqualizer.detail, 0) / scale;
        KernelRow gauss;
        for (size_t i = 0; i < gauss.size(); ++i) {
            gauss[i] = gaussian_coef(i - 2.f, 0.f, sigma);
        }
        normalize(gauss);
    
        Kernel kernel;
        for (size_t m = 0; m < kernel.size(); ++m) {
            for (size_t n = 0; n < kernel.size(); ++n) {
                kernel[m][n] = gauss[m] * gauss[n];
            }
        }
        laplacian_filter(Yorig, Ytoned, Ytoned, kernel, multiThread);

#ifdef _OPENMP
#       pragma omp parallel for if (multiThread)
#endif
        for (int y = 0; y < H; ++y) {
            for (int x = 0; x < W; ++x) {
                float f = Ytoned[y][x];
                R[y][x] *= f;
                G[y][x] *= f;
                B[y][x] *= f;
            }
        }
    }
    }
    
    rgb->normalizeFloatTo65535();
}

} // namespace rtengine
