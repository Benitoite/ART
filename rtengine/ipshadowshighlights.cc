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

#include "improcfun.h"
#include "gauss.h"
#include "sleef.c"
#include "opthelper.h"
#include "guidedfilter.h"

namespace rtengine {

// namespace {

// class LabImageAdapter {
// public:
//     LabImageAdapter(LabImage *img, const Glib::ustring &working_profile):
//         img_(img)
//     {
//         ws_ = ICCStore::getInstance()->workingSpaceMatrix(working_profile);
//         iws_ = ICCStore::getInstance()->workingSpaceInverseMatrix(working_profile);
//     }

//     inline int get_width() { return img_->W; }
//     inline int get_height() { return img_->H; }

//     inline float get_L(int y, int x)
//     {
//         return img_->L[y][x];
//     }

//     inline void get_Lab(int y, int x, float &l, float &a, float &b)
//     {
//         l = img_->L[y][x];
//         a = img_->a[y][x];
//         b = img_->b[y][x];
//     }

//     inline void set_Lab(int y, int x, float l, float a, float b)
//     {
//         img_->L[y][x] = l;
//         img_->a[y][x] = a;
//         img_->b[y][x] = b;
//     }

//     inline void get_RGB(int y, int x, float &r, float &g, float &b)
//     {
//         Color::lab2rgb(img_->L[y][x], img_->a[y][x], img_->b[y][x], r, g, b, iws_);
//     }

//     inline void set_RGB(int y, int x, float r, float g, float b)
//     {
//         Color::rgb2lab(r, g, b, img_->L[y][x], img_->a[y][x], img_->b[y][x], ws_);
//     }
    
// private:
//     LabImage *img_;
//     TMatrix ws_;
//     TMatrix iws_;
// };


// class RGBImageAdapter {
// public:
//     RGBImageAdapter(Imagefloat *img, const Glib::ustring &working_profile):
//         img_(img)
//     {
//         ws_ = ICCStore::getInstance()->workingSpaceMatrix(working_profile);
//         iws_ = ICCStore::getInstance()->workingSpaceInverseMatrix(working_profile);
//     }

//     inline int get_width() { return img_->getWidth(); }
//     inline int get_height() { return img_->getHeight(); }

//     inline float get_L(int y, int x)
//     {
//         float l, a, b;
//         Color::rgb2lab(img_->r(y, x), img_->g(y, x), img_->b(y, x), l, a, b, ws_);
//         return l;
//     }

//     inline void get_Lab(int y, int x, float &l, float &a, float &b)
//     {
//         Color::rgb2lab(img_->r(y, x), img_->g(y, x), img_->b(y, x), l, a, b, ws_);
//     }

//     inline void set_Lab(int y, int x, float l, float a, float b)
//     {
//         Color::lab2rgb(l, a, b, img_->r(y, x), img_->g(y, x), img_->b(y, x), iws_);
//     }

//     inline void get_RGB(int y, int x, float &r, float &g, float &b)
//     {
//         r = img_->r(y, x);
//         g = img_->g(y, x);
//         b = img_->b(y, x);
//     }

//     inline void set_RGB(int y, int x, float r, float g, float b)
//     {
//         img_->r(y, x) = r;
//         img_->g(y, x) = g;
//         img_->b(y, x) = b;
//     }

// private:
//     Imagefloat *img_;
//     TMatrix ws_;
//     TMatrix iws_;
// };



// template <class Img>
// class ShadowsHighlights {
// public:
//     ShadowsHighlights(const ProcParams *params, double scale, bool multithread, Img *img):
//         params_(params),
//         scale_(scale),
//         multithread_(multithread),
//         img_(img)
//     {
//         ws_ = ICCStore::getInstance()->workingSpaceMatrix(params_->icm.workingProfile);
//         iws_ = ICCStore::getInstance()->workingSpaceInverseMatrix(params_->icm.workingProfile);
//     }

//     void operator()()
//     {
//         const int width = img_->get_width();
//         const int height = img_->get_height();
//         const bool lab_mode = params_->sh.lab;

//         array2D<float> mask(width, height);
//         array2D<float> L(width, height);
//         const float radius = float(params_->sh.radius) * 10 / scale_;
//         LUTf f(lab_mode ? 32768 : 65536);

//         const auto apply =
//             [&](int amount, int tonalwidth, bool hl) -> void
//             {
//                 const float thresh = tonalwidth * 327.68f;
//                 const float scale = hl ? (thresh > 0.f ? 0.9f / thresh : 1.f) : thresh * 0.9f;

// #ifdef _OPENMP
// #               pragma omp parallel for if (multithread_)
// #endif
//                 for (int y = 0; y < height; ++y) {
//                     for (int x = 0; x < width; ++x) {
//                         float l = img_->get_L(y, x);
//                         float l1 = l / 32768.f;
//                         if (hl) {
//                             mask[y][x] = (l > thresh) ? 1.f : pow4(l * scale);
//                             L[y][x] = 1.f - l1;
//                         } else {
//                             mask[y][x] = l <= thresh ? 1.f : pow4(scale / l);
//                             L[y][x] = l1;
//                         }
//                     }
//                 }

//                 rtengine::guidedFilter(L, mask, mask, radius, 0.075, multithread_, 4);

//                 const float base = std::pow(4.f, float(std::abs(amount))/100.f);
//                 const float gamma = (hl == (amount >= 0)) ? base : 1.f / base;

//                 const float contrast = std::pow(2.f, float(std::abs(amount))/100.f);
//                 DiagonalCurve sh_contrast({
//                         DCT_NURBS,
//                             0, 0,
//                             0.125, std::pow(0.125 / 0.25, contrast) * 0.25, 
//                             0.25, 0.25,
//                             0.375, std::pow(0.375 / 0.25, contrast) * 0.25,
//                             1, 1
//                             });

//                 if (lab_mode) {
// #ifdef _OPENMP
// #                   pragma omp parallel for if (multithread_)
// #endif
//                     for (int l = 0; l < 32768; ++l) {
//                         auto base = pow_F(l / 32768.f, gamma);
//                         if (!hl) {
//                             base = sh_contrast.getVal(base);
//                         }
//                         f[l] = base * 32768.f;
//                     }
//                 } else {
// #ifdef _OPENMP
// #                   pragma omp parallel for if (multithread_)
// #endif
//                     for (int c = 0; c < 65536; ++c) {
//                         float l, a, b;
//                         float R = c, G = c, B = c;
//                         Color::rgb2lab(R, G, B, l, a, b, ws_);
//                         auto base = pow_F(l / 32768.f, gamma);
//                         if (!hl) {
//                             base = sh_contrast.getVal(base);
//                         }
//                         l = base * 32768.f;
//                         Color::lab2rgb(l, a, b, R, G, B, iws_);
//                         f[c] = G;
//                     }
//                 }

// #ifdef _OPENMP
// #               pragma omp parallel for schedule(dynamic,16) if (multithread_)
// #endif
//                 for (int y = 0; y < height; ++y) {
//                     for (int x = 0; x < width; ++x) {
//                         float blend = LIM01(mask[y][x]);
//                         float orig = 1.f - blend;
//                         if (lab_mode) {
//                             float l, a, b;
//                             img_->get_Lab(y, x, l, a, b);
//                             if (l >= 0 && l < 32768.f) {
//                                 float ll = intp(blend, f[l], l);
//                                 if (!hl && l > 1.f) {
//                                     // when pushing shadows, scale also the chromaticity
//                                     float s = max(ll / l * 0.5f, 1.f) * blend;
//                                     a = a * s + a * orig;
//                                     b = b * s + b * orig;
//                                 }
//                                 img_->set_Lab(y, x, ll, a, b);
//                             }
//                         } else {
//                             float rgb[3];
//                             img_->get_RGB(y, x, rgb[0], rgb[1], rgb[2]);
//                             for (int i = 0; i < 3; ++i) {
//                                 float c = rgb[i];
//                                 if (!OOG(c)) {
//                                     rgb[i] = intp(blend, f[c], c);
//                                 }
//                             }
//                             img_->set_RGB(y, x, rgb[0], rgb[1], rgb[2]);
//                         }
//                     }
//                 }
//             };

//         if (params_->sh.highlights) {
//             apply(params_->sh.highlights * 0.7, params_->sh.htonalwidth, true);
//         }

//         if (params_->sh.shadows) {
//             apply(params_->sh.shadows * 0.6, params_->sh.stonalwidth, false);
//         }
//     }

// private:
//     const ProcParams *params_;
//     double scale_;
//     bool multithread_;
//     Img *img_;
//     TMatrix ws_;
//     TMatrix iws_;
// };

// } // namespace


// void ImProcFunctions::shadowsHighlights(LabImage *lab)
// {
//     if (!params->sh.enabled || (!params->sh.highlights && !params->sh.shadows)){
//         return;
//     }

//     LabImageAdapter img(lab, params->icm.workingProfile);
//     ShadowsHighlights<LabImageAdapter> sh(params, scale, multiThread, &img);
//     sh();
// }


// void ImProcFunctions::shadowsHighlights(Imagefloat *rgb)
// {
//     if (!params->sh.enabled || (!params->sh.highlights && !params->sh.shadows)){
//         return;
//     }

//     RGBImageAdapter img(rgb, params->icm.workingProfile);
//     ShadowsHighlights<RGBImageAdapter> sh(params, scale, multiThread, &img);
//     sh();
// }


namespace {

void sh(array2D<float> &R, array2D<float> &G, array2D<float> &B, const SHParams &pp, const Glib::ustring &workingProfile, double scale, bool multithread)
{
    // adapted from the tone equalizer of darktable
/*
    This file is part of darktable,
    copyright (c) 2018 Aurélien Pierre.

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
    
    const int W = R.width();
    const int H = R.height();
//    array2D<float> Y(W, H);
    array2D<float> RR(W, H, R, 0);
    array2D<float> GG(W, H, G, 0);
    array2D<float> BB(W, H, B, 0);

    const int r = max(int(30 / scale), 1);
    const float epsilon = 0.001f;

    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(workingProfile);

// #ifdef _OPENMP
// #   pragma omp parallel for if (multithread)
// #endif
//     for (int y = 0; y < H; ++y) {
//         for (int x = 0; x < W; ++x) {
//             Y[y][x] = Color::rgbLuminance(R[y][x], G[y][x], B[y][x], ws) / 65535.f;
//         }
//     }
    
//    rtengine::guidedFilter(Y, Y, Y, r, epsilon, multithread);
    rtengine::guidedFilter(RR, RR, RR, r, epsilon, multithread);
    rtengine::guidedFilter(GG, GG, GG, r, epsilon, multithread);
    rtengine::guidedFilter(BB, BB, BB, r, epsilon, multithread);

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

    const auto GAUSS = [](float b, float x) -> float
                       {
                           return xexpf((-SQR(x - b) / 4.0f));
                       };

     // Build the luma channels: band-pass filters with gaussian windows of
     // std 2 EV, spaced by 2 EV
    const float centers[12] = {
        -18.0f, -16.0f, -14.0f, -12.0f, -10.0f, -8.0f, -6.0f,
        -4.0f,  -2.0f,   0.0f,  2.0f,   4.0f
    };

    const auto conv = [&](int v) -> float
                      {
                          return exp2(float(v) / 100.f * 2.f);
                      };

    const float factors[12] = {
        conv(pp.levels[0]), // -18 EV
        conv(pp.levels[0]), // -16 EV
        conv(pp.levels[0]), // -14 EV
        conv(pp.levels[0]), // -12 EV
        conv(pp.levels[0]), // -10 EV
        conv(pp.levels[0]), //  -8 EV
        conv(pp.levels[1]), //  -6 EV
        conv(pp.levels[2]), //  -4 EV
        conv(pp.levels[3]), //  -2 EV
        conv(pp.levels[4]), //   0 EV
        conv(pp.levels[4]), //   2 EV
        conv(pp.levels[4])  //   4 EV
    };

    // For every pixel luminance, the sum of the gaussian masks
    // evenly spaced by 2 EV with 2 EV std should be this constant
    float w_sum = 0.f;
    for (int i = 0; i < 12; ++i) {
        w_sum += GAUSS(centers[i], 0.f);
    }

    const auto process_pixel =
        [&](float y) -> float
        {
            // get the luminance of the pixel - just average channels
            const float luma = max(log2(max(y, 0.f)), -18.0f);

            // build the correction as the sum of the contribution of each luminance channel to current pixel
            float correction = 0.0f;
            for (int c = 0; c < 12; ++c) {
                correction += GAUSS(centers[c], luma) * factors[c];
            }
            correction /= w_sum;

            return correction;
        };
        
#ifdef _OPENMP
#   pragma omp parallel for if (multithread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            float Y = Color::rgbLuminance(RR[y][x], GG[y][x], BB[y][x], ws);// / 65535.f;
            float corr = process_pixel(Y); //Y[y][x]);
            float dR = R[y][x] - RR[y][x];
            float dG = G[y][x] - GG[y][x];
            float dB = B[y][x] - BB[y][x];
            // R[y][x] *= corr;
            // G[y][x] *= corr;
            // B[y][x] *= corr;
            R[y][x] = RR[y][x] * corr + dR * (corr + 0.25f);
            G[y][x] = GG[y][x] * corr + dG * (corr + 0.25f);
            B[y][x] = BB[y][x] * corr + dB * (corr + 0.25f);
        }
    }
}

} // namespace


void ImProcFunctions::shadowsHighlights(LabImage *lab)
{
    if (!params->sh.enabled) {
        return;
    }

    Imagefloat rgb(lab->W, lab->H);
    lab2rgb(*lab, rgb, params->icm.workingProfile);
    shadowsHighlights(&rgb);
    rgb2lab(rgb, *lab, params->icm.workingProfile);
}


void ImProcFunctions::shadowsHighlights(Imagefloat *rgb)
{
    if (!params->sh.enabled) {
        return;
    }
    
    rgb->normalizeFloatTo1();

    array2D<float> R(rgb->getWidth(), rgb->getHeight(), rgb->r.ptrs, ARRAY2D_BYREFERENCE);
    array2D<float> G(rgb->getWidth(), rgb->getHeight(), rgb->g.ptrs, ARRAY2D_BYREFERENCE);
    array2D<float> B(rgb->getWidth(), rgb->getHeight(), rgb->b.ptrs, ARRAY2D_BYREFERENCE);

    sh(R, G, B, params->sh, params->icm.workingProfile, scale, multiThread);
    
    rgb->normalizeFloatTo65535();
}

} // namespace rtengine
