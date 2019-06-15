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

// taken from darktable (src/iop/ashift.c)
/*
  This file is part of darktable,
  copyright (c) 2016 Ulrich Pegelow.

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
// Inspiration to this module comes from the program ShiftN (http://www.shiftn.de) by
// Marcus Hebel.

// Thanks to Marcus for his support when implementing part of the ShiftN functionality
// to darktable.


#include "perspectivecorrection.h"
#include "rt_math.h"
#include <string.h>
#include <math.h>

extern "C" {

#define NMS_EPSILON 1e-3                    // break criterion for Nelder-Mead simplex
#define NMS_SCALE 1.0                       // scaling factor for Nelder-Mead simplex
#define NMS_ITERATIONS 400                  // number of iterations for Nelder-Mead simplex
#define NMS_CROP_EPSILON 100.0              // break criterion for Nelder-Mead simplex on crop fitting
#define NMS_CROP_SCALE 0.5                  // scaling factor for Nelder-Mead simplex on crop fitting
#define NMS_CROP_ITERATIONS 100             // number of iterations for Nelder-Mead simplex on crop fitting
#define NMS_ALPHA 1.0                       // reflection coefficient for Nelder-Mead simplex
#define NMS_BETA 0.5                        // contraction coefficient for Nelder-Mead simplex
#define NMS_GAMMA 2.0                       // expansion coefficient for Nelder-Mead simplex
#define DEFAULT_F_LENGTH 28.0               // focal length we assume if no exif data are available


#include "ashift_lsd.c"
#include "ashift_nmsimplex.c"

} // extern "C"

namespace rtengine {

namespace {

#define MAT3SWAP(a, b) { float (*tmp)[3] = (a); (a) = (b); (b) = tmp; }
//#define CLAMP(v, lo, hi) LIM((v), (lo), (hi))

// multiply 3x3 matrix with 3x1 vector
// dst needs to be different from v
inline void mat3mulv(float *dst, const float *const mat, const float *const v)
{
    for(int k = 0; k < 3; k++)
    {
        float x = 0.0f;
        for(int i = 0; i < 3; i++) x += mat[3 * k + i] * v[i];
        dst[k] = x;
    }
}


// multiply two 3x3 matrices
// dst needs to be different from m1 and m2
inline void mat3mul(float *dst, const float *const m1, const float *const m2)
{
    for(int k = 0; k < 3; k++) {
        for(int i = 0; i < 3; i++) {
            float x = 0.0f;
            for(int j = 0; j < 3; j++) x += m1[3 * k + j] * m2[3 * j + i];
            dst[3 * k + i] = x;
        }
    }
}


inline int mat3inv(float *const dst, const float *const src)
{
    std::array<std::array<float, 3>, 3> tmpsrc;
    std::array<std::array<float, 3>, 3> tmpdst;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            tmpsrc[i][j] = src[3 * i + j];
        }
    }
    if (invertMatrix(tmpsrc, tmpdst)) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                dst[3 * i + j] = tmpdst[i][j];
            }
        }
        return 0;
    } else {
        return 1;
    }
}


void homography(float *homograph, const float angle, const float shift_v, const float shift_h,
                const float shear, const float f_length_kb, const float orthocorr, const float aspect,
                const int width, const int height, bool reverse)
{
    // calculate homograph that combines all translations, rotations
    // and warping into one single matrix operation.
    // this is heavily leaning on ShiftN where the homographic matrix expects
    // input in (y : x : 1) format. in the darktable world we want to keep the
    // (x : y : 1) convention. therefore we need to flip coordinates first and
    // make sure that output is in correct format after corrections are applied.

    const float u = width;
    const float v = height;

    const float phi = RT_PI_F * angle / 180.0f;
    const float cosi = cos(phi);
    const float sini = sin(phi);
    const float ascale = sqrt(aspect);

    // most of this comes from ShiftN
    const float f_global = f_length_kb;
    const float horifac = 1.0f - orthocorr / 100.0f;
    const float exppa_v = exp(shift_v);
    const float fdb_v = f_global / (14.4f + (v / u - 1) * 7.2f);
    const float rad_v = fdb_v * (exppa_v - 1.0f) / (exppa_v + 1.0f);
    const float alpha_v = CLAMP(atanf(rad_v), -1.5f, 1.5f);
    const float rt_v = sin(0.5f * alpha_v);
    const float r_v = fmax(0.1f, 2.0f * (horifac - 1.0f) * rt_v * rt_v + 1.0f);

    const float vertifac = 1.0f - orthocorr / 100.0f;
    const float exppa_h = exp(shift_h);
    const float fdb_h = f_global / (14.4f + (u / v - 1) * 7.2f);
    const float rad_h = fdb_h * (exppa_h - 1.0f) / (exppa_h + 1.0f);
    const float alpha_h = CLAMP(atanf(rad_h), -1.5f, 1.5f);
    const float rt_h = sin(0.5f * alpha_h);
    const float r_h = fmax(0.1f, 2.0f * (vertifac - 1.0f) * rt_h * rt_h + 1.0f);


    // three intermediate buffers for matrix calculation ...
    float m1[3][3], m2[3][3], m3[3][3];

    // ... and some pointers to handle them more intuitively
    float (*mwork)[3] = m1;
    float (*minput)[3] = m2;
    float (*moutput)[3] = m3;

    // Step 1: flip x and y coordinates (see above)
    memset(minput, 0, 9 * sizeof(float));
    minput[0][1] = 1.0f;
    minput[1][0] = 1.0f;
    minput[2][2] = 1.0f;


    // Step 2: rotation of image around its center
    memset(mwork, 0, 9 * sizeof(float));
    mwork[0][0] = cosi;
    mwork[0][1] = -sini;
    mwork[1][0] = sini;
    mwork[1][1] = cosi;
    mwork[0][2] = -0.5f * v * cosi + 0.5f * u * sini + 0.5f * v;
    mwork[1][2] = -0.5f * v * sini - 0.5f * u * cosi + 0.5f * u;
    mwork[2][2] = 1.0f;

    // multiply mwork * minput -> moutput
    mat3mul((float *)moutput, (float *)mwork, (float *)minput);


    // Step 3: apply shearing
    memset(mwork, 0, 9 * sizeof(float));
    mwork[0][0] = 1.0f;
    mwork[0][1] = shear;
    mwork[1][1] = 1.0f;
    mwork[1][0] = shear;
    mwork[2][2] = 1.0f;

    // moutput (of last calculation) -> minput
    MAT3SWAP(minput, moutput);
    // multiply mwork * minput -> moutput
    mat3mul((float *)moutput, (float *)mwork, (float *)minput);


    // Step 4: apply vertical lens shift effect
    memset(mwork, 0, 9 * sizeof(float));
    mwork[0][0] = exppa_v;
    mwork[1][0] = 0.5f * ((exppa_v - 1.0f) * u) / v;
    mwork[1][1] = 2.0f * exppa_v / (exppa_v + 1.0f);
    mwork[1][2] = -0.5f * ((exppa_v - 1.0f) * u) / (exppa_v + 1.0f);
    mwork[2][0] = (exppa_v - 1.0f) / v;
    mwork[2][2] = 1.0f;

    // moutput (of last calculation) -> minput
    MAT3SWAP(minput, moutput);
    // multiply mwork * minput -> moutput
    mat3mul((float *)moutput, (float *)mwork, (float *)minput);


    // Step 5: horizontal compression
    memset(mwork, 0, 9 * sizeof(float));
    mwork[0][0] = 1.0f;
    mwork[1][1] = r_v;
    mwork[1][2] = 0.5f * u * (1.0f - r_v);
    mwork[2][2] = 1.0f;

    // moutput (of last calculation) -> minput
    MAT3SWAP(minput, moutput);
    // multiply mwork * minput -> moutput
    mat3mul((float *)moutput, (float *)mwork, (float *)minput);


    // Step 6: flip x and y back again
    memset(mwork, 0, 9 * sizeof(float));
    mwork[0][1] = 1.0f;
    mwork[1][0] = 1.0f;
    mwork[2][2] = 1.0f;

    // moutput (of last calculation) -> minput
    MAT3SWAP(minput, moutput);
    // multiply mwork * minput -> moutput
    mat3mul((float *)moutput, (float *)mwork, (float *)minput);


    // from here output vectors would be in (x : y : 1) format

    // Step 7: now we can apply horizontal lens shift with the same matrix format as above
    memset(mwork, 0, 9 * sizeof(float));
    mwork[0][0] = exppa_h;
    mwork[1][0] = 0.5f * ((exppa_h - 1.0f) * v) / u;
    mwork[1][1] = 2.0f * exppa_h / (exppa_h + 1.0f);
    mwork[1][2] = -0.5f * ((exppa_h - 1.0f) * v) / (exppa_h + 1.0f);
    mwork[2][0] = (exppa_h - 1.0f) / u;
    mwork[2][2] = 1.0f;

    // moutput (of last calculation) -> minput
    MAT3SWAP(minput, moutput);
    // multiply mwork * minput -> moutput
    mat3mul((float *)moutput, (float *)mwork, (float *)minput);


    // Step 8: vertical compression
    memset(mwork, 0, 9 * sizeof(float));
    mwork[0][0] = 1.0f;
    mwork[1][1] = r_h;
    mwork[1][2] = 0.5f * v * (1.0f - r_h);
    mwork[2][2] = 1.0f;

    // moutput (of last calculation) -> minput
    MAT3SWAP(minput, moutput);
    // multiply mwork * minput -> moutput
    mat3mul((float *)moutput, (float *)mwork, (float *)minput);


    // Step 9: apply aspect ratio scaling
    memset(mwork, 0, 9 * sizeof(float));
    mwork[0][0] = 1.0f * ascale;
    mwork[1][1] = 1.0f / ascale;
    mwork[2][2] = 1.0f;

    // moutput (of last calculation) -> minput
    MAT3SWAP(minput, moutput);
    // multiply mwork * minput -> moutput
    mat3mul((float *)moutput, (float *)mwork, (float *)minput);


    // Step 10: find x/y offsets and apply according correction so that
    // no negative coordinates occur in output vector
    float umin = FLT_MAX, vmin = FLT_MAX;
    // visit all four corners
    for(int y = 0; y < height; y += height - 1)
        for(int x = 0; x < width; x += width - 1)
        {
            float pi[3], po[3];
            pi[0] = x;
            pi[1] = y;
            pi[2] = 1.0f;
            // moutput expects input in (x:y:1) format and gives output as (x:y:1)
            mat3mulv(po, (float *)moutput, pi);
            umin = fmin(umin, po[0] / po[2]);
            vmin = fmin(vmin, po[1] / po[2]);
        }

    memset(mwork, 0, 9 * sizeof(float));
    mwork[0][0] = 1.0f;
    mwork[1][1] = 1.0f;
    mwork[2][2] = 1.0f;
    mwork[0][2] = -umin;
    mwork[1][2] = -vmin;

    // moutput (of last calculation) -> minput
    MAT3SWAP(minput, moutput);
    // multiply mwork * minput -> moutput
    mat3mul((float *)moutput, (float *)mwork, (float *)minput);


    // on request we either keep the final matrix for forward conversions
    // or produce an inverted matrix for backward conversions
    if (!reverse) {
        // we have what we need -> copy it to the right place
        memcpy(homograph, moutput, 9 * sizeof(float));
    } else {
        // generate inverted homograph (mat3inv function defined in colorspaces.c)
        if(mat3inv((float *)homograph, (float *)moutput))
        {
            // in case of error we set to unity matrix
            memset(mwork, 0, 9 * sizeof(float));
            mwork[0][0] = 1.0f;
            mwork[1][1] = 1.0f;
            mwork[2][2] = 1.0f;
            memcpy(homograph, mwork, 9 * sizeof(float));
        }
    }
}

} // namespace


PerspectiveCorrection::PerspectiveCorrection():
    ok_(false),
    scale_(1.0),
    offx_(0.0),
    offy_(0.0),
    scalein_(1.0),
    offxin_(0.0),
    offyin_(0.0)
{
}


void PerspectiveCorrection::init(int width, int height, const procparams::PerspectiveParams &params, bool fill)
{
    const float f_length = 28.f;
    homography((float *)ihomograph_, params.angle, params.vertical / 100.0, params.horizontal / 100.0, params.shear / 100.0, f_length, 0.f, 1.f, width, height, true);

    ok_ = true;
    calc_scale(width, height, params, fill);
}


inline void PerspectiveCorrection::correct(double &x, double &y, double scale, double offx, double offy)
{
    if (ok_) {       
        float pin[3], pout[3];
        pout[0] = x;
        pout[1] = y;
        pout[0] *= scale_;
        pout[1] *= scale_;
        pout[0] += offx_;
        pout[1] += offy_;
        pout[2] = 1.f;
        mat3mulv(pin, (float *)ihomograph_, pout);
        pin[0] /= pin[2];
        pin[1] /= pin[2];
        x = pin[0] * scale + offx;
        y = pin[1] * scale + offy;
    }
}


void PerspectiveCorrection::operator()(double &x, double &y)
{
    correct(x, y, scalein_, offxin_, offyin_);
}


namespace {

std::vector<Coord2D> get_samples(int w, int h)
{
    constexpr int nsteps = 32;

    int x1 = 0, y1 = 0;
    int x2 = x1 + w - 1;
    int y2 = y1 + h - 1;

    // Build all edge points and half-way points
    std::vector<Coord2D> corners(8);
    corners[0].set(x1, y1);
    corners[1].set(x1, y2);
    corners[2].set(x2, y2);
    corners[3].set(x2, y1);
    corners[4].set((x1 + x2) / 2, y1);
    corners[5].set((x1 + x2) / 2, y2);
    corners[6].set(x1, (y1 + y2) / 2);
    corners[7].set(x2, (y1 + y2) / 2);

    // Add several steps inbetween
    int xstep = (x2 - x1) / nsteps;

    if (xstep < 1) {
        xstep = 1;
    }

    for (int i = x1 + xstep; i <= x2 - xstep; i += xstep) {
        corners.push_back(Coord2D(i, y1));
        corners.push_back(Coord2D(i, y2));
    }

    int ystep = (y2 - y1) / nsteps;

    if (ystep < 1) {
        ystep = 1;
    }

    for (int i = y1 + ystep; i <= y2 - ystep; i += ystep) {
        corners.push_back(Coord2D(x1, i));
        corners.push_back(Coord2D(x2, i));
    }

    return corners;
}

} // namespace

void PerspectiveCorrection::calc_scale(int w, int h, const procparams::PerspectiveParams &params, bool fill)
{
    double min_x = RT_INFINITY, max_x = -RT_INFINITY;
    double min_y = RT_INFINITY, max_y = -RT_INFINITY;

    auto corners = get_samples(w, h);

    float homo[3][3];
    constexpr float f_length = 28.f;
    homography((float *)homo, params.angle, params.vertical / 100.0, params.horizontal / 100.0, params.shear / 100.0, f_length, 0.f, 1.f, w, h, false);
    
    for (auto &c : corners) {
        float pin[3] = { float(c.x), float(c.y), 1.f };
        float pout[3];
        mat3mulv(pout, (float *)homo, pin);
        double x = pout[0] / pout[2];
        double y = pout[1] / pout[2];
        min_x = min(min_x, x);
        max_x = max(max_x, x);
        min_y = min(min_y, y);
        max_y = max(max_y, y);
    }

    double cw = max_x - min_x;
    double ch = max_y - min_y;

    scale_ = max(cw / double(w), ch / double(h));
    offx_ = (cw - w * scale_) * 0.5;
    offy_ = (ch - h * scale_) * 0.5;

    if (fill) {
        double scaleU = 2, scaleL = 0.001;  // upper and lower border, iterate inbetween
        do {
            double scale = (scaleU + scaleL) * 0.5;
            bool clipped = test_scale(w, h, scale);

            if (clipped) {
                scaleU = scale;
            } else {
                scaleL = scale;
            }
        } while (scaleU - scaleL > 0.001);

        scalein_ = scaleL;
        offxin_ = (w - scalein_ * w) * 0.5;
        offyin_ = (h - scalein_ * h) * 0.5;
        return;
    }
}


bool PerspectiveCorrection::test_scale(int w, int h, double scale)
{
    auto samples = get_samples(w, h);
    double offx = (w - scale * w) * 0.5;
    double offy = (h - scale * h) * 0.5;
    for (auto &s : samples) {
        int x = s.x, y = s.y;
        correct(s.x, s.y, scale, offx, offy);
        if (s.x < 0 || s.x > w - 1 || s.y < 0 || s.y > h - 1) {
            return true;
        }
    }
    return false;
}

} // namespace rtengine
