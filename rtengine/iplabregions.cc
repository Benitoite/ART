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

#ifdef _OPENMP
#include <omp.h>
#endif

#include "improcfun.h"
#include "guidedfilter.h"
//#define BENCHMARK
#include "StopWatch.h"
#include "sleef.c"
#include "coord.h"
#include "gauss.h"
#include "labmasks.h"

namespace rtengine {

bool ImProcFunctions::colorCorrection(Imagefloat *rgb)
{
    PlanarWhateverData<float> *editWhatever = nullptr;
    EditUniqueID eid = pipetteBuffer ? pipetteBuffer->getEditID() : EUID_None;

    if ((eid == EUID_LabMasks_H1 || eid == EUID_LabMasks_C1 || eid == EUID_LabMasks_L1) && pipetteBuffer->getDataProvider()->getCurrSubscriber()->getPipetteBufferType() == BT_SINGLEPLANE_FLOAT) {
        editWhatever = pipetteBuffer->getSinglePlaneBuffer();
    }
    
    if (!params->colorcorrection.enabled) {
        if (editWhatever) {
            editWhatever->fill(0.f);
        }
        return false;
    }

    if (editWhatever) {
        LabMasksEditID id = static_cast<LabMasksEditID>(int(eid) - EUID_LabMasks_H1);
        fillPipetteLabMasks(rgb, editWhatever, id, multiThread);
    }
    
    int n = params->colorcorrection.regions.size();
    int show_mask_idx = params->colorcorrection.showMask;
    if (show_mask_idx >= n) {
        show_mask_idx = -1;
    }

    std::vector<array2D<float>> abmask(n);
    std::vector<array2D<float>> Lmask(n);

    if (!generateLabMasks(rgb, params->colorcorrection.labmasks, offset_x, offset_y, full_width, full_height, scale, multiThread, show_mask_idx, &Lmask, &abmask)) {
        return true; // show mask is active, nothing more to do
    }
    
    const auto abcoord =
        [](float x) -> float
        {
            return /*12000.f **/ SGN(x) * xlog2lin(std::abs(x), 4.f);
        };

    float abca[n];
    float abcb[n];
    float rs[n];
    float slope[n];
    float offset[n];
    float power[n];
    float pivot[n];
    int channel[n];
    for (int i = 0; i < n; ++i) {
        auto &r = params->colorcorrection.regions[i];
        abca[i] = abcoord(r.a);
        abcb[i] = abcoord(r.b);
        rs[i] = 1.f + r.saturation / 100.f;
        slope[i] = r.slope;
        offset[i] = r.offset;
        power[i] = r.power;
        pivot[i] = r.pivot;
        channel[i] = r.channel;
    }

    TMatrix dws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
    float ws[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            ws[i][j] = dws[i][j];
        }
    }

    // const auto rgb2yuv =
    //     [=](float r, float g, float b, float &Y, float &u, float &v) -> void
    //     {
    //         Color::rgb2yuv(r, g, b, Y, u, v, ws);
    //     };

    // const auto yuv2rgb =
    //     [=](float Y, float u, float v, float &r, float &g, float &b) -> void
    //     {
    //         Color::yuv2rgb(Y, u, v, r, g, b, ws);
    //     };

    const auto CDL = 
        [=, &ws](float &Y, float &u, float &v, float slope, float offset, float power, float pivot, float saturation) -> void
        {
            float rgb[3];
            if (slope != 1.f || offset != 0.f || power != 1.f) {
                Color::yuv2rgb(Y, u, v, rgb[0], rgb[1], rgb[2], ws);
                for (int i = 0; i < 3; ++i) {
                    float v = (rgb[i] / 65535.f) * slope + offset;
                    if (v > 0.f) {
                        if (pivot != 1.f) {
                            v = pow_F(v / pivot, power) * pivot;
                        } else {
                            v = pow_F(v, power);
                        }
                    }
                    rgb[i] = v * 65535.f;
                }
                Color::rgb2yuv(rgb[0], rgb[1], rgb[2], Y, u, v, ws);
            }
            if (saturation != 1.f) {
                u *= saturation;
                v *= saturation;
            }
        };

    const auto chan =
        [=, &ws](float prev_Y, float prev_u, float prev_v, float &Y, float &u, float &v, int channel) -> void
        {
            if (channel >= 0) {
                float prev_rgb[3];
                Color::yuv2rgb(prev_Y, prev_u, prev_v, prev_rgb[0], prev_rgb[1], prev_rgb[2], ws);
                float rgb[3];
                Color::yuv2rgb(Y, u, v, rgb[0], rgb[1], rgb[2], ws);
                prev_rgb[channel] = rgb[channel];
                Color::rgb2yuv(prev_rgb[0], prev_rgb[1], prev_rgb[2], Y, u, v, ws);
            }
        };
    
#ifdef __SSE2__
    vfloat vws[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            vws[i][j] = F2V(float(dws[i][j]));
        }
    }

    const auto CDL_v =
        [=](vfloat &Y, vfloat &u, vfloat &v, float slope, float offset, float power, float pivot, float saturation) -> void
        {
            if (slope != 1.f || offset != 0.f || power != 1.f) {// || saturation != 1.f) {
                vfloat rgb[3];
                vfloat vslope = F2V(slope);
                vfloat voffset = F2V(offset);
                vfloat vpower = F2V(power);
                vfloat vpivot = F2V(pivot);
                vfloat v65535 = F2V(65535.f);
                Color::yuv2rgb(Y, u, v, rgb[0], rgb[1], rgb[2], vws);
                for (int i = 0; i < 3; ++i) {
                    vfloat v = (rgb[i] / v65535) * vslope + voffset;
                    if (pivot != 1.f) {
                        v = vself(vmaskf_gt(v, ZEROV), pow_F(v / vpivot, vpower) * vpivot, v);
                    } else {
                        v = vself(vmaskf_gt(v, ZEROV), pow_F(v, vpower), v);
                    }
                    rgb[i] = v * v65535;
                }
                Color::rgb2yuv(rgb[0], rgb[1], rgb[2], Y, u, v, vws);
            }
            if (saturation != 1.f) {
                vfloat vsaturation = F2V(saturation);
                u *= vsaturation;
                v *= vsaturation;
            }
        };

    const auto chan_v =
        [=](vfloat prev_Y, vfloat prev_u, vfloat prev_v, vfloat &Y, vfloat &u, vfloat &v, int channel) -> void
        {
            if (channel >= 0) {
                vfloat prev_rgb[3];
                Color::yuv2rgb(prev_Y, prev_u, prev_v, prev_rgb[0], prev_rgb[1], prev_rgb[2], vws);
                vfloat rgb[3];
                Color::yuv2rgb(Y, u, v, rgb[0], rgb[1], rgb[2], vws);
                prev_rgb[channel] = rgb[channel];
                Color::rgb2yuv(prev_rgb[0], prev_rgb[1], prev_rgb[2], Y, u, v, vws);
            }
        };
#endif

    const int H = rgb->getHeight();
    const int W = rgb->getWidth();

    rgb->setMode(Imagefloat::Mode::YUV, multiThread);

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
#ifdef _OPENMP
        #pragma omp for
#endif
        for (int y = 0; y < H; ++y) {
            int x = 0;
#ifdef __SSE2__
            for (; x < W - 3; x += 4) {
                vfloat Yv = LVF(rgb->g(y, x));
                vfloat uv = LVF(rgb->b(y, x));
                vfloat vv = LVF(rgb->r(y, x));

                for (int i = 0; i < n; ++i) {
                    vfloat blendv = LVFU(abmask[i][y][x]);
                    vfloat lblendv = LVFU(Lmask[i][y][x]);
                    vfloat fv = vmaxf(Yv, ZEROV);
                    vfloat Y_newv = Yv;
                    vfloat u_newv = uv + fv * F2V(abcb[i]);
                    vfloat v_newv = vv + fv * F2V(abca[i]);

                    CDL_v(Y_newv, u_newv, v_newv, slope[i], offset[i], power[i], pivot[i], rs[i]);
                    chan_v(Yv, uv, vv, Y_newv, u_newv, v_newv, channel[i]);
                    Yv = vintpf(lblendv, Y_newv, Yv);
                    uv = vintpf(blendv, u_newv, uv);
                    vv = vintpf(blendv, v_newv, vv);
                }
                STVF(rgb->g(y, x), Yv);
                STVF(rgb->b(y, x), uv);
                STVF(rgb->r(y, x), vv);
            }
#endif
            for (; x < W; ++x) {
                float Y = rgb->g(y, x);
                float u = rgb->b(y, x);
                float v = rgb->r(y, x);

                for (int i = 0; i < n; ++i) {
                    float blend = abmask[i][y][x];
                    float lblend = Lmask[i][y][x];

                    float f = max(Y, 0.f);
                    float Y_new = Y;
                    float u_new = u + f * abcb[i];
                    float v_new = v + f * abca[i];

                    CDL(Y_new, u_new, v_new, slope[i], offset[i], power[i], pivot[i], rs[i]);
                    chan(Y, u, v, Y_new, u_new, v_new, channel[i]);
                    Y = intp(lblend, Y_new, Y);
                    u = intp(blend, u_new, u);
                    v = intp(blend, v_new, v);
                }

                rgb->g(y, x) = Y;
                rgb->b(y, x) = u;
                rgb->r(y, x) = v;
            }
        }
    }

    return false;
}

} // namespace rtengine
