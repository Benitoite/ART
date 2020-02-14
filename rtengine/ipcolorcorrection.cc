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
#include "sleef.h"
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

    if (eid == EUID_LabMasks_DE1) {
        if (getDeltaEColor(rgb, deltaE.x, deltaE.y, offset_x, offset_y, full_width, full_height, scale, deltaE.L, deltaE.C, deltaE.H)) {
            deltaE.ok = true;
        }
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
    if (show_mask_idx >= n || (cur_pipeline != Pipeline::PREVIEW && cur_pipeline != Pipeline::OUTPUT)) {
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
            return SGN(x) * xlog2lin(std::abs(x), 4.f);
        };

    TMatrix dws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
    float ws[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            ws[i][j] = dws[i][j];
        }
    }

    float abca[n];
    float abcb[n];
    float rs[n];
    float rslope[n][3];
    float roffset[n][3];
    float rpower[n][3];
    float rpivot[n][3];
    bool rgbmode[n];
    bool enabled[n];
    for (int i = 0; i < n; ++i) {
        auto &r = params->colorcorrection.regions[i];
        rgbmode[i] = r.mode != ColorCorrectionParams::Mode::YUV;
        if (rgbmode[i]) {
            abca[i] = 0.f;
            abcb[i] = 0.f;
            rs[i] = 1.f;
        } else {
            abca[i] = abcoord(r.a);
            abcb[i] = abcoord(r.b);
            rs[i] = 1.f + r.saturation / 100.f;
        }
        enabled[i] = false;
        if (r.mode == ColorCorrectionParams::Mode::HSL) {
            float R, G, B;
            float u, v;
            for (int c = 0; c < 3; ++c) {
                float hue = (float(r.hue[c]) / 180.f) * rtengine::RT_PI;
                float sat = pow_F(float(r.sat[c]) / 100.f, 3.f);
                float f = (r.factor[c] / 100.f) + 1.f;
                Color::hsl2yuv(hue, sat, u, v);
                Color::yuv2rgb(0.5f, u, v, R, G, B, ws);
                R *= 2.f;
                G *= 2.f;
                B *= 2.f;
                switch (c) {
                case 0: // SLOPE
                    rslope[i][0] = R * f;
                    rslope[i][1] = G * f;
                    rslope[i][2] = B * f;
                    break;
                case 1: // OFFSET
                    roffset[i][0] = R + f - 2.f;
                    roffset[i][1] = G + f - 2.f;
                    roffset[i][2] = B + f - 2.f;
                    break;
                default: // POWER
                    rpower[i][0] = (2.f - R) * (2.f - f);
                    rpower[i][1] = (2.f - G) * (2.f - f);
                    rpower[i][2] = (2.f - B) * (2.f - f);
                    break;
                }
                rpivot[i][c] = 1.f;
            }
            for (int c = 0; c < 3; ++c) {
                if (rslope[i][c] != 1.f || roffset[i][c] != 0.f || rpower[i][c] != 1.f) {
                    enabled[i] = true;
                }
            }
        } else {
            for (int c = 0; c < 3; ++c) {
                int j = rgbmode[i] ? c : 0;
                rslope[i][c] = r.slope[j];
                roffset[i][c] = r.offset[j];
                rpower[i][c] = r.power[j];
                rpivot[i][c] = r.pivot[j];
                if (rslope[i][c] != 1.f || roffset[i][c] != 0.f || rpower[i][c] != 1.f) {
                    enabled[i] = true;
                }
            }
        }
    }

    const auto CDL = 
        [&](int region, float &Y, float &u, float &v) -> void
        {
            if (enabled[region]) {
                const float *slope = rslope[region];
                const float *offset = roffset[region];
                const float *power = rpower[region];
                const float *pivot = rpivot[region];
                
                if (rgbmode[region]) {
                    float rgb[3];
                    Color::yuv2rgb(Y, u, v, rgb[0], rgb[1], rgb[2], ws);
                    for (int i = 0; i < 3; ++i) {
                        float v = (rgb[i] / 65535.f) * slope[i] + offset[i]/2.f;
                        if (v > 0.f) {
                            if (pivot[i] != 1.f) {
                                v = pow_F(v / pivot[i], power[i]) * pivot[i];
                            } else {
                                v = pow_F(v, power[i]);
                            }
                        } else {
                            v = 0.f;
                        }
                        rgb[i] = v * 65535.f;
                    }
                    Color::rgb2yuv(rgb[0], rgb[1], rgb[2], Y, u, v, ws);
                } else {
                    float v = (Y / 65535.f) * slope[0] + offset[0]/2.f;
                    if (v > 0.f) {
                        if (pivot[0] != 1.f) {
                            v = pow_F(v / pivot[0], power[0]) * pivot[0];
                        } else {
                            v = pow_F(v, power[0]);
                        }
                        v *= 65535.f;
                    } else {
                        v = 0.f;
                    }
                    if (Y > 0.f) {
                        float f = v / Y;
                        Y = v;
                        u *= f;
                        v *= f;
                    } else {
                        Y = v;
                    }
                }
            }
            const float saturation = rs[region];
            if (saturation != 1.f) {
                u *= saturation;
                v *= saturation;
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
        [&](int region, vfloat &Y, vfloat &u, vfloat &v) -> void
        {
            if (enabled[region]) {
                const float *slope = rslope[region];
                const float *offset = roffset[region];
                const float *power = rpower[region];
                const float *pivot = rpivot[region];
                
                if (rgbmode[region]) {
                    vfloat rgb[3];
                    vfloat v65535 = F2V(65535.f);
                    Color::yuv2rgb(Y, u, v, rgb[0], rgb[1], rgb[2], vws);
                    for (int i = 0; i < 3; ++i) {
                        vfloat vslope = F2V(slope[i]);
                        vfloat voffset = F2V(offset[i] / 2.f);
                        vfloat vpower = F2V(power[i]);
                        vfloat vpivot = F2V(pivot[i]);
                
                        vfloat v = (rgb[i] / v65535) * vslope + voffset;
                        if (pivot[i] != 1.f) {
                            v = vself(vmaskf_gt(v, ZEROV), pow_F(v / vpivot, vpower) * vpivot, ZEROV);
                        } else {
                            v = vself(vmaskf_gt(v, ZEROV), pow_F(v, vpower), ZEROV);
                        }
                        rgb[i] = v * v65535;
                    }
                    Color::rgb2yuv(rgb[0], rgb[1], rgb[2], Y, u, v, vws);
                } else {
                    vfloat vslope = F2V(slope[0]);
                    vfloat voffset = F2V(offset[0] / 2.f);
                    vfloat vpower = F2V(power[0]);
                    vfloat vpivot = F2V(pivot[0]);
                    vfloat v65535 = F2V(65535.f);
                    vfloat YY = (Y / v65535) * vslope + voffset;
                    if (pivot[0] != 1.f) {
                        YY = vself(vmaskf_gt(YY, ZEROV), pow_F(YY / vpivot, vpower) * vpivot, ZEROV);
                    } else {
                        YY = vself(vmaskf_gt(YY, ZEROV), pow_F(YY, vpower), ZEROV);
                    }
                    YY *= v65535;
                    vfloat f = vself(vmaskf_gt(Y, ZEROV), YY / Y, F2V(1.f));
                    Y = YY;
                    u *= f;
                    v *= f;
                }
            }
            const float saturation = rs[region];
            if (saturation != 1.f) {
                vfloat vsaturation = F2V(saturation);
                u *= vsaturation;
                v *= vsaturation;
            }
        };

    vfloat vabcb[n];
    vfloat vabca[n];
    for (int i = 0; i < n; ++i) {
        vabcb[i] = F2V(abcb[i]);
        vabca[i] = F2V(abca[i]);
    }

    vfloat zerov = F2V(0.0);
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
                    if (!params->colorcorrection.labmasks[i].enabled) {
                        continue;
                    }
                    
                    vfloat blendv = LVFU(abmask[i][y][x]);
                    vfloat lblendv = LVFU(Lmask[i][y][x]);

                    vmask some_blend = vmaskf_gt(blendv, zerov);
                    vmask some_lblend = vmaskf_gt(lblendv, zerov);
                    if (_mm_movemask_ps((vfloat)some_blend) || _mm_movemask_ps((vfloat)some_lblend)) {
                        vfloat Y_newv = Yv;
                        vfloat u_newv = uv;
                        vfloat v_newv = vv;

                        CDL_v(i, Y_newv, u_newv, v_newv);
                    
                        vfloat fv = vmaxf(Y_newv, ZEROV);
                        u_newv += fv * vabcb[i];
                        v_newv += fv * vabca[i];
                    
                        Yv = vintpf(lblendv, Y_newv, Yv);
                        uv = vintpf(blendv, u_newv, uv);
                        vv = vintpf(blendv, v_newv, vv);
                    }
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
                    if (!params->colorcorrection.labmasks[i].enabled) {
                        continue;
                    }

                    float blend = abmask[i][y][x];
                    float lblend = Lmask[i][y][x];

                    if (blend > 0.f || lblend > 0.f) {
                        float Y_new = Y;
                        float u_new = u;
                        float v_new = v;

                        CDL(i, Y_new, u_new, v_new);
                    
                        float f = max(Y_new, 0.f);
                        u_new += f * abcb[i];
                        v_new += f * abca[i];
                    
                        Y = intp(lblend, Y_new, Y);
                        u = intp(blend, u_new, u);
                        v = intp(blend, v_new, v);
                    }
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
