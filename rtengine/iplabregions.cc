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

bool ImProcFunctions::colorCorrection(Imagefloat *rgb, int offset_x, int offset_y, int full_width, int full_height)
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

    LabImage tmplab(rgb->getWidth(), rgb->getHeight());
    rgb2lab(*rgb, tmplab);
    LabImage *lab = &tmplab;

    if (editWhatever) {
        LabMasksEditID id = static_cast<LabMasksEditID>(int(eid) - EUID_LabMasks_H1);
        fillPipetteLabMasks(lab, editWhatever, id, multiThread);
    }
    
    int n = params->colorcorrection.regions.size();
    int show_mask_idx = params->colorcorrection.showMask;
    if (show_mask_idx >= n) {
        show_mask_idx = -1;
    }

    std::vector<array2D<float>> abmask(n);
    std::vector<array2D<float>> Lmask(n);

    if (!generateLabMasks(lab, params->colorcorrection.labmasks, offset_x, offset_y, full_width, full_height, scale, multiThread, show_mask_idx, &Lmask, &abmask)) {
        lab2rgb(*lab, *rgb);
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
    int channel[n];
    for (int i = 0; i < n; ++i) {
        auto &r = params->colorcorrection.regions[i];
        abca[i] = abcoord(r.a);
        abcb[i] = abcoord(r.b);
        rs[i] = 1.f + r.saturation / 100.f;
        slope[i] = r.slope;
        offset[i] = r.offset;
        power[i] = r.power;
        channel[i] = r.channel;
    }

    TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
    TMatrix iws = ICCStore::getInstance()->workingSpaceInverseMatrix(params->icm.workingProfile);

    const auto CDL =
        [=](float &l, float &a, float &b, float slope, float offset, float power, float saturation) -> void
        {
            if (slope != 1.f || offset != 0.f || power != 1.f || saturation != 1.f) {
                float rgb[3];
                float x, y, z;
                Color::Lab2XYZ(l, a, b, x, y, z);
                Color::xyz2rgb(x, y, z, rgb[0], rgb[1], rgb[2], iws);
                for (int i = 0; i < 3; ++i) {
                    float v = (rgb[i] / 65535.f) * slope + offset;
                    if (v > 0.f) {
                        v = pow_F(v, power);
                    }
                    rgb[i] = v * 65535.f;
                    //rgb[i] = (pow_F(max((rgb[i] / 65535.f) * slope + offset, 0.f), power)) * 65535.f;
                }
                if (saturation != 1.f) {
                    float Y = max(Color::rgbLuminance(rgb[0], rgb[1], rgb[2], ws), 0.f);
                    for (int i = 0; i < 3; ++i) {
                        rgb[i] = max(Y + saturation * (rgb[i] - Y), 0.f);
                    }
                }
                Color::rgbxyz(rgb[0], rgb[1], rgb[2], x, y, z, ws);
                Color::XYZ2Lab(x, y, z, l, a, b);
            }
        };

    const auto chan =
        [=](float prev_l, float prev_a, float prev_b, float &l, float &a, float &b, int channel) -> void
        {
            if (channel >= 0) {
                float prev_rgb[3];
                float rgb[3];
                float x, y, z;
                Color::Lab2XYZ(l, a, b, x, y, z);
                Color::xyz2rgb(x, y, z, rgb[0], rgb[1], rgb[2], iws);
                Color::Lab2XYZ(prev_l, prev_a, prev_b, x, y, z);
                Color::xyz2rgb(x, y, z, prev_rgb[0], prev_rgb[1], prev_rgb[2], iws);
                prev_rgb[channel] = rgb[channel];
                Color::rgbxyz(prev_rgb[0], prev_rgb[1], prev_rgb[2], x, y, z, ws);
                Color::XYZ2Lab(x, y, z, l, a, b);
            }
        };

#ifdef __SSE2__
    const auto CDL_v =
        [=](vfloat &l, vfloat &a, vfloat &b, float slope, float offset, float power, float saturation) -> void
        {
            if (slope != 1.f || offset != 0.f || power != 1.f || saturation != 1.f) {
                float ll[4];
                float aa[4];
                float bb[4];
                STVFU(ll[0], l);
                STVFU(aa[0], a);
                STVFU(bb[0], b);
                for (int i = 0; i < 4; ++i) {
                    CDL(ll[i], aa[i], bb[i], slope, offset, power, saturation);
                }
                l = LVFU(ll[0]);
                a = LVFU(aa[0]);
                b = LVFU(bb[0]);
            }
        };

    const auto chan_v =
        [=](vfloat prev_l, vfloat prev_a, vfloat prev_b, vfloat &l, vfloat &a, vfloat &b, int channel) -> void
        {
            if (channel >= 0) {
                float ll[4];
                float aa[4];
                float bb[4];
                STVFU(ll[0], l);
                STVFU(aa[0], a);
                STVFU(bb[0], b);
                float prev_ll[4];
                float prev_aa[4];
                float prev_bb[4];
                STVFU(prev_ll[0], prev_l);
                STVFU(prev_aa[0], prev_a);
                STVFU(prev_bb[0], prev_b);
                for (int i = 0; i < 4; ++i) {
                    chan(prev_ll[i], prev_aa[i], prev_bb[i], ll[i], aa[i], bb[i], channel);
                }
                l = LVFU(ll[0]);
                a = LVFU(aa[0]);
                b = LVFU(bb[0]);
            }
        };
#endif

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
// #ifdef __SSE2__
//         vfloat c42000v = F2V(42000.f);
//         vfloat cm42000v = F2V(-42000.f);
// #endif
#ifdef _OPENMP
        #pragma omp for
#endif
        for (int y = 0; y < lab->H; ++y) {
            int x = 0;
#ifdef __SSE2__
            for (; x < lab->W - 3; x += 4) {
                vfloat lv = LVFU(lab->L[y][x]);
                vfloat av = LVFU(lab->a[y][x]);
                vfloat bv = LVFU(lab->b[y][x]);

                for (int i = 0; i < n; ++i) {
                    vfloat blendv = LVFU(abmask[i][y][x]);
                    vfloat lblendv = LVFU(Lmask[i][y][x]);
                    // vmask b = vorm(vmaskf_gt(blendv, ZEROV), vmaskf_gt(lblendv, ZEROV));
                    // if (!_mm_movemask_ps((vfloat)b)) {
                    //     continue;
                    // }
                    
                    vfloat l_newv = lv;
                    vfloat lfv = vmaxf(lv, ZEROV);
                    vfloat a_newv = /*vclampf*/(av + lfv * F2V(abca[i]));//, cm42000v, c42000v);
                    vfloat b_newv = /*vclampf*/(bv + lfv * F2V(abcb[i]));//, cm42000v, c42000v);
                    CDL_v(l_newv, a_newv, b_newv, slope[i], offset[i], power[i], rs[i]);
                    //l_newv = vmaxf(l_newv, ZEROV);
                    chan_v(lv, av, bv, l_newv, a_newv, b_newv, channel[i]);
                    lv = vintpf(lblendv, l_newv, lv);
                    av = vintpf(blendv, a_newv, av);
                    bv = vintpf(blendv, b_newv, bv);
                }
                STVFU(lab->L[y][x], lv);
                STVFU(lab->a[y][x], av);
                STVFU(lab->b[y][x], bv);
            }
#endif
            for (; x < lab->W; ++x) {
                float l = lab->L[y][x];
                float a = lab->a[y][x];
                float b = lab->b[y][x];

                for (int i = 0; i < n; ++i) {
                    float blend = abmask[i][y][x];
                    float lblend = Lmask[i][y][x];
                    // if (!(blend > 0.f || lblend > 0.f)) {
                    //     continue;
                    // }
                    
                    float l_new = l;
                    float lf = max(l, 0.f);
                    float a_new = /*LIM*/(a + lf * abca[i]);//, -42000.f, 42000.f);
                    float b_new = /*LIM*/(b + lf * abcb[i]);//, -42000.f, 42000.f);
                    CDL(l_new, a_new, b_new, slope[i], offset[i], power[i], rs[i]);
                    //l_new = max(l_new, 0.f);
                    chan(l, a, b, l_new, a_new, b_new, channel[i]);
                    l = intp(lblend, l_new, l);
                    a = intp(blend, a_new, a);
                    b = intp(blend, b_new, b);
                }
                lab->L[y][x] = l;
                lab->a[y][x] = a;
                lab->b[y][x] = b;
            }
        }
    }

    lab2rgb(*lab, *rgb);

    return false;
}

} // namespace rtengine
