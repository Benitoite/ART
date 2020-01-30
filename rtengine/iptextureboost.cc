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
#include "labmasks.h"
#include "array2D.h"
#include "EdgePreservingDecomposition.h"

namespace rtengine {

namespace {

void EPD(LabImage *lab, const rtengine::procparams::TextureBoostParams::Region &pp, double scale, bool multithread)
{
    unsigned int Iterates = 5;

    float stren = pp.strength; 
    float edgest = pp.edgeStopping; 
    float sca = pp.scale;
    constexpr float gamm = 1.5f;
    constexpr float rew = 1;
    constexpr float noise_floor = 1e-8f * 32768.f;
    
    //Pointers to whole data and size of it.
    float *L = lab->L[0];
    float *a = lab->a[0];
    float *b = lab->b[0];
    size_t N = lab->W * lab->H;
    EdgePreservingDecomposition epd (lab->W, lab->H);

    //Due to the taking of logarithms, L must be nonnegative. Further, scale to 0 to 1 using nominal range of L, 0 to 15 bit.
    float minL = FLT_MAX;
    float maxL = 0.f;
#ifdef _OPENMP
    #pragma omp parallel if (multithread)
#endif
    {
        float lminL = FLT_MAX;
        float lmaxL = 0.f;
#ifdef _OPENMP
        #pragma omp for
#endif
        for (size_t i = 0; i < N; i++) {
            L[i] = std::max(L[i], noise_floor);
            
            if (L[i] < lminL) {
                lminL = L[i];
            }

            if (L[i] > lmaxL) {
                lmaxL = L[i];
            }
        }

#ifdef _OPENMP
        #pragma omp critical
#endif
        {
            if (lminL < minL) {
                minL = lminL;
            }

            if (lmaxL > maxL) {
                maxL = lmaxL;
            }
        }
    }

    if (minL > 0.0f) {
        minL = 0.0f;    //Disable the shift if there are no negative numbers. I wish there were just no negative numbers to begin with.
    }

    if (maxL == 0.f) { // avoid division by zero
        maxL = 1.f;
    }

#ifdef _OPENMP
    #pragma omp parallel for if (multithread)
#endif
    for (size_t i = 0; i < N; ++i)
        //{L[i] = (L[i] - minL)/32767.0f;
    {
        L[i] = (L[i] - minL) / maxL;
        L[i] *= gamm;
    }

    //Some interpretations.
    float Compression = expf (-stren);      //This modification turns numbers symmetric around 0 into exponents.
    float DetailBoost = stren;

    if (stren < 0.0f) {
        DetailBoost = 0.0f;    //Go with effect of exponent only if uncompressing.
    }

    //Auto select number of iterates. Note that p->EdgeStopping = 0 makes a Gaussian blur.
    if (Iterates == 0) {
        Iterates = (unsigned int) (edgest * 15.0f);
    }

    /* Debuggery. Saves L for toying with outside of RT.
    char nm[64];
    sprintf(nm, "%ux%ufloat.bin", lab->W, lab->H);
    FILE *f = fopen(nm, "wb");
    fwrite(L, N, sizeof(float), f);
    fclose(f);*/

    epd.CompressDynamicRange (L, sca / scale, edgest, Compression, DetailBoost, Iterates, rew);

    //Restore past range, also desaturate a bit per Mantiuk's Color correction for tone mapping.
    float s = (1.0f + 38.7889f) * powf (Compression, 1.5856f) / (1.0f + 38.7889f * powf (Compression, 1.5856f));
#ifdef _OPENMP
    #pragma omp parallel for            // removed schedule(dynamic,10)
#endif
    for (size_t ii = 0; ii < N; ++ii) {
        a[ii] *= s;
        b[ii] *= s;
        L[ii] = L[ii] * maxL * (1.f / gamm) + minL;
    }
}


} // namespace


bool ImProcFunctions::textureBoost(Imagefloat *rgb)
{
    PlanarWhateverData<float> *editWhatever = nullptr;
    EditUniqueID eid = pipetteBuffer ? pipetteBuffer->getEditID() : EUID_None;

    if ((eid == EUID_LabMasks_H4 || eid == EUID_LabMasks_C4 || eid == EUID_LabMasks_L4) && pipetteBuffer->getDataProvider()->getCurrSubscriber()->getPipetteBufferType() == BT_SINGLEPLANE_FLOAT) {
        editWhatever = pipetteBuffer->getSinglePlaneBuffer();
    }
    
    if (eid == EUID_LabMasks_DE4) {
        if (getDeltaEColor(rgb, deltaE.x, deltaE.y, offset_x, offset_y, full_width, full_height, scale, deltaE.L, deltaE.C, deltaE.H)) {
            deltaE.ok = true;
        }
    }
    
    if (params->textureBoost.enabled) {
        if (editWhatever) {
            LabMasksEditID id = static_cast<LabMasksEditID>(int(eid) - EUID_LabMasks_H4);
            fillPipetteLabMasks(rgb, editWhatever, id, multiThread);
        }
        
        int n = params->textureBoost.regions.size();
        int show_mask_idx = params->textureBoost.showMask;
        if (show_mask_idx >= n || (cur_pipeline != Pipeline::PREVIEW && cur_pipeline != Pipeline::OUTPUT)) {
            show_mask_idx = -1;
        }
        std::vector<array2D<float>> mask(n);
        if (!generateLabMasks(rgb, params->textureBoost.labmasks, offset_x, offset_y, full_width, full_height, scale, multiThread, show_mask_idx, &mask, nullptr)) {
            return true; // show mask is active, nothing more to do
        }

        LabImage lab(rgb->getWidth(), rgb->getHeight());
        rgb2lab(*rgb, lab);

        array2D<float> L(lab.W, lab.H, lab.L, 0);
        array2D<float> a(lab.W, lab.H, lab.a, 0);
        array2D<float> b(lab.W, lab.H, lab.b, 0);

        for (int i = 0; i < n; ++i) {
            if (!params->textureBoost.labmasks[i].enabled) {
                continue;
            }
            
            auto &r = params->textureBoost.regions[i];
            EPD(&lab, r, scale, multiThread);
            const auto &blend = mask[i];

#ifdef _OPENMP
#           pragma omp parallel for if (multiThread)
#endif
            for (int y = 0; y < lab.H; ++y) {
                for (int x = 0; x < lab.W; ++x) {
                    float l = lab.L[y][x];
                    float aa = lab.a[y][x];
                    float bb = lab.b[y][x];
                    lab.L[y][x] = intp(blend[y][x], lab.L[y][x], L[y][x]);
                    lab.a[y][x] = intp(blend[y][x], lab.a[y][x], a[y][x]);
                    lab.b[y][x] = intp(blend[y][x], lab.b[y][x], b[y][x]);
                    L[y][x] = l;
                    a[y][x] = aa;
                    b[y][x] = bb;
                }
            }
        }

        lab2rgb(lab, *rgb);
    } else if (editWhatever) {
        editWhatever->fill(0.f);
    }

    return false;
}

} // namespace rtengine
