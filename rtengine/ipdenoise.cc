/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
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
#include "imagesource.h"
#include "mytime.h"

namespace rtengine {

extern const Settings *settings;

namespace {

void adjust_params(procparams::DenoiseParams &dnparams, double scale)
{
    if (scale <= 1.0) {
        return;
    }
    
    double scale_factor = 1.0 / scale;
    double noise_factor_c = std::pow(scale_factor, 0.46);//scale_factor * 1.8;
    double noise_factor_l = std::pow(scale_factor, 0.62);//scale_factor * 1.5; //std::pow(scale_factor, scale_factor);
    //dnparams.luminance *= noise_factor_l;
    dnparams.luminance *= noise_factor_l;
    //dnparams.luminanceDetail += dnparams.luminanceDetail * noise_factor_l;
    dnparams.luminanceDetail *= (1.0 + std::pow(1.0 - scale_factor, 2.2));//scale_factor * 1.2);
    // if (dnparams.chrominanceMethod == procparams::DenoiseParams::ChrominanceMethod::MANUAL) {
        dnparams.chrominance *= noise_factor_c;
        dnparams.chrominanceRedGreen *= noise_factor_c;
        dnparams.chrominanceBlueYellow *= noise_factor_c;
    // }
    // if (scale > 2) {
    //     dnparams.aggressive = false;
    // }

    if (dnparams.smoothingEnabled) {
        int j = int(dnparams.medianType) - int(1.0 / scale_factor);
        if (j < 0 && dnparams.smoothingMethod == procparams::DenoiseParams::SmoothingMethod::MEDIAN) {
            dnparams.smoothingEnabled = false;
        } else {
            dnparams.medianType = static_cast<procparams::DenoiseParams::MedianType>(j);
        }
    }

    // guided params don't need to be adjusted here, they are already scaled in ImProcFunctions::guidedSmoothing
}

} // namespace


void ImProcFunctions::denoiseComputeParams(ImageSource *imgsrc, const ColorTemp &currWB, DenoiseInfoStore &store, procparams::DenoiseParams &dnparams)
{
    float autoNR = settings->nrauto;
    float autoNRmax = settings->nrautomax;

    int tilesize;
    int overlap;

    if (settings->leveldnti == 0) {
        tilesize = 1024;
        overlap = 128;
    }

    if (settings->leveldnti == 1) {
        tilesize = 768;
        overlap = 96;
    }
    
    int numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip;
    int kall = 2;

    int widIm, heiIm;
    int tr = getCoarseBitMask(params->coarse);
    imgsrc->getFullSize(widIm, heiIm, tr);
    
    Tile_calc(tilesize, overlap, kall, widIm, heiIm, numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip);
    kall = 0;

    float min_b[9];
    float min_r[9];
    float lumL[9];
    float chromC[9];
    float ry[9];
    float sk[9];
    float pcsk[9];
    std::vector<int> centerTile_X(numtiles_W);
    std::vector<int> centerTile_Y(numtiles_H);

    for (int cX = 0; cX < numtiles_W; cX++) {
        centerTile_X[cX] = tileWskip / 2 + tileWskip * cX;
    }

    for (int cY = 0; cY < numtiles_H; cY++) {
        centerTile_Y[cY] = tileHskip / 2 + tileHskip * cY;
    }

    assert(settings->leveldnautsimpl == 0);
    
    if (!store.valid && dnparams.chrominanceMethod == procparams::DenoiseParams::ChrominanceMethod::AUTOMATIC) {
        MyTime t1aue, t2aue;
        t1aue.set();

        int crW = 100; // settings->leveldnv == 0
        int crH = 100; // settings->leveldnv == 0

        if (settings->leveldnv == 1) {
            crW = 250;
            crH = 250;
        }

        //  if(settings->leveldnv ==2) {crW=int(tileWskip/2);crH=int((tileWskip/2));}//adapted to scale of preview
        if (settings->leveldnv == 2) {
            crW = int (tileWskip / 2);
            crH = int (tileHskip / 2);
        }

        if (settings->leveldnv == 3) {
            crW = tileWskip - 10;
            crH = tileHskip - 10;
        }

        float lowdenoise = 1.f;
        int levaut = settings->leveldnaut;

        if (levaut == 1) { //Standard
            lowdenoise = 0.7f;
        }

        LUTf gamcurve(65536, 0);
        float gam, gamthresh, gamslope;
        RGB_denoise_infoGamCurve(dnparams, imgsrc->isRAW(), gamcurve, gam, gamthresh, gamslope);
        int Nb[9];

#ifdef _OPENMP
#pragma omp parallel if (multiThread)
#endif
        {
            Imagefloat *origCropPart = new Imagefloat(crW, crH); //allocate memory
            Imagefloat *provicalc = new Imagefloat((crW + 1) / 2, (crH + 1) / 2);  //for denoise curves

            int  coordW[3];//coordinate of part of image to measure noise
            int  coordH[3];
            int begW = 50;
            int begH = 50;
            coordW[0] = begW;
            coordW[1] = widIm / 2 - crW / 2;
            coordW[2] = widIm - crW - begW;
            coordH[0] = begH;
            coordH[1] = heiIm / 2 - crH / 2;
            coordH[2] = heiIm - crH - begH;

#ifdef _OPENMP
#           pragma omp for schedule(dynamic) collapse(2) nowait
#endif

            for (int wcr = 0; wcr <= 2; wcr++) {
                for (int hcr = 0; hcr <= 2; hcr++) {
                    PreviewProps ppP(coordW[wcr], coordH[hcr], crW, crH, 1);
                    imgsrc->getImage(currWB, tr, origCropPart, ppP, params->exposure, params->raw);

                    // we only need image reduced to 1/4 here
                    for (int ii = 0; ii < crH; ii += 2) {
                        for (int jj = 0; jj < crW; jj += 2) {
                            provicalc->r(ii >> 1, jj >> 1) = origCropPart->r(ii, jj);
                            provicalc->g(ii >> 1, jj >> 1) = origCropPart->g(ii, jj);
                            provicalc->b(ii >> 1, jj >> 1) = origCropPart->b(ii, jj);
                        }
                    }

                    imgsrc->convertColorSpace(provicalc, params->icm, currWB);  //for denoise luminance curve

                    float pondcorrec = 1.0f;
                    float chaut = 0.f, redaut = 0.f, blueaut = 0.f, maxredaut = 0.f, maxblueaut = 0.f, minredaut = 0.f, minblueaut = 0.f, chromina = 0.f, sigma = 0.f, lumema = 0.f, sigma_L = 0.f, redyel = 0.f, skinc = 0.f, nsknc = 0.f;
                    int nb = 0;
                    RGB_denoise_info(origCropPart, provicalc, imgsrc->isRAW(), gamcurve, gam, gamthresh, gamslope, dnparams, imgsrc->getDirPyrDenoiseExpComp(), chaut, nb, redaut, blueaut, maxredaut, maxblueaut, minredaut, minblueaut, chromina, sigma, lumema, sigma_L, redyel, skinc, nsknc);

                    //printf("DCROP skip=%d cha=%f red=%f bl=%f redM=%f bluM=%f chrom=%f sigm=%f lum=%f\n",skip, chaut,redaut,blueaut, maxredaut, maxblueaut, chromina, sigma, lumema);
                    Nb[hcr * 3 + wcr] = nb;
                    store.ch_M[hcr * 3 + wcr] = pondcorrec * chaut;
                    store.max_r[hcr * 3 + wcr] = pondcorrec * maxredaut;
                    store.max_b[hcr * 3 + wcr] = pondcorrec * maxblueaut;
                    min_r[hcr * 3 + wcr] = pondcorrec * minredaut;
                    min_b[hcr * 3 + wcr] = pondcorrec * minblueaut;
                    lumL[hcr * 3 + wcr] = lumema;
                    chromC[hcr * 3 + wcr] = chromina;
                    ry[hcr * 3 + wcr] = redyel;
                    sk[hcr * 3 + wcr] = skinc;
                    pcsk[hcr * 3 + wcr] = nsknc;

                }
            }

            delete provicalc;
            delete origCropPart;
        }
        float chM = 0.f;
        float MaxR = 0.f;
        float MaxB = 0.f;
        float MinR = 100000000000.f;
        float MinB = 100000000000.f;
        float maxr = 0.f;
        float maxb = 0.f;
        float Max_R[9] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
        float Max_B[9] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
        float Min_R[9];
        float Min_B[9];
        float MaxRMoy = 0.f;
        float MaxBMoy = 0.f;
        float MinRMoy = 0.f;
        float MinBMoy = 0.f;

        float multip = 1.f;

        if (!imgsrc->isRAW()) {
            multip = 2.f;    //take into account gamma for TIF / JPG approximate value...not good for gamma=1
        }

        float adjustr = 1.f;

        // if (params->icm.workingProfile == "ProPhoto")   {
        //     adjustr = 1.f;   //
        // } else if (params->icm.workingProfile == "Adobe RGB")  {
        //     adjustr = 1.f / 1.3f;
        // } else if (params->icm.workingProfile == "sRGB")       {
        //     adjustr = 1.f / 1.3f;
        // } else if (params->icm.workingProfile == "WideGamut")  {
        //     adjustr = 1.f / 1.1f;
        // } else if (params->icm.workingProfile == "Beta RGB")   {
        //     adjustr = 1.f / 1.2f;
        // } else if (params->icm.workingProfile == "BestRGB")    {
        //     adjustr = 1.f / 1.2f;
        // } else if (params->icm.workingProfile == "BruceRGB")   {
        //     adjustr = 1.f / 1.2f;
        // }

        float delta[9];
        int mode = 1;
        int lissage = settings->leveldnliss;

        for (int k = 0; k < 9; k++) {
            float maxmax = max(store.max_r[k], store.max_b[k]);
            calcautodn_info(store.ch_M[k], delta[k], Nb[k], levaut, maxmax, lumL[k], chromC[k], mode, lissage, ry[k], sk[k], pcsk[k]);
            //  printf("ch_M=%f delta=%f\n",ch_M[k], delta[k]);
        }

        for (int k = 0; k < 9; k++) {
            if (store.max_r[k] > store.max_b[k]) {
                Max_R[k] = (delta[k]) / ((autoNRmax * multip * adjustr * lowdenoise) / 2.f);
                Min_B[k] = - (store.ch_M[k] - min_b[k]) / (autoNRmax * multip * adjustr * lowdenoise);
                Max_B[k] = 0.f;
                Min_R[k] = 0.f;
            } else {
                Max_B[k] = (delta[k]) / ((autoNRmax * multip * adjustr * lowdenoise) / 2.f);
                Min_R[k] = - (store.ch_M[k] - min_r[k])   / (autoNRmax * multip * adjustr * lowdenoise);
                Min_B[k] = 0.f;
                Max_R[k] = 0.f;
            }
        }

        for (int k = 0; k < 9; k++) {
            //  printf("ch_M= %f Max_R=%f Max_B=%f min_r=%f min_b=%f\n",ch_M[k],Max_R[k], Max_B[k],Min_R[k], Min_B[k]);
            chM += store.ch_M[k];
            MaxBMoy += Max_B[k];
            MaxRMoy += Max_R[k];
            MinRMoy += Min_R[k];
            MinBMoy += Min_B[k];

            if (Max_R[k] > MaxR) {
                MaxR = Max_R[k];
            }

            if (Max_B[k] > MaxB) {
                MaxB = Max_B[k];
            }

            if (Min_R[k] < MinR) {
                MinR = Min_R[k];
            }

            if (Min_B[k] < MinB) {
                MinB = Min_B[k];
            }
        }

        chM /= 9;
        MaxBMoy /= 9;
        MaxRMoy /= 9;
        MinBMoy /= 9;
        MinRMoy /= 9;

        if (MaxR > MaxB) {
            maxr = MaxRMoy + (MaxR - MaxRMoy) * 0.66f; //#std Dev
            //maxb=MinB;
            maxb = MinBMoy + (MinB - MinBMoy) * 0.66f;
        } else {
            maxb = MaxBMoy + (MaxB - MaxBMoy) * 0.66f;
            maxr = MinRMoy + (MinR - MinRMoy) * 0.66f;
        }

//                  printf("DCROP skip=%d cha=%f red=%f bl=%f \n",skip, chM,maxr,maxb);
        dnparams.chrominance = chM / (autoNR * multip * adjustr);
        dnparams.chrominanceRedGreen = maxr;
        dnparams.chrominanceBlueYellow = maxb;

        dnparams.chrominance *= dnparams.chrominanceAutoFactor;
        dnparams.chrominanceRedGreen *= dnparams.chrominanceAutoFactor;
        dnparams.chrominanceBlueYellow *= dnparams.chrominanceAutoFactor;
        
        store.valid = true;

        if (settings->verbose) {
            t2aue.set();
            printf("Info denoise auto performed in %d usec:\n", t2aue.etime(t1aue));
        }
        //end evaluate noise
    }
}


void ImProcFunctions::denoise(ImageSource *imgsrc, const ColorTemp &currWB, Imagefloat *img, const DenoiseInfoStore &store, const procparams::DenoiseParams &dnparams)
{
    procparams::DenoiseParams denoiseParams = dnparams;
    NoiseCurve noiseLCurve;
    NoiseCurve noiseCCurve;

    noiseCCurve.Set({
        FCT_MinMaxCPoints,
        0.05,
        0.50,
        0.35,
        0.35,
        0.35,
        0.05,
        0.35,
        0.35
        });

    Imagefloat *calclum = nullptr;
    {
        const int fw = img->getWidth();
        const int fh = img->getHeight();
        // we only need image reduced to 1/4 here
        calclum = new Imagefloat((fw + 1) / 2, (fh + 1) / 2); //for luminance denoise curve
#ifdef _OPENMP
#       pragma omp parallel for if (multiThread)
#endif

        for (int ii = 0; ii < fh; ii += 2) {
            for (int jj = 0; jj < fw; jj += 2) {
                calclum->r(ii >> 1, jj >> 1) = img->r(ii, jj);
                calclum->g(ii >> 1, jj >> 1) = img->g(ii, jj);
                calclum->b(ii >> 1, jj >> 1) = img->b(ii, jj);
            }
        }
        imgsrc->convertColorSpace(calclum, params->icm, currWB);
    }
    
    float nresi, highresi;
    DenoiseInfoStore &dnstore = const_cast<DenoiseInfoStore &>(store);

    adjust_params(denoiseParams, scale);

//     array2D<float> Y(img->getWidth(), img->getHeight());
//     TMatrix ws = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
// #ifdef _OPENMP
// #    pragma omp parallel for if (multiThread)
// #endif
//     for (int y = 0; y < Y.height(); ++y) {
//         for (int x = 0; x < Y.width(); ++x) {
//             Y[y][x] = Color::rgbLuminance(img->r(y, x), img->g(y, x), img->b(y, x), ws);
//         }
//     }
    
    RGB_denoise(0, img, img, calclum, dnstore.ch_M, dnstore.max_r, dnstore.max_b, imgsrc->isRAW(), denoiseParams, imgsrc->getDirPyrDenoiseExpComp(), noiseLCurve, noiseCCurve, nresi, highresi);

//     LUTf lcurve(65536, LUT_CLIP_ABOVE|LUT_CLIP_BELOW);
//     {
//         FlatCurve curve({
//         FCT_MinMaxCPoints,
//             0, 1,
//             0.35, 0,
//             1, 0.75,
//             0.399383, 0.35
//             }, CURVES_MIN_POLY_POINTS / scale);
//         const bool raw = imgsrc->isRAW();
//         const float gamma = denoiseParams.gamma;
//         int m = 65535.0 / scale;
//         for (int i = 0; i <= m; ++i) {
//             double x = double(i) / double(m);
//             if (raw) {
//                 x = Color::gammanf(x, gamma);
//             }
//             lcurve[i] = curve.getVal(x);
//         }
//     }

// #ifdef _OPENMP
// #   pragma omp parallel for if (multiThread)
// #endif
//     for (int y = 0; y < Y.height(); ++y) {
//         for (int x = 0; x < Y.width(); ++x) {
//             float iY = Y[y][x];
//             float oY = Color::rgbLuminance(img->r(y, x), img->g(y, x), img->b(y, x), ws);
//             if (oY > 1e-5f) {
//                 float blend = lcurve[iY];//noiseLCurve[gammalut[CLIP(iY)]];
//                 iY = intp(blend, oY, iY);
//                 float f = iY / oY;
//                 img->r(y, x) *= f;
//                 img->g(y, x) *= f;
//                 img->b(y, x) *= f;
//             }
//         }
//     }
//     Y.free();
    
    guidedSmoothing(img);
}


} // namespace rtengine
