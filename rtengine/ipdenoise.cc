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
#include "rt_algo.h"
#include "ipdenoise.h"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace rtengine {

extern const Settings *settings;
using namespace procparams;

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


void calcautodn_info(const ProcParams *params, float &chaut, float &delta, int Nb, int levaut, float maxmax, float lumema, float chromina, int mode, int lissage, float redyel, float skinc, float nsknc)
{

    float reducdelta = 1.f;

    if (params->denoise.aggressive) {
        reducdelta = static_cast<float>(settings->nrhigh);
    }

    chaut = (chaut * Nb - maxmax) / (Nb - 1); //suppress maximum for chaut calcul

    if ((redyel > 5000.f || skinc > 1000.f) && nsknc < 0.4f  && chromina > 3000.f) {
        chaut *= 0.45f;    //reduct action in red zone, except skin for high / med chroma
    } else if ((redyel > 12000.f || skinc > 1200.f) && nsknc < 0.3f && chromina > 3000.f) {
        chaut *= 0.3f;
    }

    if (mode == 0 || mode == 2) { //Preview or Auto multizone
        if (chromina > 10000.f) {
            chaut *= 0.7f;    //decrease action for high chroma  (visible noise)
        } else if (chromina > 6000.f) {
            chaut *= 0.9f;
        } else if (chromina < 3000.f) {
            chaut *= 1.2f;    //increase action in low chroma==> 1.2  /==>2.0 ==> curve CC
        } else if (chromina < 2000.f) {
            chaut *= 1.5f;    //increase action in low chroma==> 1.5 / ==>2.7
        }

        if (lumema < 2500.f) {
            chaut *= 1.3f;    //increase action for low light
        } else if (lumema < 5000.f) {
            chaut *= 1.2f;
        } else if (lumema > 20000.f) {
            chaut *= 0.9f;    //decrease for high light
        }
    } else if (mode == 1) {//auto ==> less coefficient because interaction
        if (chromina > 10000.f) {
            chaut *= 0.8f;    //decrease action for high chroma  (visible noise)
        } else if (chromina > 6000.f) {
            chaut *= 0.9f;
        } else if (chromina < 3000.f) {
            chaut *= 1.5f;    //increase action in low chroma
        } else if (chromina < 2000.f) {
            chaut *= 2.2f;    //increase action in low chroma
        }

        if (lumema < 2500.f) {
            chaut *= 1.2f;    //increase action for low light
        } else if (lumema < 5000.f) {
            chaut *= 1.1f;
        } else if (lumema > 20000.f) {
            chaut *= 0.9f;    //decrease for high light
        }
    }

    if (levaut == 0) { //Low denoise
        if (chaut > 300.f) {
            chaut = 0.714286f * chaut + 85.71428f;
        }
    }

    delta = maxmax - chaut;
    delta *= reducdelta;

    if (lissage == 1 || lissage == 2) {
        if (chaut < 200.f && delta < 200.f) {
            delta *= 0.95f;
        } else if (chaut < 200.f && delta < 400.f) {
            delta *= 0.5f;
        } else if (chaut < 200.f && delta >= 400.f) {
            delta = 200.f;
        } else if (chaut < 400.f && delta < 400.f) {
            delta *= 0.4f;
        } else if (chaut < 400.f && delta >= 400.f) {
            delta = 120.f;
        } else if (chaut < 550.f) {
            delta *= 0.15f;
        } else if (chaut < 650.f) {
            delta *= 0.1f;
        } else { /*if (chaut >= 650.f)*/
            delta *= 0.07f;
        }

        if (mode == 0 || mode == 2) { //Preview or Auto multizone
            if (chromina < 6000.f) {
                delta *= 1.4f;    //increase maxi
            }

            if (lumema < 5000.f) {
                delta *= 1.4f;
            }
        } else if (mode == 1) { //Auto
            if (chromina < 6000.f) {
                delta *= 1.2f;    //increase maxi
            }

            if (lumema < 5000.f) {
                delta *= 1.2f;
            }
        }
    }

    if (lissage == 0) {
        if (chaut < 200.f && delta < 200.f) {
            delta *= 0.95f;
        } else if (chaut < 200.f && delta < 400.f) {
            delta *= 0.7f;
        } else if (chaut < 200.f && delta >= 400.f) {
            delta = 280.f;
        } else if (chaut < 400.f && delta < 400.f) {
            delta *= 0.6f;
        } else if (chaut < 400.f && delta >= 400.f) {
            delta = 200.f;
        } else if (chaut < 550.f) {
            delta *= 0.3f;
        } else if (chaut < 650.f) {
            delta *= 0.2f;
        } else { /*if (chaut >= 650.f)*/
            delta *= 0.15f;
        }

        if (mode == 0 || mode == 2) { //Preview or Auto multizone
            if (chromina < 6000.f) {
                delta *= 1.4f;    //increase maxi
            }

            if (lumema < 5000.f) {
                delta *= 1.4f;
            }
        } else if (mode == 1) { //Auto
            if (chromina < 6000.f) {
                delta *= 1.2f;    //increase maxi
            }

            if (lumema < 5000.f) {
                delta *= 1.2f;
            }
        }
    }

}


void RGB_denoise_infoGamCurve(const procparams::DenoiseParams & dnparams, bool isRAW, LUTf &gamcurve, float &gam, float &gamthresh, float &gamslope)
{
    gam = dnparams.gamma;
    gamthresh = 0.001f;

    if (!isRAW) {//reduce gamma under 1 for Lab mode ==> TIF and JPG
        if (gam < 1.9f) {
            gam = 1.f - (1.9f - gam) / 3.f;    //minimum gamma 0.7
        } else if (gam >= 1.9f && gam <= 3.f) {
            gam = (1.4f / 1.1f) * gam - 1.41818f;
        }
    }

    bool denoiseMethodRgb = (dnparams.colorSpace == procparams::DenoiseParams::ColorSpace::RGB);

    if (denoiseMethodRgb) {
        gamslope = exp(log(static_cast<double>(gamthresh)) / gam) / gamthresh;
        Color::gammaf2lut(gamcurve, gam, gamthresh, gamslope, 65535.f, 32768.f);
    } else {
        Color::gammanf2lut(gamcurve, gam, 65535.f, 32768.f);
    }
}


void RGB_denoise_info(ImProcData &im, Imagefloat * src, Imagefloat * provicalc, const bool isRAW, LUTf &gamcurve, float gam, float gamthresh, float gamslope, const procparams::DenoiseParams & dnparams, const double expcomp, float &chaut, int &Nb,  float &redaut, float &blueaut, float &maxredaut, float &maxblueaut, float &minredaut, float &minblueaut, float &chromina, float &sigma, float &lumema, float &sigma_L, float &redyel, float &skinc, float &nsknc)
{
    const ProcParams *params = im.params;
    double scale = im.scale;
    bool multiThread = im.multiThread;
    
    if (dnparams.chrominanceMethod != procparams::DenoiseParams::ChrominanceMethod::AUTOMATIC) {
        //nothing to do
        return;
    }

    int hei, wid;
    float** lumcalc;
    float** acalc;
    float** bcalc;
    hei = provicalc->getHeight();
    wid = provicalc->getWidth();
    TMatrix wprofi = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);

    const float wpi[3][3] = {
        {static_cast<float>(wprofi[0][0]), static_cast<float>(wprofi[0][1]), static_cast<float>(wprofi[0][2])},
        {static_cast<float>(wprofi[1][0]), static_cast<float>(wprofi[1][1]), static_cast<float>(wprofi[1][2])},
        {static_cast<float>(wprofi[2][0]), static_cast<float>(wprofi[2][1]), static_cast<float>(wprofi[2][2])}
    };

    lumcalc = new float*[hei];

    for (int i = 0; i < hei; ++i) {
        lumcalc[i] = new float[wid];
    }

    acalc = new float*[hei];

    for (int i = 0; i < hei; ++i) {
        acalc[i] = new float[wid];
    }

    bcalc = new float*[hei];

    for (int i = 0; i < hei; ++i) {
        bcalc[i] = new float[wid];
    }

#ifdef _OPENMP
    #pragma omp parallel for if (multiThread)
#endif

    for (int ii = 0; ii < hei; ++ii) {
        for (int jj = 0; jj < wid; ++jj) {
            float LLum, AAum, BBum;
            float RL = provicalc->r(ii, jj);
            float GL = provicalc->g(ii, jj);
            float BL = provicalc->b(ii, jj);
            // determine luminance for noisecurve
            float XL, YL, ZL;
            Color::rgbxyz(RL, GL, BL, XL, YL, ZL, wpi);
            Color::XYZ2Lab(XL, YL, ZL, LLum, AAum, BBum);
            lumcalc[ii][jj] = LLum;
            acalc[ii][jj] = AAum;
            bcalc[ii][jj] = BBum;
        }
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    const int imheight = src->getHeight(), imwidth = src->getWidth();

    bool denoiseMethodRgb = (dnparams.colorSpace == procparams::DenoiseParams::ColorSpace::RGB);

    const float gain = pow(2.0f, float(expcomp));

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

    //always no Tiles
    int kall = 0;
    rtengine::denoise::Tile_calc(tilesize, overlap, kall, imwidth, imheight, numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip);

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix(params->icm.workingProfile);
    const float wp[3][3] = {
        {static_cast<float>(wprof[0][0]), static_cast<float>(wprof[0][1]), static_cast<float>(wprof[0][2])},
        {static_cast<float>(wprof[1][0]), static_cast<float>(wprof[1][1]), static_cast<float>(wprof[1][2])},
        {static_cast<float>(wprof[2][0]), static_cast<float>(wprof[2][1]), static_cast<float>(wprof[2][2])}
    };

    float chau = 0.f;
    float chred = 0.f;
    float chblue = 0.f;
    float maxchred = 0.f;
    float maxchblue = 0.f;
    float minchred = 100000000.f;
    float minchblue = 100000000.f;
    int nb = 0;
    int comptlevel = 0;

    for (int tiletop = 0; tiletop < imheight; tiletop += tileHskip) {
        for (int tileleft = 0; tileleft < imwidth; tileleft += tileWskip) {

            int tileright = MIN(imwidth, tileleft + tilewidth);
            int tilebottom = MIN(imheight, tiletop + tileheight);
            int width  = tileright - tileleft;
            int height = tilebottom - tiletop;
            LabImage * labdn = new LabImage(width, height);
            float** noisevarlum = new float*[(height + 1) / 2];

            for (int i = 0; i < (height + 1) / 2; ++i) {
                noisevarlum[i] = new float[(width + 1) / 2];
            }

            float** noisevarchrom = new float*[(height + 1) / 2];

            for (int i = 0; i < (height + 1) / 2; ++i) {
                noisevarchrom[i] = new float[(width + 1) / 2];
            }

            float** noisevarhue = new float*[(height + 1) / 2];

            for (int i = 0; i < (height + 1) / 2; ++i) {
                noisevarhue[i] = new float[(width + 1) / 2];
            }

            float realred, realblue;
            float interm_med = static_cast<float>(dnparams.chrominance) / 10.0;
            float intermred, intermblue;

            if (dnparams.chrominanceRedGreen > 0.) {
                intermred = (dnparams.chrominanceRedGreen / 10.);
            } else {
                intermred = static_cast<float>(dnparams.chrominanceRedGreen) / 7.0;     //increase slower than linear for more sensit
            }

            if (dnparams.chrominanceBlueYellow > 0.) {
                intermblue = (dnparams.chrominanceBlueYellow / 10.);
            } else {
                intermblue = static_cast<float>(dnparams.chrominanceBlueYellow) / 7.0;     //increase slower than linear for more sensit
            }

            realred = interm_med + intermred;

            if (realred < 0.f) {
                realred = 0.001f;
            }

            realblue = interm_med + intermblue;

            if (realblue < 0.f) {
                realblue = 0.001f;
            }

            //fill tile from image; convert RGB to "luma/chroma"

            if (isRAW) {//image is raw; use channel differences for chroma channels
#ifdef _OPENMP
                #pragma omp parallel for if (multiThread)
#endif

                for (int i = tiletop; i < tilebottom; i += 2) {
                    int i1 = i - tiletop;
#ifdef __SSE2__
                    __m128 aNv, bNv;
                    __m128 c100v = _mm_set1_ps(100.f);
                    int j;

                    for (j = tileleft; j < tileright - 7; j += 8) {
                        int j1 = j - tileleft;
                        aNv = LVFU(acalc[i >> 1][j >> 1]);
                        bNv = LVFU(bcalc[i >> 1][j >> 1]);
                        _mm_storeu_ps(&noisevarhue[i1 >> 1][j1 >> 1], xatan2f(bNv, aNv));
                        _mm_storeu_ps(&noisevarchrom[i1 >> 1][j1 >> 1], vmaxf(vsqrtf(SQRV(aNv) + SQRV(bNv)),c100v));
                    }

                    for (; j < tileright; j += 2) {
                        int j1 = j - tileleft;
                        float aN = acalc[i >> 1][j >> 1];
                        float bN = bcalc[i >> 1][j >> 1];
                        float cN = sqrtf(SQR(aN) + SQR(bN));
                        noisevarhue[i1 >> 1][j1 >> 1] = xatan2f(bN, aN);

                        if (cN < 100.f) {
                            cN = 100.f;    //avoid divided by zero
                        }

                        noisevarchrom[i1 >> 1][j1 >> 1] = cN;
                    }

#else

                    for (int j = tileleft; j < tileright; j += 2) {
                        int j1 = j - tileleft;
                        float aN = acalc[i >> 1][j >> 1];
                        float bN = bcalc[i >> 1][j >> 1];
                        float cN = sqrtf(SQR(aN) + SQR(bN));
                        float hN = xatan2f(bN, aN);

                        if (cN < 100.f) {
                            cN = 100.f;    //avoid divided by zero
                        }

                        noisevarchrom[i1 >> 1][j1 >> 1] = cN;
                        noisevarhue[i1 >> 1][j1 >> 1] = hN;
                    }

#endif
                }

#ifdef _OPENMP
                #pragma omp parallel for if (multiThread)
#endif

                for (int i = tiletop; i < tilebottom; i += 2) {
                    int i1 = i - tiletop;

                    for (int j = tileleft; j < tileright; j += 2) {
                        int j1 = j - tileleft;
                        float Llum = lumcalc[i >> 1][j >> 1];
                        Llum = Llum < 2.f ? 2.f : Llum; //avoid divided by zero ?
                        Llum = Llum > 32768.f ? 32768.f : Llum; // not strictly necessary
                        noisevarlum[i1 >> 1][j1 >> 1] = Llum;
                    }
                }

                if (!denoiseMethodRgb) { //lab mode, modification Jacques feb 2013 and july 2014

#ifdef _OPENMP
                    #pragma omp parallel for if (multiThread)
#endif

                    for (int i = tiletop; i < tilebottom; ++i) {
                        int i1 = i - tiletop;

                        for (int j = tileleft; j < tileright; ++j) {
                            int j1 = j - tileleft;
                            float R_ = gain * src->r(i, j);
                            float G_ = gain * src->g(i, j);
                            float B_ = gain * src->b(i, j);

                            R_ = Color::denoiseIGammaTab[R_];
                            G_ = Color::denoiseIGammaTab[G_];
                            B_ = Color::denoiseIGammaTab[B_];

                            //apply gamma noise standard (slider)
                            R_ = R_ < 65535.f ? gamcurve[R_] : (Color::gammanf(R_ / 65535.f, gam) * 32768.f);
                            G_ = G_ < 65535.f ? gamcurve[G_] : (Color::gammanf(G_ / 65535.f, gam) * 32768.f);
                            B_ = B_ < 65535.f ? gamcurve[B_] : (Color::gammanf(B_ / 65535.f, gam) * 32768.f);
                            //true conversion xyz=>Lab
                            float X, Y, Z;
                            Color::rgbxyz(R_, G_, B_, X, Y, Z, wp);

                            //convert to Lab
                            float L, a, b;
                            Color::XYZ2Lab(X, Y, Z, L, a, b);

                            labdn->a[i1][j1] = a;
                            labdn->b[i1][j1] = b;
                        }
                    }
                } else { //RGB mode

                    for (int i = tiletop/*, i1=0*/; i < tilebottom; ++i/*, ++i1*/) {
                        int i1 = i - tiletop;

                        for (int j = tileleft/*, j1=0*/; j < tileright; ++j/*, ++j1*/) {
                            int j1 = j - tileleft;

                            float X = gain * src->r(i, j);
                            float Y = gain * src->g(i, j);
                            float Z = gain * src->b(i, j);

                            X = X < 65535.f ? gamcurve[X] : (Color::gammaf(X / 65535.f, gam, gamthresh, gamslope) * 32768.f);
                            Y = Y < 65535.f ? gamcurve[Y] : (Color::gammaf(Y / 65535.f, gam, gamthresh, gamslope) * 32768.f);
                            Z = Z < 65535.f ? gamcurve[Z] : (Color::gammaf(Z / 65535.f, gam, gamthresh, gamslope) * 32768.f);

                            // labdn->a[i1][j1] = (X - Y);
                            // labdn->b[i1][j1] = (Y - Z);
                            float l, u, v;
                            Color::rgb2yuv(X, Y, Z, l, u, v, wp);
                            labdn->a[i1][j1] = v;
                            labdn->b[i1][j1] = u;
                        }
                    }
                }

            } else {//image is not raw; use Lab parametrization
                for (int i = tiletop/*, i1=0*/; i < tilebottom; ++i/*, ++i1*/) {
                    int i1 = i - tiletop;

                    for (int j = tileleft/*, j1=0*/; j < tileright; ++j/*, ++j1*/) {
                        int j1 = j - tileleft;
                        float L, a, b;
                        float rLum = src->r(i, j) ; //for luminance denoise curve
                        float gLum = src->g(i, j) ;
                        float bLum = src->b(i, j) ;

                        //use gamma sRGB, not good if TIF (JPG) Output profil not with gamma sRGB  (eg : gamma =1.0, or 1.8...)
                        //very difficult to solve !
                        // solution ==> save TIF with gamma sRGB and re open
                        float rtmp = Color::igammatab_srgb[ src->r(i, j) ];
                        float gtmp = Color::igammatab_srgb[ src->g(i, j) ];
                        float btmp = Color::igammatab_srgb[ src->b(i, j) ];
                        //modification Jacques feb 2013
                        // gamma slider different from raw
                        rtmp = rtmp < 65535.f ? gamcurve[rtmp] : (Color::gammanf(rtmp / 65535.f, gam) * 32768.f);
                        gtmp = gtmp < 65535.f ? gamcurve[gtmp] : (Color::gammanf(gtmp / 65535.f, gam) * 32768.f);
                        btmp = btmp < 65535.f ? gamcurve[btmp] : (Color::gammanf(btmp / 65535.f, gam) * 32768.f);

                        float X, Y, Z;
                        Color::rgbxyz(rtmp, gtmp, btmp, X, Y, Z, wp);

                        //convert Lab
                        Color::XYZ2Lab(X, Y, Z, L, a, b);

                        if (((i1 | j1) & 1) == 0) {
                            float Llum, alum, blum;
                            float XL, YL, ZL;
                            Color::rgbxyz(rLum, gLum, bLum, XL, YL, ZL, wp);
                            Color::XYZ2Lab(XL, YL, ZL, Llum, alum, blum);
                            float kN = Llum;

                            if (kN < 2.f) {
                                kN = 2.f;
                            }

                            if (kN > 32768.f) {
                                kN = 32768.f;
                            }

                            noisevarlum[i1 >> 1][j1 >> 1] = kN;
                            float aN = alum;
                            float bN = blum;
                            float hN = xatan2f(bN, aN);
                            float cN = sqrt(SQR(aN) + SQR(bN));

                            if (cN < 100.f) {
                                cN = 100.f;    //avoid divided by zero
                            }

                            noisevarchrom[i1 >> 1][j1 >> 1] = cN;
                            noisevarhue[i1 >> 1][j1 >> 1] = hN;
                        }

                        labdn->a[i1][j1] = a;
                        labdn->b[i1][j1] = b;
                    }
                }
            }

            int datalen = labdn->W * labdn->H;

            //now perform basic wavelet denoise
            //last two arguments of wavelet decomposition are max number of wavelet decomposition levels;
            //and whether to subsample the image after wavelet filtering.  Subsampling is coded as
            //binary 1 or 0 for each level, eg subsampling = 0 means no subsampling, 1 means subsample
            //the first level only, 7 means subsample the first three levels, etc.

            wavelet_decomposition* adecomp;
            wavelet_decomposition* bdecomp;

            int schoice = 0;//shrink method

            if (dnparams.aggressive) {
                schoice = 2;
            }

            const int levwav = max(2, int(5 - std::ceil(std::log(scale))));
#ifdef _OPENMP
            #pragma omp parallel sections if (multiThread)
#endif
            {
#ifdef _OPENMP
                #pragma omp section
#endif
                {
                    adecomp = new wavelet_decomposition(labdn->data + datalen, labdn->W, labdn->H, levwav, 1);
                }
#ifdef _OPENMP
                #pragma omp section
#endif
                {
                    bdecomp = new wavelet_decomposition(labdn->data + 2 * datalen, labdn->W, labdn->H, levwav, 1);
                }
            }

            if (comptlevel == 0) {
                denoise::WaveletDenoiseAll_info(
                    levwav,
                    *adecomp,
                    *bdecomp,
                    noisevarlum,
                    noisevarchrom,
                    noisevarhue,
                    chaut,
                    Nb,
                    redaut,
                    blueaut,
                    maxredaut,
                    maxblueaut,
                    minredaut,
                    minblueaut,
                    schoice,
                    chromina,
                    sigma,
                    lumema,
                    sigma_L,
                    redyel,
                    skinc,
                    nsknc,
                    maxchred,
                    maxchblue,
                    minchred,
                    minchblue,
                    nb,
                    chau,
                    chred,
                    chblue,
                    denoiseMethodRgb
                ); // Enhance mode
            }

            comptlevel += 1;
            delete adecomp;
            delete bdecomp;
            delete labdn;

            for (int i = 0; i < (height + 1) / 2; ++i) {
                delete[] noisevarlum[i];
            }

            delete[] noisevarlum;

            for (int i = 0; i < (height + 1) / 2; ++i) {
                delete[] noisevarchrom[i];
            }

            delete[] noisevarchrom;

            for (int i = 0; i < (height + 1) / 2; ++i) {
                delete[] noisevarhue[i];
            }

            delete[] noisevarhue;

        }//end of tile row
    }//end of tile loop

    for (int i = 0; i < hei; ++i) {
        delete[] lumcalc[i];
    }

    delete[] lumcalc;

    for (int i = 0; i < hei; ++i) {
        delete[] acalc[i];
    }

    delete[] acalc;

    for (int i = 0; i < hei; ++i) {
        delete[] bcalc[i];
    }

    delete[] bcalc;

#undef TS
//#undef fTS
#undef offset
#undef epsilon

} // End of main RGB_denoise


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
    
    denoise::Tile_calc(tilesize, overlap, kall, widIm, heiIm, numtiles_W, numtiles_H, tilewidth, tileheight, tileWskip, tileHskip);
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
                    ImProcData im(params, scale, multiThread);
                    RGB_denoise_info(im, origCropPart, provicalc, imgsrc->isRAW(), gamcurve, gam, gamthresh, gamslope, dnparams, imgsrc->getDirPyrDenoiseExpComp(), chaut, nb, redaut, blueaut, maxredaut, maxblueaut, minredaut, minblueaut, chromina, sigma, lumema, sigma_L, redyel, skinc, nsknc);

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
            calcautodn_info(params, store.ch_M[k], delta[k], Nb[k], levaut, maxmax, lumL[k], chromC[k], mode, lissage, ry[k], sk[k], pcsk[k]);
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
    if (!dnparams.enabled) {
        return;
    }
    
    procparams::DenoiseParams denoiseParams = dnparams;
    NoiseCurve noiseLCurve;
    NoiseCurve noiseCCurve;

    const int W = img->getWidth();
    const int H = img->getHeight();

    Imagefloat *calclum = nullptr;
    {
        const int fw = W;
        const int fh = H;
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

    ImProcData im(params, scale, multiThread);
    denoise::RGB_denoise(im, 0, img, img, calclum, dnstore.ch_M, dnstore.max_r, dnstore.max_b, imgsrc->isRAW(), denoiseParams, imgsrc->getDirPyrDenoiseExpComp(), noiseLCurve, noiseCCurve, nresi, highresi);

    if (denoiseParams.smoothingEnabled) {
        if (denoiseParams.smoothingMethod == DenoiseParams::SmoothingMethod::MEDIAN) {
            Median median_type[3];
            int median_iterations[3];
            for (int i = 0; i < 3; ++i) {
                median_type[i] = Median(int(denoiseParams.medianType));
                median_iterations[i] = denoiseParams.medianIterations;
            }
            if (denoiseParams.medianMethod != DenoiseParams::MedianMethod::RGB) {
                if (denoiseParams.colorSpace == DenoiseParams::ColorSpace::LAB) {
                    img->setMode(Imagefloat::Mode::LAB, multiThread);
                } else {
                    img->setMode(Imagefloat::Mode::YUV, multiThread);
                }
                switch (denoiseParams.medianMethod) {
                case DenoiseParams::MedianMethod::LUMINANCE:
                    median_iterations[1] = median_iterations[2] = 0;
                    break;
                case DenoiseParams::MedianMethod::CHROMINANCE:
                    median_iterations[0] = 0;
                    break;
                case DenoiseParams::MedianMethod::LAB_WEIGHTED:
                    switch (median_type[0]) {
                    case Median::TYPE_3X3_SOFT:
                        break;
                    case Median::TYPE_3X3_STRONG:
                    case Median::TYPE_5X5_SOFT:
                        median_type[0] = Median::TYPE_3X3_SOFT;
                        break;
                    case Median::TYPE_5X5_STRONG:
                    case Median::TYPE_7X7:
                        median_type[0] = Median::TYPE_3X3_STRONG;
                        break;
                    case Median::TYPE_9X9:
                        median_type[0] = Median::TYPE_5X5_SOFT;
                        break;
                    }
                    break;
                default:
                    break;
                }
            }
            int num_threads = 1;
#ifdef _OPENMP
            if (multiThread) {
                num_threads = omp_get_max_threads();
            }
#endif
            array2D<float> buf(W, H);
            float **channels[3] = { img->g.ptrs, img->r.ptrs, img->b.ptrs };
            for (int c = 0; c < 3; ++c) {
                if (median_iterations[c] > 0) {
                    denoise::Median_Denoise(channels[c], channels[c], W, H, median_type[c], median_iterations[c], num_threads, buf);
                }
            }
            img->setMode(Imagefloat::Mode::RGB, multiThread);
        } else {
            denoise::denoiseGuidedSmoothing(im, img);
        }
    }
}


} // namespace rtengine
