/* -*- C++ -*-
 *
 *  This file is part of ART.
 *
 *  Copyright 2019 Alberto Griggio <alberto.griggio@gmail.com>
 *
 *  ART is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ART is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ART.  If not, see <http://www.gnu.org/licenses/>.
 */
// extracted and datapted from ImProcFunctions::rgbProc (improcfun.cc) of
// RawTherapee

#include "improcfun.h"
#include "curves.h"
#include "color.h"

namespace rtengine {

/** @brief Compute the B&W constants for the Black and White processing and its GUI
 * @param setting main mode
 * @param filter string of the filter effect to use
 * @param algo choice between linear and special for OYCPM colors
 * @param mixerRed red channel value of the channel mixer [-100 ; +200]
 * @param mixerGreen green channel value of the channel mixer [-100 ; +200]
 * @param mixerBlue blue channel value of the channel mixer [-100 ; +200]
 * @param mixerOrange orange channel value of the channel mixer [-100 ; +200]
 * @param mixerYellow yellow channel value of the channel mixer [-100 ; +200]
 * @param mixerCyan cyan channel value of the channel mixer [-100 ; +200]
 * @param mixerPurple purple channel value of the channel mixer [-100 ; +200]
 * @param mixerMagenta magenta channel value of the channel mixer [-100 ; +200]
 * @param autoc automatic mode of the channel mixer
 * @param complement adjust complementary channel
 * @param kcorec in absolute mode, value to correct the mixer [1 ; 3], usually near 1 (return value)
 * @param rrm red channel of the mixer (return value)
 * @param ggm green channel of the mixer (return value)
 * @param bbm blue channel of the mixer (return value)
 */
void computeBWMixerConstants(const Glib::ustring &setting, const Glib::ustring &filter, const Glib::ustring &algo,
                             float &filcor, float &mixerRed, float &mixerGreen,
                             float &mixerBlue, 
                             float &kcorec, double &rrm, double &ggm, double &bbm)
{
    float somm;
    float som = mixerRed + mixerGreen + mixerBlue;

    if (som >= 0.f && som < 1.f) {
        som = 1.f;
    }

    if (som < 0.f && som > -1.f) {
        som = -1.f;
    }


    // rM = mixerRed, gM = mixerGreen, bM = mixerBlue !
    //presets
    if (setting == "RGB-Abs" || setting == "ROYGCBPM-Abs") {
        kcorec = som / 100.f;
    }

    if (setting == "NormalContrast") {
        mixerRed = 43.f ;
        mixerGreen = 33.f;
        mixerBlue = 30.f;
    } else if (setting == "Panchromatic")      {
        mixerRed = 33.3f;
        mixerGreen = 33.3f;
        mixerBlue = 33.3f;
    } else if (setting == "HyperPanchromatic") {
        mixerRed = 41.f ;
        mixerGreen = 25.f;
        mixerBlue = 34.f;
    } else if (setting == "LowSensitivity")    {
        mixerRed = 27.f ;
        mixerGreen = 27.f;
        mixerBlue = 46.f;
    } else if (setting == "HighSensitivity")   {
        mixerRed = 30.f ;
        mixerGreen = 28.f;
        mixerBlue = 42.f;
    } else if (setting == "Orthochromatic")    {
        mixerRed = 0.f  ;
        mixerGreen = 42.f;
        mixerBlue = 58.f;
    } else if (setting == "HighContrast")      {
        mixerRed = 40.f ;
        mixerGreen = 34.f;
        mixerBlue = 60.f;
    } else if (setting == "Luminance")         {
        mixerRed = 30.f ;
        mixerGreen = 59.f;
        mixerBlue = 11.f;
    } else if (setting == "Landscape")         {
        mixerRed = 66.f ;
        mixerGreen = 24.f;
        mixerBlue = 10.f;
    } else if (setting == "Portrait")          {
        mixerRed = 54.f ;
        mixerGreen = 44.f;
        mixerBlue = 12.f;
    } else if (setting == "InfraRed")          {
        mixerRed = -40.f;
        mixerGreen = 200.f;
        mixerBlue = -17.f;
    }

    rrm = mixerRed;
    ggm = mixerGreen;
    bbm = mixerBlue;

    somm = mixerRed + mixerGreen + mixerBlue;

    if(somm >= 0.f && somm < 1.f) {
        somm = 1.f;
    }

    if(somm < 0.f && somm > -1.f) {
        somm = -1.f;
    }

    mixerRed = mixerRed / somm;
    mixerGreen = mixerGreen / somm;
    mixerBlue = mixerBlue / somm;

    //Color filters
    float filred, filgreen, filblue;
    filred = 1.f;
    filgreen = 1.f;
    filblue = 1.f;
    filcor = 1.f;

    if          (filter == "None")        {
        filred = 1.f;
        filgreen = 1.f;
        filblue = 1.f;
        filcor = 1.f;
    } else if     (filter == "Red")         {
        filred = 1.f;
        filgreen = 0.05f;
        filblue = 0.f;
        filcor = 1.08f;
    } else if     (filter == "Orange")      {
        filred = 1.f;
        filgreen = 0.6f;
        filblue = 0.f;
        filcor = 1.35f;
    } else if     (filter == "Yellow")      {
        filred = 1.f;
        filgreen = 1.f;
        filblue = 0.05f;
        filcor = 1.23f;
    } else if     (filter == "YellowGreen") {
        filred = 0.6f;
        filgreen = 1.f;
        filblue = 0.3f;
        filcor = 1.32f;
    } else if     (filter == "Green")       {
        filred = 0.2f;
        filgreen = 1.f;
        filblue = 0.3f;
        filcor = 1.41f;
    } else if     (filter == "Cyan")        {
        filred = 0.05f;
        filgreen = 1.f;
        filblue = 1.f;
        filcor = 1.23f;
    } else if     (filter == "Blue")        {
        filred = 0.f;
        filgreen = 0.05f;
        filblue = 1.f;
        filcor = 1.20f;
    } else if     (filter == "Purple")      {
        filred = 1.f;
        filgreen = 0.05f;
        filblue = 1.f;
        filcor = 1.23f;
    }

    mixerRed   = mixerRed * filred;
    mixerGreen = mixerGreen * filgreen;
    mixerBlue  = mixerBlue  * filblue;

    if(mixerRed + mixerGreen + mixerBlue == 0) {
        mixerRed += 1.f;
    }

    mixerRed   = filcor * mixerRed   / (mixerRed + mixerGreen + mixerBlue);
    mixerGreen = filcor * mixerGreen / (mixerRed + mixerGreen + mixerBlue);
    mixerBlue  = filcor * mixerBlue  / (mixerRed + mixerGreen + mixerBlue);

    if(filter != "None") {
        som = mixerRed + mixerGreen + mixerBlue;

        if(som >= 0.f && som < 1.f) {
            som = 1.f;
        }

        if(som < 0.f && som > -1.f) {
            som = -1.f;
        }

        if(setting == "RGB-Abs") {
            kcorec = kcorec * som;
        }
    }
}

namespace {

/**
 * @brief Used by Black and White to correct gamma for each channel red, green and blue channel
 * @param r red channel input and output value [0 ; 65535]
 * @param g green channel input and output value [0 ; 65535]
 * @param b blue channel input and output value [0 ; 65535]
 * @param gammabwr gamma value for red channel [>0]
 * @param gammabwg gamma value for red channel [>0]
 * @param gammabwb gamma value for red channel [>0]
 */
#ifdef __SSE2__
void trcGammaBW(float &r, float &g, float &b, float gammabwr, float gammabwg, float gammabwb)
{
    // correct gamma for black and white image : pseudo TRC curve of ICC profile
    vfloat rgbv = _mm_set_ps(0.f, r, r, r); // input channel is always r
    vfloat gammabwv = _mm_set_ps(0.f, gammabwb, gammabwg, gammabwr);
    vfloat c65535v = F2V(65535.f);
    rgbv /= c65535v;
    rgbv = vmaxf(rgbv, ZEROV);
    rgbv = pow_F(rgbv, gammabwv);
    rgbv *= c65535v;
    float temp[4] ALIGNED16;
    STVF(temp[0], rgbv);
    r = temp[0];
    g = temp[1];
    b = temp[2];
}

void trcGammaBWRow(float *r, float *g, float *b, int width, float gammabwr, float gammabwg, float gammabwb)
{
    // correct gamma for black and white image : pseudo TRC curve of ICC profile
    vfloat c65535v = F2V(65535.f);
    vfloat gammabwrv = F2V(gammabwr);
    vfloat gammabwgv = F2V(gammabwg);
    vfloat gammabwbv = F2V(gammabwb);
    int i = 0;
    for(; i < width - 3; i += 4 ) {
        vfloat inv = _mm_loadu_ps(&r[i]); // input channel is always r
        inv /= c65535v;
        inv = vmaxf(inv, ZEROV);
        vfloat rv = pow_F(inv, gammabwrv);
        vfloat gv = pow_F(inv, gammabwgv);
        vfloat bv = pow_F(inv, gammabwbv);
        rv *= c65535v;
        gv *= c65535v;
        bv *= c65535v;
        _mm_storeu_ps(&r[i], rv);
        _mm_storeu_ps(&g[i], gv);
        _mm_storeu_ps(&b[i], bv);
    }
    for(; i < width; i++) {
        trcGammaBW(r[i], g[i], b[i], gammabwr, gammabwg, gammabwb);
    }
}

#else // __SSE2__
void trcGammaBW(float &r, float &g, float &b, float gammabwr, float gammabwg, float gammabwb)
{
    // correct gamma for black and white image : pseudo TRC curve of ICC profile
    float in = r; // input channel is always r
    in /= 65535.0f;
    in = max(in, 0.f);
    b = pow_F (in, gammabwb);
    b *= 65535.0f;
    r = pow_F (in, gammabwr);
    r *= 65535.0f;
    g = pow_F (in, gammabwg);
    g *= 65535.0f;
}
#endif // __SSE2__

} // namespace


void ImProcFunctions::blackAndWhite(Imagefloat *img)
{
    if (!params->blackwhite.enabled) {
        return;
    }

    img->setMode(Imagefloat::Mode::RGB, multiThread);

    float bwr = float(params->blackwhite.mixerRed);
    float bwg = float(params->blackwhite.mixerGreen);
    float bwb = float(params->blackwhite.mixerBlue);
    float bwrgam = float(params->blackwhite.gammaRed);
    float bwggam = float(params->blackwhite.gammaGreen);
    float bwbgam = float(params->blackwhite.gammaBlue);

    float gamvalr = 125.f;
    float gamvalg = 125.f;
    float gamvalb = 125.f;

    if (bwrgam < 0) {
        gamvalr = 100.f;
    }

    if (bwggam < 0) {
        gamvalg = 100.f;
    }

    if (bwbgam < 0) {
        gamvalb = 100.f;
    }

    float gammabwr = 1.f;
    float gammabwg = 1.f;
    float gammabwb = 1.f;
    {
        gammabwr = 1.f - bwrgam / gamvalr;
        gammabwg = 1.f - bwggam / gamvalg;
        gammabwb = 1.f - bwbgam / gamvalb;
    }
    bool hasgammabw = gammabwr != 1.f || gammabwg != 1.f || gammabwb != 1.f;
    
    const int W = img->getWidth();
    const int H = img->getHeight();

    float kcorec = 1.f;
    float filcor;
    double rrm, ggm, bbm;
    computeBWMixerConstants(params->blackwhite.setting, params->blackwhite.filter, "", filcor, bwr, bwg, bwb, kcorec, rrm, ggm, bbm);

#ifdef _OPENMP
#   pragma omp parallel for if (multiThread)
#endif
    for (int y = 0; y < H; ++y) {
#ifdef __SSE2__
        if (hasgammabw) {
            trcGammaBWRow(img->r(y), img->g(y), img->b(y), W, gammabwr, gammabwg, gammabwb);
        }
#endif
        for (int x = 0; x < W; ++x) {
#ifndef __SSE2__
            if (hasgammabw) {
                trcGammaBW(img->r(y, x), img->g(y, x), img->b(y, x), gammabwr, gammabwg, gammabwb);
            }
#endif
            img->r(y, x) = img->g(y, x) = img->b(y, x) = ((bwr * img->r(y, x) + bwg * img->g(y, x) + bwb * img->b(y, x)) * kcorec);
        }
    }

    if (params->blackwhite.colorCast.getBottom() > 0) {
        // apply color cast
        float s = float(params->blackwhite.colorCast.getBottom()) / 100.f;
        float h = float(params->blackwhite.colorCast.getTop()) / 180.f * rtengine::RT_PI;
        float u, v;
        Color::hsl2yuv(h, s, u, v);
        img->setMode(Imagefloat::Mode::YUV, multiThread);

#ifdef __SSE2__
        vfloat uv = F2V(u);
        vfloat vv = F2V(v);
#endif

#ifdef _OPENMP
#       pragma omp parallel for if (multiThread)
#endif
        for (int y = 0; y < H; ++y) {
            int x = 0;
#ifdef __SSE2__
            for (; x < W - 3; x += 4) {
                vfloat Yv = LVF(img->g(y, x));
                STVF(img->b(y, x), LVF(img->b(y, x)) + Yv * uv);
                STVF(img->r(y, x), LVF(img->r(y, x)) + Yv * vv);
            }
#endif
            for (; x < W; ++x) {
                float Y = img->g(y, x);
                img->b(y, x) += Y * u;
                img->r(y, x) += Y * v;
            }
        }
    }
}

} // namespace rtengine
