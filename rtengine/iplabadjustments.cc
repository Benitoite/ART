/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright 2019 Alberto Griggio <alberto.griggio@gmail.com>
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
#include "curves.h"
#include "settings.h"
#include "mytime.h"

namespace rtengine {

extern const Settings* settings;

namespace {

void complexLCurve(double br, double contr, const std::vector<double>& curvePoints,
                   const LUTu & histogram, LUTf & outCurve,
                   LUTu & outBeforeCCurveHistogram, int skip, bool & utili)
{

    utili = false;
    // clear array that stores histogram valid before applying the custom curve
    if (outBeforeCCurveHistogram) {
        outBeforeCCurveHistogram.clear();
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // tone curve base. a: slope (from exp.comp.), b: black, def_mul: max. x value (can be>1), hr,sr: highlight,shadow recovery
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // check if brightness curve is needed
    if (br > 0.00001 || br < -0.00001) {
        utili = true;

        std::vector<double> brightcurvePoints;
        brightcurvePoints.resize(9);
        brightcurvePoints.at(0) = double(DCT_NURBS);

        brightcurvePoints.at(1) = 0.; // black point.  Value in [0 ; 1] range
        brightcurvePoints.at(2) = 0.; // black point.  Value in [0 ; 1] range

        if (br > 0) {
            brightcurvePoints.at(3) = 0.1; // toe point
            brightcurvePoints.at(4) = 0.1 + br / 150.0; //value at toe point

            brightcurvePoints.at(5) = 0.7; // shoulder point
            brightcurvePoints.at(6) = min(1.0, 0.7 + br / 300.0); //value at shoulder point
        } else {
            brightcurvePoints.at(3) = 0.1 - br / 150.0; // toe point
            brightcurvePoints.at(4) = 0.1; // value at toe point

            brightcurvePoints.at(5) = min(1.0, 0.7 - br / 300.0); // shoulder point
            brightcurvePoints.at(6) = 0.7; // value at shoulder point
        }

        brightcurvePoints.at(7) = 1.; // white point
        brightcurvePoints.at(8) = 1.; // value at white point

        DiagonalCurve brightcurve(brightcurvePoints, CURVES_MIN_POLY_POINTS / skip);
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        // Applying brightness curve
        for (int i = 0; i < 32768; i++) { // L values range up to 32767, higher values are for highlight overflow

            // change to [0,1] range
            float val = (float)i / 32767.0;

            // apply brightness curve
            val = brightcurve.getVal (val);

            // store result in a temporary array
            outCurve[i] = LIM01(val);
        }

    } else {
        outCurve.makeIdentity(32767.f);
    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // check if contrast curve is needed
    if (contr > 0.00001 || contr < -0.00001) {
        utili = true;

        // compute mean luminance of the image with the curve applied
        int sum = 0;
        float avg = 0;

        for (int i = 0; i < 32768; i++) {
            avg += outCurve[i] * histogram[i];
            sum += histogram[i];
        }

        std::vector<double> contrastcurvePoints;

        if(sum) {
            avg /= sum;

            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            contrastcurvePoints.resize(9);
            contrastcurvePoints.at(0) = double(DCT_NURBS);

            contrastcurvePoints.at(1) = 0.; // black point.  Value in [0 ; 1] range
            contrastcurvePoints.at(2) = 0.; // black point.  Value in [0 ; 1] range

            contrastcurvePoints.at(3) = avg - avg * (0.6 - contr / 250.0); // toe point
            contrastcurvePoints.at(4) = avg - avg * (0.6 + contr / 250.0); // value at toe point

            contrastcurvePoints.at(5) = avg + (1 - avg) * (0.6 - contr / 250.0); // shoulder point
            contrastcurvePoints.at(6) = avg + (1 - avg) * (0.6 + contr / 250.0); // value at shoulder point

            contrastcurvePoints.at(7) = 1.; // white point
            contrastcurvePoints.at(8) = 1.; // value at white point

            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        } else {
            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            // sum has an invalid value (next to 0, producing a division by zero, so we create a fake contrast curve, producing a white image
            contrastcurvePoints.resize(5);
            contrastcurvePoints.at(0) = double(DCT_NURBS);

            contrastcurvePoints.at(1) = 0.; // black point.  Value in [0 ; 1] range
            contrastcurvePoints.at(2) = 1.; // black point.  Value in [0 ; 1] range

            contrastcurvePoints.at(3) = 1.; // white point
            contrastcurvePoints.at(4) = 1.; // value at white point

            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        }

        DiagonalCurve contrastcurve(contrastcurvePoints, CURVES_MIN_POLY_POINTS / skip);

        // apply contrast enhancement
        for (int i = 0; i < 32768; i++) {
            outCurve[i] = contrastcurve.getVal (outCurve[i]);
        }

    }

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    // create a curve if needed
    std::unique_ptr<DiagonalCurve> tcurve;
    bool histNeeded = false;

    if (!curvePoints.empty() && curvePoints[0] != 0) {
        tcurve = std::unique_ptr<DiagonalCurve>(new DiagonalCurve (curvePoints, CURVES_MIN_POLY_POINTS / skip));

        if (outBeforeCCurveHistogram) {
            histNeeded = true;
        }
    }

    if (tcurve && tcurve->isIdentity()) {
        tcurve = nullptr;
    }

    if (tcurve) {
        utili = true; //if active

        // L values go up to 32767, last stop is for highlight overflow
        for (int i = 0; i < 32768; i++) {
            float val;

            if (histNeeded) {
                float hval = outCurve[i];
                int hi = (int)(255.f * hval);
                outBeforeCCurveHistogram[hi] += histogram[i] ;
            }

            // apply custom/parametric/NURBS curve, if any
            val = tcurve->getVal (outCurve[i]);

            outCurve[i] = (32767.f * val);
        }
    } else {

        // Skip the slow getval method if no curve is used (or an identity curve)
        // L values go up to 32767, last stop is for highlight overflow
        if(histNeeded) {
            histogram.compressTo(outBeforeCCurveHistogram, 32768, outCurve);
        }

        outCurve *= 32767.f;

    }

    for (int i = 32768; i < 32770; i++) { // set last two elements of lut to 32768 and 32769 to allow linear interpolation
        outCurve[i] = (float)i;
    }
}


// add curve Lab : C=f(L)
void curveCL( bool & clcutili, const std::vector<double>& clcurvePoints, LUTf & clCurve, int skip)
{
    clcutili = false;
    std::unique_ptr<DiagonalCurve> dCurve;

    if (!clcurvePoints.empty() && clcurvePoints[0] != 0) {
        dCurve = std::unique_ptr<DiagonalCurve>(new DiagonalCurve(clcurvePoints, CURVES_MIN_POLY_POINTS / skip));

        if (dCurve && !dCurve->isIdentity()) {
            clcutili = true;
        }
    }

    CurveFactory::fillCurveArray(dCurve.get(), clCurve, skip, clcutili);
}


void complexsgnCurve (bool & autili,  bool & butili, bool & ccutili, bool & cclutili,
                                    const std::vector<double>& acurvePoints, const std::vector<double>& bcurvePoints, const std::vector<double>& cccurvePoints,
                                    const std::vector<double>& lccurvePoints, LUTf & aoutCurve, LUTf & boutCurve, LUTf & satCurve, LUTf & lhskCurve,
                                    int skip)
{

    autili = butili = ccutili = cclutili = false;
    std::unique_ptr<DiagonalCurve> dCurve;

    // create a curve if needed
    if (!acurvePoints.empty() && acurvePoints[0] != 0) {
        dCurve = std::unique_ptr<DiagonalCurve>(new DiagonalCurve(acurvePoints, CURVES_MIN_POLY_POINTS / skip));

        if (dCurve && !dCurve->isIdentity()) {
            autili = true;
        }
    }

    CurveFactory::fillCurveArray(dCurve.get(), aoutCurve, skip, autili);

    dCurve = nullptr;

    //-----------------------------------------------------

    if (!bcurvePoints.empty() && bcurvePoints[0] != 0) {
        dCurve = std::unique_ptr<DiagonalCurve>(new DiagonalCurve(bcurvePoints, CURVES_MIN_POLY_POINTS / skip));

        if (dCurve && !dCurve->isIdentity()) {
            butili = true;
        }
    }

    CurveFactory::fillCurveArray(dCurve.get(), boutCurve, skip, butili);

    dCurve = nullptr;

    //-----------------------------------------------

    if (!cccurvePoints.empty() && cccurvePoints[0] != 0) {
        dCurve = std::unique_ptr<DiagonalCurve>(new DiagonalCurve(cccurvePoints, CURVES_MIN_POLY_POINTS / skip));

        if (dCurve && !dCurve->isIdentity()) {
            ccutili = true;
        }
    }

    CurveFactory::fillCurveArray(dCurve.get(), satCurve, skip, ccutili);

    dCurve = nullptr;

    //----------------------------

    if (!lccurvePoints.empty() && lccurvePoints[0] != 0) {
        dCurve = std::unique_ptr<DiagonalCurve>(new DiagonalCurve(lccurvePoints, CURVES_MIN_POLY_POINTS / skip));

        if (dCurve && !dCurve->isIdentity()) {
            cclutili = true;
        }
    }

    CurveFactory::fillCurveArray(dCurve.get(), lhskCurve, skip, cclutili);
}


} // namespace

void ImProcFunctions::labAdjustments(LabImage *lab, LUTu *histCCurve, LUTu *histLCurve)
{
    LUTu dummy;
    LUTu hist16;
    LUTf lumacurve;
    LUTf clcurve;
    LUTf satcurve;
    LUTf lhskcurve;
    LUTf curve1;
    LUTf curve2;

    hist16(65536);
    lumacurve(32770, 0); // lumacurve[32768] and lumacurve[32769] will be set to 32768 and 32769 later to allow linear interpolation
    clcurve(65536, 0);
    satcurve(65536, 0);
    lhskcurve(65536, 0);
    curve1(65536);
    curve2(65536);

    if (params->labCurve.contrast != 0) { //only use hist16 for contrast
        hist16.clear();
        int fh = lab->H;
        int fw = lab->W;
        
#ifdef _OPENMP
#       pragma omp parallel if (multiThread)
#endif
        {
            LUTu hist16thr (hist16.getSize());  // one temporary lookup table per thread
            hist16thr.clear();
#ifdef _OPENMP
#           pragma omp for schedule(static) nowait
#endif
            for (int i = 0; i < fh; i++) {
                for (int j = 0; j < fw; j++) {
                    hist16thr[(int)((lab->L[i][j]))]++;
                }
            }
            
#ifdef _OPENMP
#           pragma omp critical
#endif
            {
                hist16 += hist16thr;
            }
        }
    }  

    bool utili;
    complexLCurve(params->labCurve.brightness, params->labCurve.contrast, params->labCurve.lcurve, hist16, lumacurve, dummy, scale == 1 ? 1 : 16, utili);

    bool clcutili;
    curveCL(clcutili, params->labCurve.clcurve, clcurve, scale == 1 ? 1 : 16);

    bool autili, butili;
    bool ccutili, cclutili;
    complexsgnCurve(autili, butili, ccutili, cclutili, params->labCurve.acurve, params->labCurve.bcurve, params->labCurve.cccurve,
                                  params->labCurve.lccurve, curve1, curve2, satcurve, lhskcurve, scale == 1 ? 1 : 16);

    chromiLuminanceCurve(lab, lab, curve1, curve2, satcurve, lhskcurve, clcurve, lumacurve, utili, autili, butili, ccutili, cclutili, clcutili, histCCurve, histLCurve);
}


void ImProcFunctions::chromiLuminanceCurve(LabImage* lold, LabImage* lnew, LUTf & acurve, LUTf & bcurve, LUTf & satcurve, LUTf & lhskcurve, LUTf & clcurve, LUTf & curve, bool utili, bool autili, bool butili, bool ccutili, bool cclutili, bool clcutili, LUTu *histCCurve, LUTu *histLCurve)
{
    const auto get_hue_val =
        [](FlatCurve *c, float x) -> float
        {
            constexpr float mid = 0.5f;
            constexpr float range = 0.5f;
            constexpr float pivot = 0.5f;
            constexpr float base = 10.f;

            float y = c->getVal(Color::huelab_to_huehsv2(x));
            float yy = 0.f;
    
            if (y >= mid) {
                float v = (y - mid) / range;
                yy = pivot + (pow_F(base, v) - 1.f) / (base - 1.f) * pivot;
            } else {
                float v = (mid - y) / range;
                yy = pivot - (pow_F(base, v) - 1.f) / (base - 1.f) * pivot;
            }
            return yy;
        };

    int W = lold->W;
    int H = lold->H;

    PlanarWhateverData<float>* editWhatever = nullptr;
    EditUniqueID editID = EUID_None;
    bool editPipette = false;

    if (pipetteBuffer) {
        editID = pipetteBuffer->getEditID();

        if (editID != EUID_None) {

            switch  (pipetteBuffer->getDataProvider()->getCurrSubscriber()->getPipetteBufferType()) {
                case (BT_IMAGEFLOAT):
                    break;

                case (BT_LABIMAGE):
                    break;

                case (BT_SINGLEPLANE_FLOAT):
                    editPipette = true;
                    editWhatever = pipetteBuffer->getSinglePlaneBuffer();
                    break;
            }
        }
    }

    if (!params->labCurve.enabled) {
        if (editPipette && (editID == EUID_Lab_LCurve || editID == EUID_Lab_aCurve || editID == EUID_Lab_bCurve || editID == EUID_Lab_LHCurve || editID == EUID_Lab_CHCurve || editID == EUID_Lab_HHCurve || editID == EUID_Lab_CLCurve || editID == EUID_Lab_CCurve || editID == EUID_Lab_LCCurve)) {
            // fill pipette buffer with zeros to avoid crashes
            editWhatever->fill(0.f);
        }
        if (params->blackwhite.enabled) {
            for (int i = 0; i < lnew->H; ++i) {
                for (int j = 0; j < lnew->W; ++j) {
                    lnew->a[i][j] = lnew->b[i][j] = 0.f;
                }
            }
        }
        return;
    }

    // lhskcurve.dump("lh_curve");
    //init Flatcurve for C=f(H)


    FlatCurve* chCurve = nullptr;// curve C=f(H)
    bool chutili = false;

    if (params->labCurve.chromaticity > -100) {
        chCurve = new FlatCurve (params->labCurve.chcurve);

        if (!chCurve || chCurve->isIdentity()) {
            if (chCurve) {
                delete chCurve;
                chCurve = nullptr;
            }
        }//do not use "Munsell" if Chcurve not used
        else {
            chutili = true;
        }
    }

    FlatCurve* lhCurve = nullptr;//curve L=f(H)
    bool lhutili = false;

    if (params->labCurve.chromaticity > -100) {
        lhCurve = new FlatCurve (params->labCurve.lhcurve);

        if (!lhCurve || lhCurve->isIdentity()) {
            if (lhCurve) {
                delete lhCurve;
                lhCurve = nullptr;
            }
        }//do not use "Munsell" if Chcurve not used
        else {
            lhutili = true;
        }
    }

    FlatCurve* hhCurve = nullptr;//curve H=f(H)
    bool hhutili = false;

    if (params->labCurve.chromaticity > -100) {
        hhCurve = new FlatCurve (params->labCurve.hhcurve);

        if (!hhCurve || hhCurve->isIdentity()) {
            if (hhCurve) {
                delete hhCurve;
                hhCurve = nullptr;
            }
        }//do not use "Munsell" if Chcurve not used
        else {
            hhutili = true;
        }
    }


    const float histLFactor = histLCurve ? histLCurve->getSize() / 100.f : 1.f;
    const float histCFactor = histCCurve ? histCCurve->getSize() / 48000.f : 1.f;

#ifdef _DEBUG
    MyTime t1e, t2e;
    t1e.set();
    // init variables to display Munsell corrections
    MunsellDebugInfo* MunsDebugInfo = new MunsellDebugInfo();
#endif

    float adjustr = 1.0f;

//  if(params->labCurve.avoidclip ){
    // parameter to adapt curve C=f(C) to gamut

    // if      (params->icm.workingProfile == "ProPhoto")   {
    //     adjustr = 1.2f;   // 1.2 instead 1.0 because it's very rare to have C>170..
    // } else if (params->icm.workingProfile == "Adobe RGB")  {
    //     adjustr = 1.8f;
    // } else if (params->icm.workingProfile == "sRGB")       {
    //     adjustr = 2.0f;
    // } else if (params->icm.workingProfile == "WideGamut")  {
    //     adjustr = 1.2f;
    // } else if (params->icm.workingProfile == "Beta RGB")   {
    //     adjustr = 1.4f;
    // } else if (params->icm.workingProfile == "BestRGB")    {
    //     adjustr = 1.4f;
    // } else if (params->icm.workingProfile == "BruceRGB")   {
    //     adjustr = 1.8f;
    // }

    // reference to the params structure has to be done outside of the parallelization to avoid CPU cache problem
    const bool highlight = params->exposure.enabled && params->exposure.hrmode != procparams::ExposureParams::HR_OFF;
    const int chromaticity = params->labCurve.chromaticity;
    const float chromapro = (chromaticity + 100.0f) / 100.0f;
    const bool bwonly = params->blackwhite.enabled;
    bool bwq = false;
//  if(params->ppVersion > 300  && params->labCurve.chromaticity == - 100) bwq = true;
    // const bool bwToning = params->labCurve.chromaticity == - 100  /*|| params->blackwhite.method=="Ch" || params->blackwhite.enabled */ || bwonly;
    const bool bwToning = bwq  /*|| params->blackwhite.method=="Ch" || params->blackwhite.enabled */ || bwonly;
    //if(chromaticity==-100) chromaticity==-99;
    const bool LCredsk = params->labCurve.lcredsk;
    const bool ccut = ccutili;
    const bool clut = clcutili;
    const double rstprotection = 100. - params->labCurve.rstprotection; // Red and Skin Tones Protection
    // avoid color shift is disabled when bwToning is activated and enabled if gamut is true in colorappearanace
    const bool avoidColorShift = params->labCurve.avoidcolorshift && !bwToning ;
    const float protectRed = (float)settings->protectred;
    const double protectRedH = settings->protectredh;
    float protect_red, protect_redh;
    protect_red = protectRed;//default=60  chroma: one can put more or less if necessary...in 'option'  40...160

    if (protect_red < 20.0f) {
        protect_red = 20.0;    // avoid too low value
    }

    if (protect_red > 180.0f) {
        protect_red = 180.0;    // avoid too high value
    }

    protect_redh = float (protectRedH); //default=0.4 rad : one can put more or less if necessary...in 'option'  0.2 ..1.0

    if (protect_redh < 0.1f) {
        protect_redh = 0.1f;    //avoid divide by 0 and negatives values
    }

    if (protect_redh > 1.0f) {
        protect_redh = 1.0f;    //avoid too big values
    }

    float protect_redhcur = protectRedH;//default=0.4 rad : one can put more or less if necessary...in 'option'  0.2 ..1

    if (protect_redhcur < 0.1f) {
        protect_redhcur = 0.1f;    //avoid divide by 0 and negatives values:minimal protection for transition
    }

    if (protect_redhcur > 3.5f) {
        protect_redhcur = 3.5f;    //avoid too big values
    }

    //increase saturation after denoise : ...approximation
    float factnoise = 1.f;

    if (params->denoise.enabled) {
        factnoise = (1.f + params->denoise.chrominance / 500.f); //levels=5
//      if(yyyy) factnoise=(1.f+params->denoise.chrominance/100.f);//levels=7
    }

    const float scaleConst = 100.0f / 100.1f;


    const bool gamutLch = settings->gamutLch;
    const float amountchroma = (float) settings->amchroma;

    TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix (params->icm.workingProfile);
    const double wip[3][3] = {
        {wiprof[0][0], wiprof[0][1], wiprof[0][2]},
        {wiprof[1][0], wiprof[1][1], wiprof[1][2]},
        {wiprof[2][0], wiprof[2][1], wiprof[2][2]}
    };

    TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix (params->icm.workingProfile);
    const double wp[3][3] = {
        {wprof[0][0], wprof[0][1], wprof[0][2]},
        {wprof[1][0], wprof[1][1], wprof[1][2]},
        {wprof[2][0], wprof[2][1], wprof[2][2]}
    };

#ifdef _OPENMP
#ifdef _DEBUG
    #pragma omp parallel default(shared) firstprivate(lold, lnew, MunsDebugInfo) if (multiThread)
#else
    #pragma omp parallel if (multiThread)
#endif
#endif
    {
#ifdef __SSE2__
        float HHBuffer[W] ALIGNED16;
        float CCBuffer[W] ALIGNED16;
#endif
#ifdef _OPENMP
        #pragma omp for schedule(dynamic, 16)
#endif

        for (int i = 0; i < H; i++) {
            if (avoidColorShift) {
                bool hrenabled = params->exposure.enabled && params->exposure.hrmode != procparams::ExposureParams::HR_OFF;
                // only if user activate Lab adjustments
                if (autili || butili || ccutili ||  cclutili || chutili || lhutili || hhutili || clcutili || utili || chromaticity) {
                    Color::LabGamutMunsell (lold->L[i], lold->a[i], lold->b[i], W, /*corMunsell*/true, /*lumaMuns*/false, hrenabled, /*gamut*/true, wip);
                }
            }

#ifdef __SSE2__

            // precalculate some values using SSE
            if (bwToning || (!autili && !butili)) {
                __m128 c327d68v = _mm_set1_ps (327.68f);
                __m128 av, bv;
                int k;

                for (k = 0; k < W - 3; k += 4) {
                    av = LVFU (lold->a[i][k]);
                    bv = LVFU (lold->b[i][k]);
                    STVF (HHBuffer[k], xatan2f (bv, av));
                    STVF (CCBuffer[k], vsqrtf (SQRV (av) + SQRV (bv)) / c327d68v);
                }

                for (; k < W; k++) {
                    HHBuffer[k] = xatan2f (lold->b[i][k], lold->a[i][k]);
                    CCBuffer[k] = sqrt (SQR (lold->a[i][k]) + SQR (lold->b[i][k])) / 327.68f;
                }
            }

#endif // __SSE2__

            for (int j = 0; j < W; j++) {
                const float Lin = lold->L[i][j];
                float LL = Lin / 327.68f;
                float CC;
                float HH;
                float Chprov;
                float Chprov1;
                float memChprov;
                float2 sincosval;

                if (bwToning) { // this values will be also set when bwToning is false some lines down
#ifdef __SSE2__
                    // use precalculated values from above
                    HH = HHBuffer[j];
                    CC = CCBuffer[j];
#else
                    HH = xatan2f (lold->b[i][j], lold->a[i][j]);
                    CC = sqrt (SQR (lold->a[i][j]) + SQR (lold->b[i][j])) / 327.68f;
#endif

                    // According to mathematical laws we can get the sin and cos of HH by simple operations
                    if (CC == 0.0f) {
                        sincosval.y = 1.0f;
                        sincosval.x = 0.0f;
                    } else {
                        sincosval.y = lold->a[i][j] / (CC * 327.68f);
                        sincosval.x = lold->b[i][j] / (CC * 327.68f);
                    }

                    Chprov = CC;
                    Chprov1 = CC;
                    memChprov = Chprov;
                }

                if (editPipette && editID == EUID_Lab_LCurve) {
                    editWhatever->v (i, j) = LIM01<float> (Lin / 32768.0f);  // Lab L pipette
                }

                lnew->L[i][j] = curve[Lin];

                float Lprov1 = (lnew->L[i][j]) / 327.68f;

                if (editPipette) {
                    if (editID == EUID_Lab_aCurve) { // Lab a pipette
                        float chromapipa = lold->a[i][j] + (32768.f * 1.28f);
                        editWhatever->v (i, j) = LIM01<float> ((chromapipa) / (65536.f * 1.28f));
                    } else if (editID == EUID_Lab_bCurve) { //Lab b pipette
                        float chromapipb = lold->b[i][j] + (32768.f * 1.28f);
                        editWhatever->v (i, j) = LIM01<float> ((chromapipb) / (65536.f * 1.28f));
                    }
                }

                float atmp, btmp;

                atmp = lold->a[i][j];

                if (autili) {
                    atmp = acurve[atmp + 32768.0f] - 32768.0f;    // curves Lab a
                }

                btmp = lold->b[i][j];

                if (butili) {
                    btmp = bcurve[btmp + 32768.0f] - 32768.0f;    // curves Lab b
                }

                if (!bwToning) { //take into account modification of 'a' and 'b'
#ifdef __SSE2__
                    if (!autili && !butili) {
                        // use precalculated values from above
                        HH = HHBuffer[j];
                        CC = CCBuffer[j];
                    } else {
                        CC = sqrt (SQR (atmp) + SQR (btmp)) / 327.68f;
                        HH = xatan2f (btmp, atmp);
                    }

#else
                    CC = sqrt (SQR (atmp) + SQR (btmp)) / 327.68f;
                    HH = xatan2f (btmp, atmp);
#endif

                    // According to mathematical laws we can get the sin and cos of HH by simple operations
                    //float2  sincosval;
                    if (CC == 0.f) {
                        sincosval.y = 1.f;
                        sincosval.x = 0.f;
                    } else {
                        sincosval.y = atmp / (CC * 327.68f);
                        sincosval.x = btmp / (CC * 327.68f);
                    }

                    Chprov = CC;
                    Chprov1 = CC;
                    memChprov = Chprov;
                } // now new values of lold with 'a' and 'b'

                if (editPipette)
                    if (editID == EUID_Lab_LHCurve || editID == EUID_Lab_CHCurve || editID == EUID_Lab_HHCurve) {//H pipette
                        float valpar = Color::huelab_to_huehsv2 (HH);
                        editWhatever->v (i, j) = valpar;
                    }

                if (lhutili) {  // L=f(H)
                    const float ClipLevel = 65535.f;
                    float l_r;//Luminance Lab in 0..1
                    l_r = Lprov1 / 100.f;
                    {
                        float valparam = float ((get_hue_val(lhCurve, HH) - 0.5f)); //get l_r=f(H)
                        float valparamneg;
                        valparamneg = valparam;
                        float kcc = (CC / amountchroma); //take Chroma into account...40 "middle low" of chromaticity (arbitrary and simple), one can imagine other algorithme
                        //reduce action for low chroma and increase action for high chroma
                        valparam *= 2.f * kcc;
                        valparamneg *= kcc; //slightly different for negative

                        if (valparam > 0.f) {
                            l_r = (1.f - valparam) * l_r + valparam * (1.f - SQR (((SQR (1.f - min (l_r, 1.0f))))));
                        } else
                            //for negative
                        {
                            float khue = 1.9f; //in reserve in case of!
                            l_r *= (1.f + khue * valparamneg);
                        }
                    }

                    Lprov1 = l_r * 100.f;

                    float Chprov2 = sqrt (SQR (atmp) + SQR (btmp)) / 327.68f;
                    //Gamut control especially for negative values slightly different from gamutlchonly
                    bool inRGB;

                    do {
                        inRGB = true;
                        float aprov1 = Chprov2 * sincosval.y;
                        float bprov1 = Chprov2 * sincosval.x;

                        float fy = (Color::c1By116 * Lprov1 ) + Color::c16By116;
                        float fx = (0.002f * aprov1) + fy;
                        float fz = fy - (0.005f * bprov1);

                        float x_ = 65535.0f * Color::f2xyz (fx) * Color::D50x;
                        float z_ = 65535.0f * Color::f2xyz (fz) * Color::D50z;
                        float y_ = (Lprov1 > Color::epskap) ? 65535.0 * fy * fy * fy : 65535.0 * Lprov1 / Color::kappa;
                        float R, G, B;
                        Color::xyz2rgb (x_, y_, z_, R, G, B, wip);

                        if (R < 0.0f || G < 0.0f || B < 0.0f) {
                            if (Lprov1 < 0.1f) {
                                Lprov1 = 0.1f;
                            }

                            Chprov2 *= 0.95f;
                            inRGB = false;
                        } else if (!highlight && (R > ClipLevel || G > ClipLevel || B > ClipLevel)) {
                            if (Lprov1 > 99.98f) {
                                Lprov1 = 99.98f;
                            }

                            Chprov2 *= 0.95f;
                            inRGB = false;
                        }
                    } while (!inRGB);

                    atmp = 327.68f * Chprov2 * sincosval.y;
                    btmp = 327.68f * Chprov2 * sincosval.x;
                }

//          calculate C=f(H)
                if (chutili) {
                    float chparam = float ((get_hue_val(chCurve, HH) - 0.5f) * 2.0f); //get C=f(H)
                    float chromaChfactor = 1.0f + chparam;
                    atmp *= chromaChfactor;//apply C=f(H)
                    btmp *= chromaChfactor;
                }

                if (hhutili) {  // H=f(H)
                    //hue Lab in -PI +PI
                    float valparam = float ((get_hue_val(hhCurve, HH) - 0.5f) * 1.7f) + HH; //get H=f(H)  1.7 optimisation !
                    HH = valparam;
                    sincosval = xsincosf (HH);
                }

                if (!bwToning) {
                    float factorskin, factorsat, factorskinext;

                    if (chromapro > 1.f) {
                        float scale = scaleConst;//reduction in normal zone
                        float scaleext = 1.f;//reduction in transition zone
                        Color::scalered ( rstprotection, chromapro, 0.0, HH, protect_redh, scale, scaleext);//1.0
                        float interm = (chromapro - 1.f);
                        factorskin = 1.f + (interm * scale);
                        factorskinext = 1.f + (interm * scaleext);
                    } else {
                        factorskin = chromapro ; // +(chromapro)*scale;
                        factorskinext = chromapro ;// +(chromapro)*scaleext;
                    }

                    factorsat = chromapro * factnoise;

                    //simulate very approximative gamut f(L) : with pyramid transition
                    float dred /*=55.f*/;//C red value limit

                    if     (Lprov1 < 25.f) {
                        dred = 40.f;
                    } else if (Lprov1 < 30.f) {
                        dred = 3.f * Lprov1 - 35.f;
                    } else if (Lprov1 < 70.f) {
                        dred = 55.f;
                    } else if (Lprov1 < 75.f) {
                        dred = -3.f * Lprov1 + 265.f;
                    } else {
                        dred = 40.f;
                    }

                    // end pyramid

                    // Test if chroma is in the normal range first
                    Color::transitred ( HH, Chprov1, dred, factorskin, protect_red, factorskinext, protect_redh, factorsat, factorsat);
                    atmp *= factorsat;
                    btmp *= factorsat;

                    if (editPipette && editID == EUID_Lab_CLCurve) {
                        editWhatever->v (i, j) = LIM01<float> (LL / 100.f);  // Lab C=f(L) pipette
                    }

                    if (clut && LL > 0.f) { // begin C=f(L)
                        float factorskin, factorsat, factor, factorskinext;
                        float chromaCfactor = (clcurve[LL * 655.35f]) / (LL * 655.35f); //apply C=f(L)
                        float curf = 0.7f; //empirical coeff because curve is more progressive
                        float scale = 100.0f / 100.1f; //reduction in normal zone for curve C
                        float scaleext = 1.0f; //reduction in transition zone for curve C
                        float protect_redcur, protect_redhcur; //perhaps the same value than protect_red and protect_redh
                        float deltaHH;//HH value transition for C curve
                        protect_redcur = curf * protectRed; //default=60  chroma: one can put more or less if necessary...in 'option'  40...160==> curf =because curve is more progressive

                        if (protect_redcur < 20.0f) {
                            protect_redcur = 20.0;    // avoid too low value
                        }

                        if (protect_redcur > 180.0f) {
                            protect_redcur = 180.0;    // avoid too high value
                        }

                        protect_redhcur = curf * float (protectRedH); //default=0.4 rad : one can put more or less if necessary...in 'option'  0.2 ..1.0 ==> curf =because curve is more progressive

                        if (protect_redhcur < 0.1f) {
                            protect_redhcur = 0.1f;    //avoid divide by 0 and negatives values
                        }

                        if (protect_redhcur > 1.0f) {
                            protect_redhcur = 1.0f;    //avoid too big values
                        }

                        deltaHH = protect_redhcur; //transition hue

                        if (chromaCfactor > 0.0) {
                            Color::scalered ( rstprotection, chromaCfactor, 0.0, HH, deltaHH, scale, scaleext);    //1.0
                        }

                        if (chromaCfactor > 1.0) {
                            float interm = (chromaCfactor - 1.0f) * 100.0f;
                            factorskin = 1.0f + (interm * scale) / 100.0f;
                            factorskinext = 1.0f + (interm * scaleext) / 100.0f;
                        } else {
                            factorskin = chromaCfactor; // +(1.0f-chromaCfactor)*scale;
                            factorskinext = chromaCfactor ; //+(1.0f-chromaCfactor)*scaleext;
                        }

                        factorsat = chromaCfactor;
                        factor = factorsat;
                        Color::transitred ( HH, Chprov1, dred, factorskin, protect_redcur, factorskinext, deltaHH, factorsat, factor);
                        atmp = LIM(atmp * factor, min(-42000.f, atmp), max(42000.f, atmp));
                        btmp = LIM(btmp * factor, min(-42000.f, btmp), max(42000.f, btmp));
                    }

                    // end C=f(L)
                    //  if (editID == EUID_Lab_CLCurve)
                    //      editWhatever->v(i,j) = LIM01<float>(Lprov2/100.f);// Lab C=f(L) pipette

                    // I have placed C=f(C) after all C treatments to assure maximum amplitude of "C"
                    if (editPipette && editID == EUID_Lab_CCurve) {
                        float chromapip = sqrt (SQR (atmp) + SQR (btmp) + 0.001f);
                        editWhatever->v (i, j) = LIM01<float> ((chromapip) / (48000.f));
                    }//Lab C=f(C) pipette

                    if (ccut) {
                        float factorskin, factorsat, factor, factorskinext;
                        float chroma = sqrt (SQR (atmp) + SQR (btmp) + 0.001f);
                        float chromaCfactor = (satcurve[chroma * adjustr]) / (chroma * adjustr); //apply C=f(C)
                        float curf = 0.7f; //empirical coeff because curve is more progressive
                        float scale = 100.0f / 100.1f; //reduction in normal zone for curve CC
                        float scaleext = 1.0f; //reduction in transition zone for curve CC
                        float protect_redcur, protect_redhcur; //perhaps the same value than protect_red and protect_redh
                        float deltaHH;//HH value transition for CC curve
                        protect_redcur = curf * protectRed; //default=60  chroma: one can put more or less if necessary...in 'option'  40...160==> curf =because curve is more progressive

                        if (protect_redcur < 20.0f) {
                            protect_redcur = 20.0;    // avoid too low value
                        }

                        if (protect_redcur > 180.0f) {
                            protect_redcur = 180.0;    // avoid too high value
                        }

                        protect_redhcur = curf * float (protectRedH); //default=0.4 rad : one can put more or less if necessary...in 'option'  0.2 ..1.0 ==> curf =because curve is more progressive

                        if (protect_redhcur < 0.1f) {
                            protect_redhcur = 0.1f;    //avoid divide by 0 and negatives values
                        }

                        if (protect_redhcur > 1.0f) {
                            protect_redhcur = 1.0f;    //avoid too big values
                        }

                        deltaHH = protect_redhcur; //transition hue

                        if (chromaCfactor > 0.0) {
                            Color::scalered ( rstprotection, chromaCfactor, 0.0, HH, deltaHH, scale, scaleext);    //1.0
                        }

                        if (chromaCfactor > 1.0) {
                            float interm = (chromaCfactor - 1.0f) * 100.0f;
                            factorskin = 1.0f + (interm * scale) / 100.0f;
                            factorskinext = 1.0f + (interm * scaleext) / 100.0f;
                        } else {
                            //factorskin= chromaCfactor*scale;
                            //factorskinext=chromaCfactor*scaleext;
                            factorskin = chromaCfactor; // +(1.0f-chromaCfactor)*scale;
                            factorskinext = chromaCfactor ; //+(1.0f-chromaCfactor)*scaleext;

                        }

                        factorsat = chromaCfactor;
                        factor = factorsat;
                        Color::transitred ( HH, Chprov1, dred, factorskin, protect_redcur, factorskinext, deltaHH, factorsat, factor);
                        atmp *= factor;
                        btmp *= factor;
                    }
                }

                // end chroma C=f(C)

                //update histogram C
                if (histCCurve) { //only with improccoordinator
                    (*histCCurve)[histCFactor * sqrt(atmp * atmp + btmp * btmp)]++;
                }

                if (editPipette && editID == EUID_Lab_LCCurve) {
                    float chromapiplc = sqrt (SQR (atmp) + SQR (btmp) + 0.001f);
                    editWhatever->v (i, j) = LIM01<float> ((chromapiplc) / (48000.f));
                }//Lab L=f(C) pipette


                if (cclutili && !bwToning) {    //apply curve L=f(C) for skin and rd...but also for extended color ==> near green and blue (see 'curf')

                    const float xx = 0.25f; //soft : between 0.2 and 0.4
                    float skdeltaHH;

                    skdeltaHH = protect_redhcur; //transition hue

                    float skbeg = -0.05f; //begin hue skin
                    float skend = 1.60f; //end hue skin
                    const float chrmin = 50.0f; //to avoid artifact, because L curve is not a real curve for luminance
                    float aa, bb;
                    float zz = 0.0f;
                    float yy = 0.0f;

                    if (Chprov1 < chrmin) {
                        yy = SQR (Chprov1 / chrmin) * xx;
                    } else {
                        yy = xx;    //avoid artifact for low C
                    }

                    if (!LCredsk) {
                        skbeg = -3.1415;
                        skend = 3.14159;
                        skdeltaHH = 0.001f;
                    }

                    if (HH > skbeg && HH < skend ) {
                        zz = yy;
                    } else if (HH > skbeg - skdeltaHH && HH <= skbeg) { //transition
                        aa = yy / skdeltaHH;
                        bb = -aa * (skbeg - skdeltaHH);
                        zz = aa * HH + bb;
                    } else if (HH >= skend && HH < skend + skdeltaHH) { //transition
                        aa = -yy / skdeltaHH;
                        bb = -aa * (skend + skdeltaHH);
                        zz = aa * HH + bb;
                    }

                    float chroma = sqrt (SQR (atmp) + SQR (btmp) + 0.001f);
                    float Lc = (lhskcurve[chroma * adjustr]) / (chroma * adjustr); //apply L=f(C)
                    Lc = (Lc - 1.0f) * zz + 1.0f; //reduct action
                    Lprov1 *= Lc; //adjust luminance
                }

                //update histo LC
                if (histLCurve) { //only with improccoordinator
                    (*histLCurve)[Lprov1 * histLFactor]++;
                }

                Chprov1 = sqrt (SQR (atmp) + SQR (btmp)) / 327.68f;

                // labCurve.bwtoning option allows to decouple modulation of a & b curves by saturation
                // with bwtoning enabled the net effect of a & b curves is visible
                if (bwToning) {
                    atmp -= lold->a[i][j];
                    btmp -= lold->b[i][j];
                }

                if (avoidColorShift) {
                    //gamutmap Lch ==> preserve Hue,but a little slower than gamutbdy for high values...and little faster for low values
                    if (gamutLch) {
                        float R, G, B;

#ifdef _DEBUG
                        bool neg = false;
                        bool more_rgb = false;
                        //gamut control : Lab values are in gamut
                        Color::gamutLchonly (HH, sincosval, Lprov1, Chprov1, R, G, B, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
                        //gamut control : Lab values are in gamut
                        Color::gamutLchonly (HH, sincosval, Lprov1, Chprov1, R, G, B, wip, highlight, 0.15f, 0.96f);
#endif
                        lnew->L[i][j] = Lprov1 * 327.68f;
//                  float2 sincosval = xsincosf(HH);
                        lnew->a[i][j] = 327.68f * Chprov1 * sincosval.y;
                        lnew->b[i][j] = 327.68f * Chprov1 * sincosval.x;
                    } else {
                        //use gamutbdy
                        //Luv limiter
                        float Y, u, v;
                        Color::Lab2Yuv (lnew->L[i][j], atmp, btmp, Y, u, v);
                        //Yuv2Lab includes gamut restriction map
                        Color::Yuv2Lab (Y, u, v, lnew->L[i][j], lnew->a[i][j], lnew->b[i][j], wp);
                    }

                    if (utili || autili || butili || ccut || clut || cclutili || chutili || lhutili || hhutili || clcutili || chromaticity) {
                        float correctionHue = 0.f; // Munsell's correction
                        float correctlum = 0.f;

                        Lprov1 = lnew->L[i][j] / 327.68f;
                        Chprov = sqrt (SQR (lnew->a[i][j]) + SQR (lnew->b[i][j])) / 327.68f;

#ifdef _DEBUG
                        Color::AllMunsellLch (/*lumaMuns*/true, Lprov1, LL, HH, Chprov, memChprov, correctionHue, correctlum, MunsDebugInfo);
#else
                        Color::AllMunsellLch (/*lumaMuns*/true, Lprov1, LL, HH, Chprov, memChprov, correctionHue, correctlum);
#endif

                        if (correctionHue != 0.f || correctlum != 0.f) {
                            if (fabs (correctionHue) < 0.015f) {
                                HH += correctlum;    // correct only if correct Munsell chroma very little.
                            }

                            /*      if((HH>0.0f && HH < 1.6f)   && memChprov < 70.0f) HH+=correctlum;//skin correct
                                    else if(fabs(correctionHue) < 0.3f) HH+=0.08f*correctlum;
                                    else if(fabs(correctionHue) < 0.2f) HH+=0.25f*correctlum;
                                    else if(fabs(correctionHue) < 0.1f) HH+=0.35f*correctlum;
                                    else if(fabs(correctionHue) < 0.015f) HH+=correctlum;   // correct only if correct Munsell chroma very little.
                            */
                            sincosval = xsincosf (HH + correctionHue);
                        }

                        lnew->a[i][j] = 327.68f * Chprov * sincosval.y; // apply Munsell
                        lnew->b[i][j] = 327.68f * Chprov * sincosval.x;
                    }
                } else {
//              if(Lprov1 > maxlp) maxlp=Lprov1;
//              if(Lprov1 < minlp) minlp=Lprov1;
                    if (!bwToning) {
                        lnew->L[i][j] = Lprov1 * 327.68f;
//                  float2 sincosval = xsincosf(HH);
                        lnew->a[i][j] = 327.68f * Chprov1 * sincosval.y;
                        lnew->b[i][j] = 327.68f * Chprov1 * sincosval.x;
                    } else {
                        //Luv limiter only
                        lnew->a[i][j] = atmp;
                        lnew->b[i][j] = btmp;
                    }
                }
            }
        }
    } // end of parallelization

#ifdef _DEBUG

    if (settings->verbose) {
        t2e.set();
        printf ("Color::AllMunsellLch (correction performed in %d usec):\n", t2e.etime (t1e));
        printf ("   Munsell chrominance: MaxBP=%1.2frad MaxRY=%1.2frad MaxGY=%1.2frad MaxRP=%1.2frad  dep=%u\n", MunsDebugInfo->maxdhue[0],    MunsDebugInfo->maxdhue[1],    MunsDebugInfo->maxdhue[2],    MunsDebugInfo->maxdhue[3],    MunsDebugInfo->depass);
        printf ("   Munsell luminance  : MaxBP=%1.2frad MaxRY=%1.2frad MaxGY=%1.2frad MaxRP=%1.2frad  dep=%u\n", MunsDebugInfo->maxdhuelum[0], MunsDebugInfo->maxdhuelum[1], MunsDebugInfo->maxdhuelum[2], MunsDebugInfo->maxdhuelum[3], MunsDebugInfo->depassLum);
    }

    delete MunsDebugInfo;
#endif

    if (chCurve) {
        delete chCurve;
    }

    if (lhCurve) {
        delete lhCurve;
    }

    if (hhCurve) {
        delete hhCurve;
    }

    //  t2e.set();
    //  printf("Chromil took %d nsec\n",t2e.etime(t1e));
}

} // namespace rtengine
