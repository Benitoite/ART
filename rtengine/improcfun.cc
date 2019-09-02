/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
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
#include <cmath>
#include <glib.h>
#include <glibmm.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "alignedbuffer.h"
#include "rtengine.h"
#include "improcfun.h"
#include "curves.h"
#include "mytime.h"
#include "iccstore.h"
#include "imagesource.h"
#include "rtthumbnail.h"
#include "utils.h"
#include "iccmatrices.h"
#include "color.h"
#include "calc_distort.h"
#include "rt_math.h"
#include "EdgePreservingDecomposition.h"
#include "improccoordinator.h"
#include "clutstore.h"
#include "StopWatch.h"
#include "../rtgui/ppversion.h"
#include "../rtgui/guiutils.h"

namespace rtengine {

using namespace procparams;

extern const Settings* settings;

ImProcFunctions::ImProcFunctions(const ProcParams* iparams, bool imultiThread):
    monitorTransform(nullptr),
    params(iparams),
    scale(1),
    multiThread(imultiThread),
    dcpProf(nullptr),
    dcpApplyState(nullptr),
    pipetteBuffer(nullptr),
    lumimul{},
    offset_x(0),
    offset_y(0),
    full_width(-1),
    full_height(-1),
    hist_tonecurve(nullptr),
    hist_ccurve(nullptr),
    hist_lcurve(nullptr),
    show_sharpening_mask(false)
{
}


ImProcFunctions::~ImProcFunctions ()
{
    if (monitorTransform) {
        cmsDeleteTransform (monitorTransform);
    }
}

void ImProcFunctions::setScale (double iscale)
{
    scale = iscale;
}


void ImProcFunctions::updateColorProfiles (const Glib::ustring& monitorProfile, RenderingIntent monitorIntent, bool softProof, bool gamutCheck)
{
    // set up monitor transform
    if (monitorTransform) {
        cmsDeleteTransform (monitorTransform);
    }
    gamutWarning.reset(nullptr);

    monitorTransform = nullptr;

    cmsHPROFILE monitor = nullptr;

    if (!monitorProfile.empty()) {
#if !defined(__APPLE__) // No support for monitor profiles on OS X, all data is sRGB
        monitor = ICCStore::getInstance()->getProfile (monitorProfile);
#else
        monitor = ICCStore::getInstance()->getProfile (options.rtSettings.srgb);
#endif
    }

    if (monitor) {
        MyMutex::MyLock lcmsLock (*lcmsMutex);

        cmsUInt32Number flags;
        cmsHPROFILE iprof  = cmsCreateLab4Profile (nullptr);
        cmsHPROFILE gamutprof = nullptr;
        cmsUInt32Number gamutbpc = 0;
        RenderingIntent gamutintent = RI_RELATIVE;

        bool softProofCreated = false;

        if (softProof) {
            cmsHPROFILE oprof = nullptr;
            RenderingIntent outIntent;
            
            flags = cmsFLAGS_SOFTPROOFING | cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE;

            if (!settings->printerProfile.empty()) {
                oprof = ICCStore::getInstance()->getProfile (settings->printerProfile);
                if (settings->printerBPC) {
                    flags |= cmsFLAGS_BLACKPOINTCOMPENSATION;
                }
                outIntent = settings->printerIntent;
            } else {
                oprof = ICCStore::getInstance()->getProfile(params->icm.outputProfile);
                if (params->icm.outputBPC) {
                    flags |= cmsFLAGS_BLACKPOINTCOMPENSATION;
                }
                outIntent = params->icm.outputIntent;
            }

            if (oprof) {
                // NOCACHE is for thread safety, NOOPTIMIZE for precision

                // if (gamutCheck) {
                //     flags |= cmsFLAGS_GAMUTCHECK;
                // }

                const auto make_gamma_table =
                    [](cmsHPROFILE prof, cmsTagSignature tag) -> void
                    {
                        cmsToneCurve *tc = static_cast<cmsToneCurve *>(cmsReadTag(prof, tag));
                        if (tc) {
                            const cmsUInt16Number *table = cmsGetToneCurveEstimatedTable(tc);
                            cmsToneCurve *tc16 = cmsBuildTabulatedToneCurve16(nullptr, cmsGetToneCurveEstimatedTableEntries(tc), table);
                            if (tc16) {
                                cmsWriteTag(prof, tag, tc16);
                                cmsFreeToneCurve(tc16);
                            }
                        }
                    };

                cmsHPROFILE softproof = ProfileContent(oprof).toProfile();
                if (softproof) {
                    make_gamma_table(softproof, cmsSigRedTRCTag);
                    make_gamma_table(softproof, cmsSigGreenTRCTag);
                    make_gamma_table(softproof, cmsSigBlueTRCTag);
                }

                monitorTransform = cmsCreateProofingTransform (
                                       iprof, TYPE_Lab_FLT,
                                       monitor, TYPE_RGB_FLT,
                                       softproof, //oprof,
                                       monitorIntent, outIntent,
                                       flags
                                   );

                if (softproof) {
                    cmsCloseProfile(softproof);
                }

                if (monitorTransform) {
                    softProofCreated = true;
                }

                if (gamutCheck) {
                    gamutprof = oprof;
                    if (params->icm.outputBPC) {
                        gamutbpc = cmsFLAGS_BLACKPOINTCOMPENSATION;
                    }
                    gamutintent = outIntent;
                }
            }
        } else if (gamutCheck) {
            // flags = cmsFLAGS_GAMUTCHECK | cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE;
            // if (settings->monitorBPC) {
            //     flags |= cmsFLAGS_BLACKPOINTCOMPENSATION;
            // }

            // monitorTransform = cmsCreateProofingTransform(iprof, TYPE_Lab_FLT, monitor, TYPE_RGB_8, monitor, monitorIntent, monitorIntent, flags);

            // if (monitorTransform) {
            //     softProofCreated = true;
            // }
            gamutprof = monitor;
            if (settings->monitorBPC) {
                gamutbpc = cmsFLAGS_BLACKPOINTCOMPENSATION;
            }
            gamutintent = monitorIntent;
        }

        if (!softProofCreated) {
            flags = cmsFLAGS_NOOPTIMIZE | cmsFLAGS_NOCACHE;

            if (settings->monitorBPC) {
                flags |= cmsFLAGS_BLACKPOINTCOMPENSATION;
            }

            monitorTransform = cmsCreateTransform (iprof, TYPE_Lab_FLT, monitor, TYPE_RGB_FLT, monitorIntent, flags);
        }

        if (gamutCheck && gamutprof) {
            gamutWarning.reset(new GamutWarning(iprof, gamutprof, gamutintent, gamutbpc));
        }

        cmsCloseProfile (iprof);
    }
}

void ImProcFunctions::firstAnalysis (const Imagefloat* const original, const ProcParams &params, LUTu & histogram)
{

    TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix (params.icm.workingProfile);

    lumimul[0] = wprof[1][0];
    lumimul[1] = wprof[1][1];
    lumimul[2] = wprof[1][2];
    int W = original->getWidth();
    int H = original->getHeight();

    float lumimulf[3] = {static_cast<float> (lumimul[0]), static_cast<float> (lumimul[1]), static_cast<float> (lumimul[2])};

    // calculate histogram of the y channel needed for contrast curve calculation in exposure adjustments
    histogram.clear();

    if (multiThread) {

#ifdef _OPENMP
        const int numThreads = min (max (W * H / (int)histogram.getSize(), 1), omp_get_max_threads());
        #pragma omp parallel num_threads(numThreads) if(numThreads>1)
#endif
        {
            LUTu hist (histogram.getSize());
            hist.clear();
#ifdef _OPENMP
            #pragma omp for nowait
#endif

            for (int i = 0; i < H; i++) {
                for (int j = 0; j < W; j++) {

                    float r = original->r (i, j);
                    float g = original->g (i, j);
                    float b = original->b (i, j);

                    int y = (lumimulf[0] * r + lumimulf[1] * g + lumimulf[2] * b);
                    hist[y]++;
                }
            }

#ifdef _OPENMP
            #pragma omp critical
#endif
            histogram += hist;

        }
#ifdef _OPENMP
        static_cast<void> (numThreads); // to silence cppcheck warning
#endif
    } else {
        for (int i = 0; i < H; i++) {
            for (int j = 0; j < W; j++) {

                float r = original->r (i, j);
                float g = original->g (i, j);
                float b = original->b (i, j);

                int y = (lumimulf[0] * r + lumimulf[1] * g + lumimulf[2] * b);
                histogram[y]++;
            }
        }
    }
}


#if 0

// Process RGB image and convert to LAB space
void ImProcFunctions::rgbProc(Imagefloat *working)
{
    BENCHFUN

    working->setMode(Imagefloat::Mode::RGB, multiThread);
        
    constexpr int TS = 112;
    
    LUTf hltonecurve(65536);
    LUTf shtonecurve(65536);
    LUTf rCurve, gCurve, bCurve;

    double expcomp = params->exposure.enabled ? params->exposure.expcomp : 0.0;
    int hlcompr = params->exposure.enabled ? params->exposure.hlcompr : 0;
    int hlcomprthresh = params->exposure.hlcomprthresh;

    {
        int black = params->exposure.enabled ? params->exposure.black : 0;
        int shcompr = params->exposure.enabled ? params->exposure.shcompr : 0;
        
        LUTf tonecurve(65536);
        LUTu vhist16(65536), histToneCurve(256);
        ToneCurve customToneCurve1, customToneCurve2;
        
        CurveFactory::complexCurve(expcomp, black / 65535.0,
                                   hlcompr, hlcomprthresh,
                                   shcompr, 0, 0, 
                                   { DCT_Linear }, { DCT_Linear },
                                   vhist16, hltonecurve, shtonecurve, tonecurve,
                                   histToneCurve, customToneCurve1,
                                   customToneCurve2, scale);
    }

    // std::unique_ptr<Imagefloat> workimage(working->copy());
    // working = workimage.get();
    
    Imagefloat *tmpImage = nullptr;

    Imagefloat* editImgFloat = nullptr;
    PlanarWhateverData<float>* editWhatever = nullptr;
    EditUniqueID editID = pipetteBuffer ? pipetteBuffer->getEditID() : EUID_None;

    if (editID != EUID_None) {
        switch  (pipetteBuffer->getDataProvider()->getCurrSubscriber()->getPipetteBufferType()) {
            case (BT_IMAGEFLOAT):
                editImgFloat = pipetteBuffer->getImgFloatBuffer();
                break;

            case (BT_LABIMAGE):
                break;

            case (BT_SINGLEPLANE_FLOAT):
                editWhatever = pipetteBuffer->getSinglePlaneBuffer();
                break;
        }
    }

    if (params->rgbCurves.enabled) {
        CurveFactory::RGBCurve(params->rgbCurves.rcurve, rCurve, scale);
        CurveFactory::RGBCurve(params->rgbCurves.gcurve, gCurve, scale);
        CurveFactory::RGBCurve(params->rgbCurves.bcurve, bCurve, scale);
    }

    TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix (params->icm.workingProfile);
    TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix (params->icm.workingProfile);

    float toxyz[3][3] = {
        {
            static_cast<float> ( wprof[0][0] / Color::D50x),
            static_cast<float> ( wprof[0][1] / Color::D50x),
            static_cast<float> ( wprof[0][2] / Color::D50x)
        }, {
            static_cast<float> ( wprof[1][0]),
            static_cast<float> ( wprof[1][1]),
            static_cast<float> ( wprof[1][2])
        }, {
            static_cast<float> ( wprof[2][0] / Color::D50z),
            static_cast<float> ( wprof[2][1] / Color::D50z),
            static_cast<float> ( wprof[2][2] / Color::D50z)
        }
    };
    float maxFactorToxyz = max (toxyz[1][0], toxyz[1][1], toxyz[1][2]);
    float equalR = maxFactorToxyz / toxyz[1][0];
    float equalG = maxFactorToxyz / toxyz[1][1];
    float equalB = maxFactorToxyz / toxyz[1][2];

    //inverse matrix user select
    double wip[3][3] = {
        {wiprof[0][0], wiprof[0][1], wiprof[0][2]},
        {wiprof[1][0], wiprof[1][1], wiprof[1][2]},
        {wiprof[2][0], wiprof[2][1], wiprof[2][2]}
    };

    double wp[3][3] = {
        {wprof[0][0], wprof[0][1], wprof[0][2]},
        {wprof[1][0], wprof[1][1], wprof[1][2]},
        {wprof[2][0], wprof[2][1], wprof[2][2]}
    };

    bool mixchannels = params->chmixer.enabled &&
        (params->chmixer.red[0] != 100 || params->chmixer.red[1] != 0     || params->chmixer.red[2] != 0   ||
                        params->chmixer.green[0] != 0 || params->chmixer.green[1] != 100 || params->chmixer.green[2] != 0 ||
                        params->chmixer.blue[0] != 0  || params->chmixer.blue[1] != 0    || params->chmixer.blue[2] != 100);

    FlatCurve* bwlCurve = nullptr;

    FlatCurveType bwlCurveType = (FlatCurveType)params->blackwhite.luminanceCurve.at (0);
    bool bwlCurveEnabled = bwlCurveType > FCT_Linear;

    if (bwlCurveEnabled) {
        bwlCurve = new FlatCurve (params->blackwhite.luminanceCurve);

        if (bwlCurve->isIdentity()) {
            delete bwlCurve;
            bwlCurve = nullptr;
            bwlCurveEnabled = false;
        }
    }

    const float exp_scale = pow (2.0, expcomp);
    const float comp = (max (0.0, expcomp) + 1.0) * hlcompr / 100.0;
    const float shoulder = ((65536.0 / max (1.0f, exp_scale)) * (hlcomprthresh / 200.0)) + 0.1;
    const float hlrange = 65536.0 - shoulder;
    const bool isProPhoto = (params->icm.workingProfile == "ProPhoto");
    bool highlight = params->exposure.enabled && params->exposure.hrmode != procparams::ExposureParams::HR_OFF;

    float chMixRR = float (params->chmixer.red[0])/10.f;
    float chMixRG = float (params->chmixer.red[1])/10.f;
    float chMixRB = float (params->chmixer.red[2])/10.f;
    float chMixGR = float (params->chmixer.green[0])/10.f;
    float chMixGG = float (params->chmixer.green[1])/10.f;
    float chMixGB = float (params->chmixer.green[2])/10.f;
    float chMixBR = float (params->chmixer.blue[0])/10.f;
    float chMixBG = float (params->chmixer.blue[1])/10.f;
    float chMixBB = float (params->chmixer.blue[2])/10.f;

    bool blackwhite = params->blackwhite.enabled;
    float bwr = float (params->blackwhite.mixerRed);
    float bwg = float (params->blackwhite.mixerGreen);
    float bwb = float (params->blackwhite.mixerBlue);
    float bwrgam = float (params->blackwhite.gammaRed);
    float bwggam = float (params->blackwhite.gammaGreen);
    float bwbgam = float (params->blackwhite.gammaBlue);
    int algm = 0;

    if     (params->blackwhite.method == "Desaturation") {
        algm = 0;
    } else if (params->blackwhite.method == "LumEqualizer") {
        algm = 1;
    } else if (params->blackwhite.method == "ChannelMixer") {
        algm = 2;
    }

    float kcorec = 1.f;
    //gamma correction of each channel
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
    //if     (params->blackwhite.setting=="Ma" || params->blackwhite.setting=="Mr" || params->blackwhite.setting=="Fr" || params->blackwhite.setting=="Fa")  {
    {
        gammabwr = 1.f - bwrgam / gamvalr;
        gammabwg = 1.f - bwggam / gamvalg;
        gammabwb = 1.f - bwbgam / gamvalb;
    }
    bool hasgammabw = gammabwr != 1.f || gammabwg != 1.f || gammabwb != 1.f;

    if (blackwhite) {// || (params->dirpyrequalizer.cbdlMethod == "bef" && params->dirpyrequalizer.enabled)) {
        tmpImage = new Imagefloat (working->getWidth(), working->getHeight());
    }

    if (mixchannels) {
#ifdef _OPENMP
#       pragma omp parallel for if (multiThread)
#endif
        for (int y = 0; y < working->getHeight(); ++y) {
            for (int x = 0; x < working->getWidth(); ++x) {
                float r = working->r(y, x);
                float g = working->g(y, x);
                float b = working->b(y, x);

                float rmix = (r * chMixRR + g * chMixRG + b * chMixRB) / 100.f;
                float gmix = (r * chMixGR + g * chMixGG + b * chMixGB) / 100.f;
                float bmix = (r * chMixBR + g * chMixBG + b * chMixBB) / 100.f;

                working->r(y, x) = rmix;
                working->g(y, x) = gmix;
                working->b(y, x) = bmix;
            }
        }
    }

    highlightToneCurve(hltonecurve, working, exp_scale, comp, hlrange, multiThread);
    shadowToneCurve(shtonecurve, working, multiThread);
    toneEqualizer(working);

#ifdef _OPENMP
    #pragma omp parallel if (multiThread)
#endif
    {
        size_t perChannelSizeBytes = padToAlignment(sizeof (float) * TS * TS + 4 * 64);
        AlignedBuffer<float> buffer(3 * perChannelSizeBytes);
        char *editIFloatBuffer = nullptr;
        char *editWhateverBuffer = nullptr;

        float *rtemp = buffer.data;
        float *gtemp = &rtemp[perChannelSizeBytes / sizeof(float)];
        float *btemp = &gtemp[perChannelSizeBytes / sizeof(float)];
        int istart;
        int jstart;
        int tW;
        int tH;

        // zero out the buffers
        memset(rtemp, 0, 3 * perChannelSizeBytes);

        float *editIFloatTmpR = nullptr, *editIFloatTmpG = nullptr, *editIFloatTmpB = nullptr, *editWhateverTmp = nullptr;

        if (editImgFloat) {
            editIFloatBuffer = (char *) malloc (3 * sizeof (float) * TS * TS + 20 * 64 + 63);
            char *data = (char*) ( ( uintptr_t (editIFloatBuffer) + uintptr_t (63)) / 64 * 64);

            editIFloatTmpR = (float (*))data;
            editIFloatTmpG = (float (*))         ((char*)editIFloatTmpR + sizeof (float) * TS * TS + 4 * 64);
            editIFloatTmpB = (float (*))         ((char*)editIFloatTmpG + sizeof (float) * TS * TS + 8 * 64);
        }

        if (editWhatever) {
            editWhateverBuffer = (char *) malloc (sizeof (float) * TS * TS + 20 * 64 + 63);
            char *data = (char*) ( ( uintptr_t (editWhateverBuffer) + uintptr_t (63)) / 64 * 64);

            editWhateverTmp = (float (*))data;
        }

        // float out_rgbx[4 * TS] ALIGNED16; // Line buffer for CLUT
        // float clutr[TS] ALIGNED16;
        // float clutg[TS] ALIGNED16;
        // float clutb[TS] ALIGNED16;

#ifdef _OPENMP
        #pragma omp for schedule(dynamic) collapse(2)
#endif
        for (int ii = 0; ii < working->getHeight(); ii += TS) {
            for (int jj = 0; jj < working->getWidth(); jj += TS) {
                istart = ii;
                jstart = jj;
                tH = min (ii + TS, working->getHeight());
                tW = min (jj + TS, working->getWidth());


                for (int i = istart, ti = 0; i < tH; i++, ti++) {
                    for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                        rtemp[ti * TS + tj] = working->r (i, j);
                        gtemp[ti * TS + tj] = working->g (i, j);
                        btemp[ti * TS + tj] = working->b (i, j);
                    }
                }
                
                if (editID == EUID_RGB_R) {
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            editWhateverTmp[ti * TS + tj] = Color::gamma2curve[rtemp[ti * TS + tj]] / 65536.f;
                        }
                    }
                } else if (editID == EUID_RGB_G) {
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            editWhateverTmp[ti * TS + tj] = Color::gamma2curve[gtemp[ti * TS + tj]] / 65536.f;
                        }
                    }
                } else if (editID == EUID_RGB_B) {
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            editWhateverTmp[ti * TS + tj] = Color::gamma2curve[btemp[ti * TS + tj]] / 65536.f;
                        }
                    }
                }

                if (params->rgbCurves.enabled && (rCurve || gCurve || bCurve)) { // if any of the RGB curves is engaged
                    if (!params->rgbCurves.lumamode) { // normal RGB mode

                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                // individual R tone curve
                                if (rCurve) {
                                    setUnlessOOG(rtemp[ti * TS + tj], rCurve[ rtemp[ti * TS + tj] ]);
                                }

                                // individual G tone curve
                                if (gCurve) {
                                    setUnlessOOG(gtemp[ti * TS + tj], gCurve[ gtemp[ti * TS + tj] ]);
                                }

                                // individual B tone curve
                                if (bCurve) {
                                    setUnlessOOG(btemp[ti * TS + tj], bCurve[ btemp[ti * TS + tj] ]);
                                }
                            }
                        }
                    } else { //params->rgbCurves.lumamode==true (Luminosity mode)
                        // rCurve.dump("r_curve");//debug

                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                // rgb values before RGB curves
                                float r = rtemp[ti * TS + tj] ;
                                float g = gtemp[ti * TS + tj] ;
                                float b = btemp[ti * TS + tj] ;
                                //convert to Lab to get a&b before RGB curves
                                float x = toxyz[0][0] * r + toxyz[0][1] * g + toxyz[0][2] * b;
                                float y = toxyz[1][0] * r + toxyz[1][1] * g + toxyz[1][2] * b;
                                float z = toxyz[2][0] * r + toxyz[2][1] * g + toxyz[2][2] * b;

                                float fx = x < MAXVALF ? Color::cachef[x] : 327.68f * std::cbrt (x / MAXVALF);
                                float fy = y < MAXVALF ? Color::cachef[y] : 327.68f * std::cbrt (y / MAXVALF);
                                float fz = z < MAXVALF ? Color::cachef[z] : 327.68f * std::cbrt (z / MAXVALF);

                                float a_1 = 500.0f * (fx - fy);
                                float b_1 = 200.0f * (fy - fz);

                                // rgb values after RGB curves
                                if (rCurve) {
                                    float rNew = rCurve[r];
                                    r += (rNew - r) * equalR;
                                }

                                if (gCurve) {
                                    float gNew = gCurve[g];
                                    g += (gNew - g) * equalG;
                                }

                                if (bCurve) {
                                    float bNew = bCurve[b];
                                    b += (bNew - b) * equalB;
                                }

                                // Luminosity after
                                // only Luminance in Lab
                                float newy = toxyz[1][0] * r + toxyz[1][1] * g + toxyz[1][2] * b;
                                float L_2 = newy <= MAXVALF ? Color::cachefy[newy] : 327.68f * (116.f * xcbrtf(newy / MAXVALF) - 16.f);

                                //gamut control
                                if (settings->rgbcurveslumamode_gamut) {
                                    float Lpro = L_2 / 327.68f;
                                    float Chpro = sqrtf (SQR (a_1) + SQR (b_1)) / 327.68f;
                                    float HH = NAN; // we set HH to NAN, because then it will be calculated in Color::gamutLchonly only if needed
//                                    float HH = xatan2f(b_1, a_1);
                                    // According to mathematical laws we can get the sin and cos of HH by simple operations even if we don't calculate HH
                                    float2 sincosval;

                                    if (Chpro == 0.0f) {
                                        sincosval.y = 1.0f;
                                        sincosval.x = 0.0f;
                                    } else {
                                        sincosval.y = a_1 / (Chpro * 327.68f);
                                        sincosval.x = b_1 / (Chpro * 327.68f);
                                    }

#ifdef _DEBUG
                                    bool neg = false;
                                    bool more_rgb = false;
                                    //gamut control : Lab values are in gamut
                                    Color::gamutLchonly (HH, sincosval, Lpro, Chpro, r, g, b, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
                                    //gamut control : Lab values are in gamut
                                    Color::gamutLchonly (HH, sincosval, Lpro, Chpro, r, g, b, wip, highlight, 0.15f, 0.96f);
#endif
                                    //end of gamut control
                                } else {
                                    float x_, y_, z_;
                                    //calculate RGB with L_2 and old value of a and b
                                    Color::Lab2XYZ (L_2, a_1, b_1, x_, y_, z_) ;
                                    Color::xyz2rgb (x_, y_, z_, r, g, b, wip);
                                }

                                setUnlessOOG(rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj], r, g, b);
                            }
                        }
                    }
                }

                if (editID == EUID_HSV_H || editID == EUID_HSV_S || editID == EUID_HSV_V) {
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            float h, s, v;
                            Color::rgb2hsv (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj], h, s, v);
                            editWhateverTmp[ti * TS + tj] = h;
                        }
                    }
                }

                if (isProPhoto) { // this is a hack to avoid the blue=>black bug (Issue 2141)
                    proPhotoBlue(rtemp, gtemp, btemp, istart, tH, jstart, tW, TS);
                }

                // filling the pipette buffer
                if (editID == EUID_BlackWhiteLuminance) {
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            float X, Y, Z, L, aa, bb;
                            //rgb=>lab
                            Color::rgbxyz (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj], X, Y, Z, wp);
                            //convert Lab
                            Color::XYZ2Lab (X, Y, Z, L, aa, bb);
                            //end rgb=>lab
                            float HH = xatan2f (bb, aa); // HH hue in -3.141  +3.141

                            editWhateverTmp[ti * TS + tj] = float (Color::huelab_to_huehsv2 (HH));
                        }
                    }
                }

                //black and white
                if (blackwhite) {
                    if (algm == 0) { //lightness
                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {

                                float r = rtemp[ti * TS + tj];
                                float g = gtemp[ti * TS + tj];
                                float b = btemp[ti * TS + tj];

                                // --------------------------------------------------

                                // Method 1: Luminosity (code taken from Gimp)
                                /*
                                float maxi = max(r, g, b);
                                float mini = min(r, g, b);
                                r = g = b = (maxi+mini)/2;
                                */

                                // Method 2: Luminance (former RT code)
                                r = g = b = (0.299f * r + 0.587f * g + 0.114f * b);

                                // --------------------------------------------------

#ifndef __SSE2__

                                //gamma correction: pseudo TRC curve
                                if (hasgammabw) {
                                    Color::trcGammaBW (r, g, b, gammabwr, gammabwg, gammabwb);
                                }

#endif
                                rtemp[ti * TS + tj] = r;
                                gtemp[ti * TS + tj] = g;
                                btemp[ti * TS + tj] = b;
                            }

#ifdef __SSE2__

                            if (hasgammabw) {
                                //gamma correction: pseudo TRC curve
                                Color::trcGammaBWRow (&rtemp[ti * TS], &gtemp[ti * TS], &btemp[ti * TS], tW - jstart, gammabwr, gammabwg, gammabwb);
                            }

#endif

                        }
                    } else if (algm == 1) { //Luminance mixer in Lab mode to avoid artifacts
                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                                //rgb => xyz
                                float X, Y, Z;
                                Color::rgbxyz (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj], X, Y, Z, wp);
                                //xyz => Lab
                                float L, aa, bb;
                                Color::XYZ2Lab (X, Y, Z, L, aa, bb);
                                float CC = sqrtf (SQR (aa) + SQR (bb)) / 327.68f; //CC chromaticity in 0..180 or more
                                float HH = xatan2f (bb, aa); // HH hue in -3.141  +3.141
                                float2 sincosval;

                                if (CC == 0.0f) {
                                    sincosval.y = 1.f;
                                    sincosval.x = 0.0f;
                                } else {
                                    sincosval.y = aa / (CC * 327.68f);
                                    sincosval.x = bb / (CC * 327.68f);
                                }

                                if (bwlCurveEnabled) {
                                    L /= 32768.f;
                                    double hr = Color::huelab_to_huehsv2 (HH);
                                    float valparam = float ((bwlCurve->getVal (hr) - 0.5f) * 2.0f); //get l_r=f(H)
                                    float kcc = (CC / 70.f); //take Chroma into account...70 "middle" of chromaticity (arbitrary and simple), one can imagine other algorithme
                                    //reduct action for low chroma and increase action for high chroma
                                    valparam *= kcc;

                                    if (valparam > 0.f) {
                                        L = (1.f - valparam) * L + valparam * (1.f - SQR (SQR (SQR (SQR (1.f - min (L, 1.0f)))))); // SQR (SQR((SQR)  to increase action in low light
                                    } else {
                                        L *= (1.f + valparam);    //for negative
                                    }

                                    L *= 32768.f;
                                }

                                float RR, GG, BB;
                                L /= 327.68f;
#ifdef _DEBUG
                                bool neg = false;
                                bool more_rgb = false;
                                //gamut control : Lab values are in gamut
                                Color::gamutLchonly (HH, sincosval, L, CC, RR, GG, BB, wip, highlight, 0.15f, 0.96f, neg, more_rgb);
#else
                                //gamut control : Lab values are in gamut
                                Color::gamutLchonly (HH, sincosval, L, CC, RR, GG, BB, wip, highlight, 0.15f, 0.96f);
#endif
                                L *= 327.68f;
                                //convert l => rgb
                                Color::L2XYZ (L, X, Y, Z);
                                float newRed; // We use the red channel for bw
                                Color::xyz2r (X, Y, Z, newRed, wip);
                                rtemp[ti * TS + tj] = gtemp[ti * TS + tj] = btemp[ti * TS + tj] = newRed;
#ifndef __SSE2__

                                if (hasgammabw) {
                                    //gamma correction: pseudo TRC curve
                                    Color::trcGammaBW (rtemp[ti * TS + tj], gtemp[ti * TS + tj], btemp[ti * TS + tj], gammabwr, gammabwg, gammabwb);
                                }

#endif
                            }

#ifdef __SSE2__

                            if (hasgammabw) {
                                //gamma correction: pseudo TRC curve
                                Color::trcGammaBWRow (&rtemp[ti * TS], &gtemp[ti * TS], &btemp[ti * TS], tW - jstart, gammabwr, gammabwg, gammabwb);
                            }

#endif
                        }
                    }
                }

                if (!blackwhite) {
                    if (editImgFloat || editWhatever) {
                        for (int i = istart, ti = 0; i < tH; i++, ti++) {
                            for (int j = jstart, tj = 0; j < tW; j++, tj++) {

                                // filling the pipette buffer by the content of the temp pipette buffers
                                if (editImgFloat) {
                                    editImgFloat->r (i, j) = editIFloatTmpR[ti * TS + tj];
                                    editImgFloat->g (i, j) = editIFloatTmpG[ti * TS + tj];
                                    editImgFloat->b (i, j) = editIFloatTmpB[ti * TS + tj];
                                } else if (editWhatever) {
                                    editWhatever->v (i, j) = editWhateverTmp[ti * TS + tj];
                                }
                            }
                        }
                    }
                    // ready, fill lab
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            int idx = ti * TS + tj;
                            working->r(i, j) = rtemp[idx];
                            working->g(i, j) = gtemp[idx];
                            working->b(i, j) = btemp[idx];
                        }
                    }
                } else { // black & white
                    // Auto channel mixer needs whole image, so we now copy to tmpImage and close the tiled processing
                    for (int i = istart, ti = 0; i < tH; i++, ti++) {
                        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
                            // filling the pipette buffer by the content of the temp pipette buffers
                            if (editImgFloat) {
                                editImgFloat->r (i, j) = editIFloatTmpR[ti * TS + tj];
                                editImgFloat->g (i, j) = editIFloatTmpG[ti * TS + tj];
                                editImgFloat->b (i, j) = editIFloatTmpB[ti * TS + tj];
                            } else if (editWhatever) {
                                editWhatever->v (i, j) = editWhateverTmp[ti * TS + tj];
                            }

                            tmpImage->r (i, j) = rtemp[ti * TS + tj];
                            tmpImage->g (i, j) = gtemp[ti * TS + tj];
                            tmpImage->b (i, j) = btemp[ti * TS + tj];
                        }
                    }
                }
            }
        }

        if (editIFloatBuffer) {
            free (editIFloatBuffer);
        }

        if (editWhateverBuffer) {
            free (editWhateverBuffer);
        }
    }

    // starting a new tile processing with a 'reduction' clause for the auto mixer computing
    if (blackwhite) {//channel-mixer
        int tW = working->getWidth();
        int tH = working->getHeight();

        if (algm == 2) { //channel-mixer
            float filcor;
            double rrm, ggm, bbm;
            Color::computeBWMixerConstants (params->blackwhite.setting, params->blackwhite.filter, "", filcor,
                                            bwr, bwg, bwb, kcorec, rrm, ggm, bbm);

#ifdef _OPENMP
            #pragma omp parallel for schedule(dynamic, 16)
#endif

            for (int i = 0; i < tH; i++) {
                for (int j = 0; j < tW; j++) {

                    //mix channel
                    tmpImage->r (i, j) = tmpImage->g (i, j) = tmpImage->b (i, j) = /*CLIP*/ ((bwr * tmpImage->r (i, j) + bwg * tmpImage->g (i, j) + bwb * tmpImage->b (i, j)) * kcorec);

#ifndef __SSE2__

                    //gamma correction: pseudo TRC curve
                    if (hasgammabw) {
                        Color::trcGammaBW (tmpImage->r (i, j), tmpImage->g (i, j), tmpImage->b (i, j), gammabwr, gammabwg, gammabwb);
                    }

#endif
                }

#ifdef __SSE2__

                if (hasgammabw) {
                    //gamma correction: pseudo TRC curve
                    Color::trcGammaBWRow (tmpImage->r (i), tmpImage->g (i), tmpImage->b (i), tW, gammabwr, gammabwg, gammabwb);
                }

#endif
            }
        }

        // ready, fill lab (has to be the same code than the "fill lab" above!)
#ifdef _OPENMP
        #pragma omp parallel for schedule(dynamic, 5)
#endif

        for (int i = 0; i < tH; i++) {
            for (int j = 0; j < tW; j++) {
                working->r(i, j) = tmpImage->r(i, j);
                working->g(i, j) = tmpImage->g(i, j);
                working->b(i, j) = tmpImage->b(i, j);
            }
        }
    }

    if (tmpImage) {
        delete tmpImage;
    }

    hslEqualizer(working);
}

#else // if 0

namespace {

void proPhotoBlue(Imagefloat *rgb, bool multiThread)
{
    // TODO
    const int W = rgb->getWidth();
    const int H = rgb->getHeight();
#ifdef _OPENMP
#   pragma omp parallel for if (multiThread)
#endif
    for (int y = 0; y < H; ++y) {
        int x = 0;
#ifdef __SSE2__
        for (; x < W - 3; x += 4) {
            vfloat rv = LVF(rgb->r(y, x));
            vfloat gv = LVF(rgb->g(y, x));
            vmask zeromask = vorm(vmaskf_eq(rv, ZEROV), vmaskf_eq(gv, ZEROV));
            if (_mm_movemask_ps((vfloat)zeromask)) {
                for (int k = 0; k < 4; ++k) {
                    float r = rgb->r(y, x+k);
                    float g = rgb->g(y, x+k);
                    float b = rgb->b(y, x+k);
                    
                    if ((r == 0.0f || g == 0.0f) && rtengine::min(r, g, b) >= 0.f) {
                        float h, s, v;
                        Color::rgb2hsv(r, g, b, h, s, v);
                        s *= 0.99f;
                        Color::hsv2rgb(h, s, v, rgb->r(y, x+k), rgb->g(y, x+k), rgb->b(y, x+k));
                    }
                }
            }
        }
#endif
        for (; x < W; ++x) {
            float r = rgb->r(y, x);
            float g = rgb->g(y, x);
            float b = rgb->b(y, x);

            if ((r == 0.0f || g == 0.0f) && rtengine::min(r, g, b) >= 0.f) {
                float h, s, v;
                Color::rgb2hsv(r, g, b, h, s, v);
                s *= 0.99f;
                Color::hsv2rgb(h, s, v, rgb->r(y, x), rgb->g(y, x), rgb->b(y, x));
            }
        }
    }
}

} // namespace


void ImProcFunctions::rgbProc(Imagefloat *rgb)
{
    channelMixer(rgb);
    exposure(rgb);
    hslEqualizer(rgb);
    toneEqualizer(rgb);
    rgbCurves(rgb);
    if (params->icm.workingProfile == "ProPhoto") {
        proPhotoBlue(rgb, multiThread);
    }
    blackAndWhite(rgb);
}

#endif

void ImProcFunctions::impulsedenoise(Imagefloat *rgb)
{

    if (params->impulseDenoise.enabled && rgb->getWidth() >= 8 && rgb->getHeight() >= 8)

    {
        rgb->setMode(Imagefloat::Mode::LAB, multiThread);
        impulse_nr(rgb, (float)params->impulseDenoise.thresh / 20.0 );
    }
}

void ImProcFunctions::defringe(Imagefloat *rgb)
{
    if (params->defringe.enabled && rgb->getWidth() >= 8 && rgb->getHeight() >= 8)

    {
        rgb->setMode(Imagefloat::Mode::LAB, multiThread);
        PF_correct_RT(rgb, params->defringe.radius, params->defringe.threshold);
    }
}

void ImProcFunctions::getAutoExp  (const LUTu &histogram, int histcompr, double clip,
                                   double& expcomp, int& bright, int& contr, int& black, int& hlcompr, int& hlcomprthresh)
{

    float scale = 65536.0f;
    float midgray = 0.1842f; //middle gray in linear gamma =1 50% luminance

    int imax = 65536 >> histcompr;
    int overex = 0;
    float sum = 0.f, hisum = 0.f, losum = 0.f;
    float ave = 0.f, hidev = 0.f, lodev = 0.f;

    //find average luminance
    histogram.getSumAndAverage (sum, ave);

    //find median of luminance
    size_t median = 0, count = histogram[0];

    while (count < sum / 2) {
        median++;
        count += histogram[median];
    }

    if (median == 0 || ave < 1.f) { //probably the image is a blackframe
        expcomp = 0.;
        black = 0;
        bright = 0;
        contr = 0;
        hlcompr = 0;
        hlcomprthresh = 0;
        return;
    }

    // compute std dev on the high and low side of median
    // and octiles of histogram
    float octile[8] = {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f}, ospread = 0.f;
    count = 0;

    int i = 0;

    for (; i < min ((int)ave, imax); i++) {
        if (count < 8) {
            octile[count] += histogram[i];

            if (octile[count] > sum / 8.f || (count == 7 && octile[count] > sum / 16.f)) {
                octile[count] = xlog (1. + (float)i) / log (2.f);
                count++;// = min(count+1,7);
            }
        }

        //lodev += SQR(ave-i)*histogram[i];
        lodev += (xlog (ave + 1.f) - xlog ((float)i + 1.)) * histogram[i];
        losum += histogram[i];
    }

    for (; i < imax; i++) {
        if (count < 8) {
            octile[count] += histogram[i];

            if (octile[count] > sum / 8.f || (count == 7 && octile[count] > sum / 16.f)) {
                octile[count] = xlog (1. + (float)i) / log (2.f);
                count++;// = min(count+1,7);
            }
        }

        //hidev += SQR(i-ave)*histogram[i];
        hidev += (xlog ((float)i + 1.) - xlog (ave + 1.f)) * histogram[i];
        hisum += histogram[i];

    }

    if (losum == 0 || hisum == 0) { //probably the image is a blackframe
        expcomp = 0.;
        black = 0;
        bright = 0;
        contr = 0;
        hlcompr = 0;
        hlcomprthresh = 0;
        return;
    }

//    lodev = (lodev / (log(2.f) * losum));
//    hidev = (hidev / (log(2.f) * hisum));

    if (octile[6] > log1p ((float)imax) / log2 (2.f)) { //if very overxposed image
        octile[6] = 1.5f * octile[5] - 0.5f * octile[4];
        overex = 2;
    }

    if (octile[7] > log1p ((float)imax) / log2 (2.f)) { //if overexposed
        octile[7] = 1.5f * octile[6] - 0.5f * octile[5];
        overex = 1;
    }

    // store values of octile[6] and octile[7] for calculation of exposure compensation
    // if we don't do this and the pixture is underexposed, calculation of exposure compensation assumes
    // that it's overexposed and calculates the wrong direction
    float oct6, oct7;
    oct6 = octile[6];
    oct7 = octile[7];


    for (int i = 1; i < 8; i++) {
        if (octile[i] == 0.0f) {
            octile[i] = octile[i - 1];
        }
    }

    // compute weighted average separation of octiles
    // for future use in contrast setting
    for (int i = 1; i < 6; i++) {
        ospread += (octile[i + 1] - octile[i]) / max (0.5f, (i > 2 ? (octile[i + 1] - octile[3]) : (octile[3] - octile[i])));
    }

    ospread /= 5.f;

    if (ospread <= 0.f) { //probably the image is a blackframe
        expcomp = 0.;
        black = 0;
        bright = 0;
        contr = 0;
        hlcompr = 0;
        hlcomprthresh = 0;
        return;
    }


    // compute clipping points based on the original histograms (linear, without exp comp.)
    unsigned int clipped = 0;
    int rawmax = (imax) - 1;

    while (histogram[rawmax] + clipped <= 0 && rawmax > 1) {
        clipped += histogram[rawmax];
        rawmax--;
    }

    //compute clipped white point
    unsigned int clippable = (int) (sum * clip / 100.f );
    clipped = 0;
    int whiteclip = (imax) - 1;

    while (whiteclip > 1 && (histogram[whiteclip] + clipped) <= clippable) {
        clipped += histogram[whiteclip];
        whiteclip--;
    }

    //compute clipped black point
    clipped = 0;
    int shc = 0;

    while (shc < whiteclip - 1 && histogram[shc] + clipped <= clippable) {
        clipped += histogram[shc];
        shc++;
    }

    //rescale to 65535 max
    rawmax <<= histcompr;
    whiteclip <<= histcompr;
    ave = ave * (1 << histcompr);
    median <<= histcompr;
    shc <<= histcompr;

//    //prevent division by 0
//    if (lodev == 0.f) {
//        lodev = 1.f;
//    }

    //compute exposure compensation as geometric mean of the amount that
    //sets the mean or median at middle gray, and the amount that sets the estimated top
    //of the histogram at or near clipping.
    //float expcomp1 = (log(/*(median/ave)*//*(hidev/lodev)*/midgray*scale/(ave-shc+midgray*shc))+log((hidev/lodev)))/log(2.f);
    float expcomp1 = (log (/*(median/ave)*//*(hidev/lodev)*/midgray * scale / (ave - shc + midgray * shc))) / log (2.f);
    float expcomp2;

    if (overex == 0) { // image is not overexposed
        expcomp2 = 0.5f * ( (15.5f - histcompr - (2.f * oct7 - oct6)) + log (scale / rawmax) / log (2.f) );
    } else {
        expcomp2 = 0.5f * ( (15.5f - histcompr - (2.f * octile[7] - octile[6])) + log (scale / rawmax) / log (2.f) );
    }

    if (fabs (expcomp1) - fabs (expcomp2) > 1.f) { //for great expcomp
        expcomp = (expcomp1 * fabs (expcomp2) + expcomp2 * fabs (expcomp1)) / (fabs (expcomp1) + fabs (expcomp2));
    } else {
        expcomp = 0.5 * (double)expcomp1 + 0.5 * (double) expcomp2; //for small expcomp
    }

    float gain = exp ((float)expcomp * log (2.f));

    float corr = sqrt (gain * scale / rawmax);
    black = (int) shc * corr;


    //now tune hlcompr to bring back rawmax to 65535
    hlcomprthresh = 0;
    //this is a series approximation of the actual formula for comp,
    //which is a transcendental equation
    float comp = (gain * ((float)whiteclip) / scale - 1.f) * 2.3f; // 2.3 instead of 2 to increase slightly comp
    hlcompr = (int) (100.*comp / (max (0.0, expcomp) + 1.0));
    hlcompr = max (0, min (100, hlcompr));

    //now find brightness if gain didn't bring ave to midgray using
    //the envelope of the actual 'control cage' brightness curve for simplicity
    float midtmp = gain * sqrt (median * ave) / scale;

    if (midtmp < 0.1f) {
        bright = (midgray - midtmp) * 15.0 / (midtmp);
    } else {
        bright = (midgray - midtmp) * 15.0 / (0.10833 - 0.0833 * midtmp);
    }

    bright = 0.25 */*(median/ave)*(hidev/lodev)*/max (0, bright);

    //compute contrast that spreads the average spacing of octiles
    contr = (int) 50.0f * (1.1f - ospread);
    contr = max (0, min (100, contr));
    //take gamma into account
    double whiteclipg = (int) (CurveFactory::gamma2 (whiteclip * corr / 65536.0) * 65536.0);

    float gavg = 0.;

    float val = 0.f;
    float increment = corr * (1 << histcompr);

    for (int i = 0; i < 65536 >> histcompr; i++) {
        gavg += histogram[i] * Color::gamma2curve[val];
        val += increment;
    }

    gavg /= sum;

    if (black < gavg) {
        int maxwhiteclip = (gavg - black) * 4 / 3 + black; // don't let whiteclip be such large that the histogram average goes above 3/4

        if (whiteclipg < maxwhiteclip) {
            whiteclipg = maxwhiteclip;
        }
    }

    whiteclipg = CurveFactory::igamma2 ((float) (whiteclipg / 65535.0)) * 65535.0; //need to inverse gamma transform to get correct exposure compensation parameter

    //corection with gamma
    black = (int) ((65535 * black) / whiteclipg);
    //expcomp = log(65535.0 / (whiteclipg)) / log(2.0);

    //diagnostics
    //printf ("**************** AUTO LEVELS ****************\n");
    /*
    if (settings->verbose) {
        printf("sum=%i clip=%f clippable=%i  clipWh=%i  clipBl=%i\n",somm, clip, clippable,clipwh, clipbl);
        printf ("expcomp1= %f   expcomp2= %f gain= %f  expcomp=%f\n",expcomp1,expcomp2,gain,expcomp);
        printf ("expo=%f\n",expo);
        printf ("median: %i  average: %f    median/average: %f\n",median,ave, median/ave);
        printf ("average: %f\n",ave);
        printf("comp=%f hlcomp: %i\n",comp, hlcompr);
        printf ("median/average: %f\n",median/ave);
        printf ("lodev: %f   hidev: %f      hidev/lodev: %f\n",lodev,hidev,hidev/lodev);
        printf ("lodev: %f\n",lodev);
        printf ("hidev: %f\n",hidev);
        printf ("imax=%d rawmax= %d  whiteclip= %d  gain= %f\n",imax,rawmax,whiteclip,gain);
        printf ("octiles: %f %f %f %f %f %f %f %f\n",octile[0],octile[1],octile[2],octile[3],octile[4],octile[5],octile[6],octile[7]);
        printf ("ospread= %f\n",ospread);
        printf ("overexp= %i\n",overex);
    }
    */
    /*
     // %%%%%%%%%% LEGACY AUTOEXPOSURE CODE %%%%%%%%%%%%%
     // black point selection is based on the linear result (yielding better visual results)
     black = (int)(shc * corr);
     // compute the white point of the exp. compensated gamma corrected image
     double whiteclipg = (int)(CurveFactory::gamma2 (whiteclip * corr / 65536.0) * 65536.0);

     // compute average intensity of the exp compensated, gamma corrected image
     double gavg = 0;
     for (int i=0; i<65536>>histcompr; i++)
     gavg += histogram[i] * CurveFactory::gamma2((int)(corr*(i<<histcompr)<65535 ? corr*(i<<histcompr) : 65535)) / sum;

     if (black < gavg) {
     int maxwhiteclip = (gavg - black) * 4 / 3 + black; // don't let whiteclip be such large that the histogram average goes above 3/4
     //double mavg = 65536.0 / (whiteclipg-black) * (gavg - black);
     if (whiteclipg < maxwhiteclip)
     whiteclipg = maxwhiteclip;
     }

     whiteclipg = CurveFactory::igamma2 ((float)(whiteclipg/65535.0)) * 65535.0; //need to inverse gamma transform to get correct exposure compensation parameter

     black = (int)((65535*black)/whiteclipg);
     expcomp = log(65535.0 / (whiteclipg)) / log(2.0);

     if (expcomp<0.0)   expcomp = 0.0;*/
    if (expcomp < -5.0) {
        expcomp = -5.0;
    }

    if (expcomp > 12.0) {
        expcomp = 12.0;
    }

    bright = max (-100, min (bright, 100));

}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

double ImProcFunctions::getAutoDistor  (const Glib::ustring &fname, int thumb_size)
{
    if (fname != "") {
        int w_raw = -1, h_raw = thumb_size;
        int w_thumb = -1, h_thumb = thumb_size;

        eSensorType sensorType = rtengine::ST_NONE;
        Thumbnail* thumb = rtengine::Thumbnail::loadQuickFromRaw (fname, sensorType, w_thumb, h_thumb, 1, FALSE);

        if (!thumb) {
            return 0.0;
        }

        Thumbnail* raw =   rtengine::Thumbnail::loadFromRaw(fname, sensorType, w_raw, h_raw, 1, 1.0, FALSE);

        if (!raw) {
            delete thumb;
            return 0.0;
        }

        if (h_thumb != h_raw) {
            delete thumb;
            delete raw;
            return 0.0;
        }

        int width;

        if (w_thumb > w_raw) {
            width = w_raw;
        } else {
            width = w_thumb;
        }

        unsigned char* thumbGray;
        unsigned char* rawGray;
        thumbGray = thumb->getGrayscaleHistEQ (width);
        rawGray = raw->getGrayscaleHistEQ (width);

        if (!thumbGray || !rawGray) {
            if (thumbGray) {
                delete thumbGray;
            }

            if (rawGray) {
                delete rawGray;
            }

            delete thumb;
            delete raw;
            return 0.0;
        }

        double dist_amount;
        int dist_result = calcDistortion (thumbGray, rawGray, width, h_thumb, 1, dist_amount);

        if (dist_result == -1) { // not enough features found, try increasing max. number of features by factor 4
            calcDistortion (thumbGray, rawGray, width, h_thumb, 4, dist_amount);
        }

        delete thumbGray;
        delete rawGray;
        delete thumb;
        delete raw;
        return dist_amount;
    } else {
        return 0.0;
    }
}

void ImProcFunctions::rgb2lab (Imagefloat &src, LabImage &dst, const Glib::ustring &workingSpace)
{
    src.assignColorSpace(workingSpace);
    src.toLab(dst, true);
}

void ImProcFunctions::lab2rgb (const LabImage &src, Imagefloat &dst, const Glib::ustring &workingSpace)
{
    dst.assignColorSpace(workingSpace);
    dst.assignMode(Imagefloat::Mode::RGB);
    
    TMatrix wiprof = ICCStore::getInstance()->workingSpaceInverseMatrix ( workingSpace );
    const float wip[3][3] = {
        {static_cast<float> (wiprof[0][0]), static_cast<float> (wiprof[0][1]), static_cast<float> (wiprof[0][2])},
        {static_cast<float> (wiprof[1][0]), static_cast<float> (wiprof[1][1]), static_cast<float> (wiprof[1][2])},
        {static_cast<float> (wiprof[2][0]), static_cast<float> (wiprof[2][1]), static_cast<float> (wiprof[2][2])}
    };

    const int W = dst.getWidth();
    const int H = dst.getHeight();
#ifdef __SSE2__
    vfloat wipv[3][3];

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            wipv[i][j] = F2V (wiprof[i][j]);
        }
    }

#endif

#ifdef _OPENMP
    #pragma omp parallel for schedule(dynamic,16)
#endif

    for (int i = 0; i < H; i++) {
        int j = 0;
#ifdef __SSE2__

        for (; j < W - 3; j += 4) {
            vfloat X, Y, Z;
            vfloat R, G, B;
            Color::Lab2XYZ (LVFU (src.L[i][j]), LVFU (src.a[i][j]), LVFU (src.b[i][j]), X, Y, Z);
            Color::xyz2rgb (X, Y, Z, R, G, B, wipv);
            STVFU (dst.r (i, j), R);
            STVFU (dst.g (i, j), G);
            STVFU (dst.b (i, j), B);
        }

#endif

        for (; j < W; j++) {
            float X, Y, Z;
            Color::Lab2XYZ (src.L[i][j], src.a[i][j], src.b[i][j], X, Y, Z);
            Color::xyz2rgb (X, Y, Z, dst.r (i, j), dst.g (i, j), dst.b (i, j), wip);
        }
    }
}

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


void ImProcFunctions::setViewport(int ox, int oy, int fw, int fh)
{
    offset_x = ox;
    offset_y = oy;
    full_width = fw;
    full_height = fh;
}


void ImProcFunctions::setOutputHistograms(LUTu *histToneCurve, LUTu *histCCurve, LUTu *histLCurve)
{
    hist_tonecurve = histToneCurve;
    hist_ccurve = histCCurve;
    hist_lcurve = histLCurve;
}


void ImProcFunctions::setShowSharpeningMask(bool yes)
{
    show_sharpening_mask = yes;
}


bool ImProcFunctions::process(Pipeline pipeline, Stage stage, Imagefloat *img)
{
    bool stop = false;
    switch (stage) {
    case Stage::STAGE_0:
        dehaze(img);
        dynamicRangeCompression(img);
        break;
    case Stage::STAGE_1:
        rgbProc(img);
        break;
    case Stage::STAGE_2:
        if (pipeline == Pipeline::OUTPUT ||
            (pipeline == Pipeline::PREVIEW && scale == 1)) {
            stop = sharpening(img, params->sharpening, show_sharpening_mask);
            if (!stop) {
                impulsedenoise(img);
                defringe(img);
            }
        }
        stop = stop || colorCorrection(img, offset_x, offset_y, full_width, full_height);
        stop = stop || guidedSmoothing(img, offset_x, offset_y, full_width, full_height);
        break;
    case Stage::STAGE_3:
        logEncoding(img, hist_tonecurve);
        labAdjustments(img, hist_ccurve, hist_lcurve);
        stop = stop || textureBoost(img, offset_x, offset_y, full_width, full_height);
        if (pipeline != Pipeline::THUMBNAIL) {
            stop = stop || contrastByDetailLevels(img, offset_x, offset_y, full_width, full_height);
        }
        if (!stop) {
            softLight(img);
            localContrast(img);
            filmGrain(img, offset_x, offset_y, full_width, full_height);
        }
        break;
    }
    return stop;
}

} // namespace rtengine
