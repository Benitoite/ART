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
#include "rtengine.h"
#include "colortemp.h"
#include "imagesource.h"
#include "improcfun.h"
#include "curves.h"
#include "iccstore.h"
#include "clutstore.h"
#include "processingjob.h"
#include <glibmm.h>
#include "../rtgui/options.h"
#include "rawimagesource.h"
#include "../rtgui/multilangmgr.h"
#include "mytime.h"
#include "rescale.h"
#undef THREAD_PRIORITY_NORMAL

namespace rtengine
{
extern const Settings* settings;

namespace
{

template <typename T>
void adjust_radius (const T &default_param, double scale_factor, T &param)
{
    const double delta = (param - default_param) * scale_factor;
    param = default_param + delta;
}


class ImageProcessor
{
public:
    ImageProcessor(
        ProcessingJob* pjob,
        int& errorCode,
        ProgressListener* pl,
        bool flush
    ) :
        job (static_cast<ProcessingJobImpl*> (pjob)),
        errorCode (errorCode),
        pl (pl),
        flush (flush),
        // internal state
        ii(nullptr),
        imgsrc(nullptr),
        fw(0),
        fh(0),
        oX(0),
        oY(0),
        oW(0),
        oH(0),
        tr(0),
        pp(0, 0, 0, 0, 0),
        tilesize(0),
        overlap(0),
        dnstore(),
        expcomp(0.0),
        bright(0),
        contr(0),
        black(0),
        hlcompr(0),
        hlcomprthresh(0),
        baseImg(nullptr),
        labView(nullptr),
        autili(false),
        butili(false)
    {
    }

    Imagefloat *operator()()
    {
        if (!job->fast) {
            return normal_pipeline();
        } else {
            return fast_pipeline();
        }
    }

private:
    Imagefloat *normal_pipeline()
    {
        if (settings->verbose) {
            std::cout << "Processing with the normal pipeline" << std::endl;
        }
        
        if (!stage_init(false)) {
            return nullptr;
        }

        stage_denoise();
        stage_transform();
        return stage_finish();
    }

    Imagefloat *fast_pipeline()
    {
        if (!job->pparams.resize.enabled) {
            return normal_pipeline();
        }

        if (settings->verbose) {
            std::cout << "Processing with the fast pipeline" << std::endl;
        }

        pl = nullptr;

        if (!stage_init(true)) {
            return nullptr;
        }

        stage_transform();
        stage_early_resize();
        stage_denoise();
        return stage_finish();
    }

    bool stage_init(bool is_fast)
    {
        errorCode = 0;

        if (pl) {
            pl->setProgressStr ("PROGRESSBAR_PROCESSING");
            pl->setProgress (0.0);
        }

        ii = job->initialImage;

        if (!ii) {
            ii = InitialImage::load (job->fname, job->isRaw, &errorCode);

            if (errorCode) {
                delete job;
                return false; //return nullptr;
            }
        }

        procparams::ProcParams& params = job->pparams;

        // acquire image from imagesource
        imgsrc = ii->getImageSource ();

        tr = getCoarseBitMask (params.coarse);
        if(imgsrc->getSensorType() == ST_BAYER) {
            if(params.raw.bayersensor.method!= RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::PIXELSHIFT)) {
                imgsrc->setBorder(params.raw.bayersensor.border);
            } else {
                imgsrc->setBorder(std::max(params.raw.bayersensor.border, 2));
            }
        }
        imgsrc->getFullSize (fw, fh, tr);
        oW = fw;
        oH = fh;

        // check the crop params
        if (params.crop.x > fw || params.crop.y > fh) {
            // the crop is completely out of the image, so we disable the crop
            params.crop.enabled = false;
            // and we set the values to the defaults
            params.crop.x = 0;
            params.crop.y = 0;
            params.crop.w = fw;
            params.crop.h = fh;
        } else {
            if (params.crop.x < 0) {
                params.crop.x = 0;
            }

            if (params.crop.y < 0) {
                params.crop.y = 0;
            }

            if ((params.crop.x + params.crop.w) > fw) {
                // crop overflow in the width dimension ; we trim it
                params.crop.w = fw - params.crop.x;
            }

            if ((params.crop.y + params.crop.h) > fh) {
                // crop overflow in the height dimension ; we trim it
                params.crop.h = fh - params.crop.y;
            }
        }

//    MyTime t1,t2;
//    t1.set();

        ipf_p.reset (new ImProcFunctions (&params, true));
        ImProcFunctions &ipf = * (ipf_p.get());
        double scale_factor = 1.0;
        if (is_fast) {
            int imw, imh;
            scale_factor = ipf.resizeScale (&params, fw, fh, imw, imh);
            adjust_procparams(scale_factor);
        }
        

        imgsrc->setCurrentFrame (params.raw.bayersensor.imageNum);
        imgsrc->preprocess ( params.raw, params.lensProf, params.coarse, params.denoise.enabled);

        if (pl) {
            pl->setProgress (0.20);
        }
        bool autoContrast = imgsrc->getSensorType() == ST_BAYER ? params.raw.bayersensor.dualDemosaicAutoContrast : params.raw.xtranssensor.dualDemosaicAutoContrast;
        double contrastThreshold = imgsrc->getSensorType() == ST_BAYER ? params.raw.bayersensor.dualDemosaicContrast : params.raw.xtranssensor.dualDemosaicContrast;

        imgsrc->demosaic (params.raw, autoContrast, contrastThreshold);


        if (pl) {
            pl->setProgress (0.30);
        }
        pp = PreviewProps (0, 0, fw, fh, 1);

        if (params.retinex.enabled) { //enabled Retinex
            LUTf cdcurve (65536, 0);
            LUTf mapcurve (65536, 0);
            LUTu dummy;
            RetinextransmissionCurve dehatransmissionCurve;
            RetinexgaintransmissionCurve dehagaintransmissionCurve;
            bool dehacontlutili = false;
            bool mapcontlutili = false;
            bool useHsl = false;
//        multi_array2D<float, 3> conversionBuffer(1, 1);
            multi_array2D<float, 4> conversionBuffer (1, 1);
            imgsrc->retinexPrepareBuffers (params.icm, params.retinex, conversionBuffer, dummy);
            imgsrc->retinexPrepareCurves (params.retinex, cdcurve, mapcurve, dehatransmissionCurve, dehagaintransmissionCurve, dehacontlutili, mapcontlutili, useHsl, dummy, dummy );
            float minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax;
            imgsrc->retinex ( params.icm, params.retinex, params.toneCurve, cdcurve, mapcurve, dehatransmissionCurve, dehagaintransmissionCurve, conversionBuffer, dehacontlutili, mapcontlutili, useHsl, minCD, maxCD, mini, maxi, Tmean, Tsigma, Tmin, Tmax, dummy);
        }

        if (pl) {
            pl->setProgress (0.40);
        }

        imgsrc->HLRecovery_Global ( params.toneCurve );


        if (pl) {
            pl->setProgress (0.45);
        }

        // set the color temperature
        currWB = ColorTemp (params.wb.temperature, params.wb.green, params.wb.equal, params.wb.method);

        if (!params.wb.enabled) {
            currWB = ColorTemp();
        } else if (params.wb.method == "Camera") {
            currWB = imgsrc->getWB ();
        } else if (params.wb.method == "Auto") {
            double rm, gm, bm;
            imgsrc->getAutoWBMultipliers (rm, gm, bm);
            currWB.update (rm, gm, bm, params.wb.equal, params.wb.tempBias);
        }

        if (params.denoise.enabled) {
            ipf.denoiseComputeParams(imgsrc, currWB, dnstore, params.denoise);
        }
        
        baseImg = new Imagefloat (fw, fh);
        imgsrc->getImage (currWB, tr, baseImg, pp, params.toneCurve, params.raw);

        if (pl) {
            pl->setProgress (0.50);
        }

//  LUTf Noisecurve (65536,0);
//!!!// auto exposure!!!
        expcomp = params.toneCurve.expcomp;
        bright = params.toneCurve.brightness;
        contr = params.toneCurve.contrast;
        black = params.toneCurve.black;
        hlcompr = params.toneCurve.hlcompr;
        hlcomprthresh = params.toneCurve.hlcomprthresh;


        if (params.toneCurve.autoexp) {
            LUTu aehist;
            int aehistcompr;
            imgsrc->getAutoExpHistogram (aehist, aehistcompr);
            ipf.getAutoExp (aehist, aehistcompr, params.toneCurve.clip, expcomp, bright, contr, black, hlcompr, hlcomprthresh);
        }
        if (params.toneCurve.histmatching) {
            if (!params.toneCurve.fromHistMatching) {
                imgsrc->getAutoMatchedToneCurve(params.icm, params.toneCurve.curve);
            }

            if (params.toneCurve.autoexp) {
                params.toneCurve.expcomp = 0.0;
            }

            params.toneCurve.autoexp = false;
            params.toneCurve.curveMode = ToneCurveParams::TcMode::FILMLIKE;
            params.toneCurve.curve2 = { 0 };
            params.toneCurve.brightness = 0;
            params.toneCurve.contrast = 0;
            params.toneCurve.black = 0;
        }
        if (params.logenc.enabled && params.logenc.autocompute) {
            ipf.getAutoLog(imgsrc, params.logenc);
        }

        // at this stage, we can flush the raw data to free up quite an important amount of memory
        // commented out because it makes the application crash when batch processing...
        // TODO: find a better place to flush rawData and rawRGB
        if (flush) {
            imgsrc->flushRawData();
            imgsrc->flushRGB();
        }

        return true;
    }

    void stage_denoise()
    {
        procparams::ProcParams& params = job->pparams;
        ImProcFunctions &ipf = *(ipf_p.get());

        if (params.denoise.enabled) {
            ipf.denoise(2, imgsrc, currWB, baseImg, dnstore, params.denoise);
        }
    }

    void stage_transform()
    {
        procparams::ProcParams& params = job->pparams;
        //ImProcFunctions ipf (&params, true);
        ImProcFunctions &ipf = * (ipf_p.get());

        imgsrc->convertColorSpace (baseImg, params.icm, currWB);

        // perform first analysis
        hist16 (65536);

        ipf.firstAnalysis (baseImg, params, hist16);

        ipf.dehaze(baseImg);
        ipf.dynamicRangeCompression(baseImg);

        // perform transform (excepted resizing)
        if (ipf.needsTransform()) {
            Imagefloat* trImg = nullptr;
            if (ipf.needsLuminanceOnly()) {
                trImg = baseImg;
            } else {
                trImg = new Imagefloat (fw, fh);
            }
            ipf.transform (baseImg, trImg, 0, 0, 0, 0, fw, fh, fw, fh,
                           imgsrc->getMetaData(), imgsrc->getRotateDegree(), true);
            if(trImg != baseImg) {
                delete baseImg;
                baseImg = trImg;
            }
        }
    }

    Imagefloat *stage_finish()
    {
        procparams::ProcParams& params = job->pparams;
        //ImProcFunctions ipf (&params, true);
        ImProcFunctions &ipf = * (ipf_p.get());

        if (params.dirpyrequalizer.cbdlMethod == "bef" && params.dirpyrequalizer.enabled && !params.colorappearance.enabled) {
            const int W = baseImg->getWidth();
            const int H = baseImg->getHeight();
            LabImage labcbdl (W, H);
            ipf.rgb2lab (*baseImg, labcbdl, params.icm.workingProfile);
            ipf.dirpyrequalizer (&labcbdl, 1);
            ipf.lab2rgb (labcbdl, *baseImg, params.icm.workingProfile);
        }

        //gamma TRC working
        if (params.icm.workingTRC == "Custom") { //exec TRC IN free
            const Glib::ustring profile = params.icm.workingProfile;

            if (profile == "sRGB" || profile == "Adobe RGB" || profile == "ProPhoto" || profile == "WideGamut" || profile == "BruceRGB" || profile == "Beta RGB" || profile == "BestRGB" || profile == "Rec2020" || profile == "ACESp0" || profile == "ACESp1") {
                const int cw = baseImg->getWidth();
                const int ch = baseImg->getHeight();
                cmsHTRANSFORM dummy = nullptr;
                // put gamma TRC to 1
                ipf.workingtrc(baseImg, baseImg, cw, ch, -5, params.icm.workingProfile, 2.4, 12.92310, dummy, true, false, false);
                //adjust TRC
                ipf.workingtrc(baseImg, baseImg, cw, ch, 5, params.icm.workingProfile, params.icm.workingTRCGamma, params.icm.workingTRCSlope, dummy, false, true, false);
            }
        }

        // RGB processing

        curve1 (65536);
        curve2 (65536);
        curve (65536, 0);
        satcurve (65536, 0);
        lhskcurve (65536, 0);
        lumacurve (32770, 0); // lumacurve[32768] and lumacurve[32769] will be set to 32768 and 32769 later to allow linear interpolation
        clcurve (65536, 0);
        wavclCurve (65536, 0);

        //if(params.blackwhite.enabled) params.toneCurve.hrenabled=false;

        CurveFactory::complexCurve (expcomp, black / 65535.0, hlcompr, hlcomprthresh, params.toneCurve.shcompr, bright, contr,
                                    params.toneCurve.curve, params.toneCurve.curve2,
                                    hist16, curve1, curve2, curve, dummy, customToneCurve1, customToneCurve2 );

        CurveFactory::RGBCurve (params.rgbCurves.rcurve, rCurve, 1);
        CurveFactory::RGBCurve (params.rgbCurves.gcurve, gCurve, 1);
        CurveFactory::RGBCurve (params.rgbCurves.bcurve, bCurve, 1);

        bool opautili = false;

        if (params.colorToning.enabled) {
            TMatrix wprof = ICCStore::getInstance()->workingSpaceMatrix (params.icm.workingProfile);
            double wp[3][3] = {
                {wprof[0][0], wprof[0][1], wprof[0][2]},
                {wprof[1][0], wprof[1][1], wprof[1][2]},
                {wprof[2][0], wprof[2][1], wprof[2][2]}
            };
            params.colorToning.getCurves (ctColorCurve, ctOpacityCurve, wp, opautili);
            clToningcurve (65536, 0);
            CurveFactory::curveToning (params.colorToning.clcurve, clToningcurve, 1);
            cl2Toningcurve (65536, 0);
            CurveFactory::curveToning (params.colorToning.cl2curve, cl2Toningcurve, 1);
        }

        labView = new LabImage (fw, fh);

        if (params.blackwhite.enabled) {
            CurveFactory::curveBW (params.blackwhite.beforeCurve, params.blackwhite.afterCurve, hist16, dummy, customToneCurvebw1, customToneCurvebw2, 1);
        }

        double rrm, ggm, bbm;
        float autor, autog, autob;
        float satLimit = float (params.colorToning.satProtectionThreshold) / 100.f * 0.7f + 0.3f;
        float satLimitOpacity = 1.f - (float (params.colorToning.saturatedOpacity) / 100.f);

        if (params.colorToning.enabled  && params.colorToning.autosat && params.colorToning.method != "LabGrid") { //for colortoning evaluation of saturation settings
            float moyS = 0.f;
            float eqty = 0.f;
            ipf.moyeqt (baseImg, moyS, eqty);//return image : mean saturation and standard dev of saturation
            float satp = ((moyS + 1.5f * eqty) - 0.3f) / 0.7f; //1.5 sigma ==> 93% pixels with high saturation -0.3 / 0.7 convert to Hombre scale

            if (satp >= 0.92f) {
                satp = 0.92f;    //avoid values too high (out of gamut)
            }

            if (satp <= 0.15f) {
                satp = 0.15f;    //avoid too low values
            }

            satLimit = 100.f * satp;

            satLimitOpacity = 100.f * (moyS - 0.85f * eqty); //-0.85 sigma==>20% pixels with low saturation
        }

        autor = -9000.f; // This will ask to compute the "auto" values for the B&W tool (have to be inferior to -5000)
        DCPProfile::ApplyState as;
        DCPProfile *dcpProf = imgsrc->getDCP (params.icm, as);

        LUTu histToneCurve;

        ipf.setDCPProfile(dcpProf, as);
        ipf.rgbProc (baseImg, labView, nullptr, curve1, curve2, curve, params.toneCurve.saturation, rCurve, gCurve, bCurve, satLimit, satLimitOpacity, ctColorCurve, ctOpacityCurve, opautili, clToningcurve, cl2Toningcurve, customToneCurve1, customToneCurve2, customToneCurvebw1, customToneCurvebw2, rrm, ggm, bbm, autor, autog, autob, expcomp, hlcompr, hlcomprthresh, histToneCurve);

        if (settings->verbose) {
            printf ("Output image / Auto B&W coefs:   R=%.2f   G=%.2f   B=%.2f\n", autor, autog, autob);
        }

        // if clut was used and size of clut cache == 1 we free the memory used by the clutstore (default clut cache size = 1 for 32 bit OS)
        if ( params.filmSimulation.enabled && !params.filmSimulation.clutFilename.empty() && options.clutCacheSize == 1) {
            CLUTStore::getInstance().clearCache();
        }

        // freeing up some memory
        customToneCurve1.Reset();
        customToneCurve2.Reset();
        ctColorCurve.Reset();
        ctOpacityCurve.Reset();
        customToneCurvebw1.Reset();
        customToneCurvebw2.Reset();

        // Freeing baseImg because not used anymore
        delete baseImg;
        baseImg = nullptr;

        if (pl) {
            pl->setProgress (0.55);
        }

        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        // start tile processing...???


        if (params.labCurve.contrast != 0) { //only use hist16 for contrast
            hist16.clear();

#ifdef _OPENMP
            #pragma omp parallel
#endif
            {
                LUTu hist16thr (hist16.getSize());  // one temporary lookup table per thread
                hist16thr.clear();
#ifdef _OPENMP
                #pragma omp for schedule(static) nowait
#endif

                for (int i = 0; i < fh; i++)
                    for (int j = 0; j < fw; j++) {
                        hist16thr[ (int) ((labView->L[i][j]))]++;
                    }

                #pragma omp critical
                {
                    hist16 += hist16thr;
                }
            }
        }

        bool utili;
        CurveFactory::complexLCurve (params.labCurve.brightness, params.labCurve.contrast, params.labCurve.lcurve, hist16, lumacurve, dummy, 1, utili);

        bool clcutili;
        CurveFactory::curveCL (clcutili, params.labCurve.clcurve, clcurve, 1);

        bool ccutili, cclutili;
        CurveFactory::complexsgnCurve (autili, butili, ccutili, cclutili, params.labCurve.acurve, params.labCurve.bcurve, params.labCurve.cccurve,
                                       params.labCurve.lccurve, curve1, curve2, satcurve, lhskcurve, 1);

        ipf.chromiLuminanceCurve (nullptr, 1, labView, labView, curve1, curve2, satcurve, lhskcurve, clcurve, lumacurve, utili, autili, butili, ccutili, cclutili, clcutili, dummy, dummy);

        ipf.vibrance (labView);
        ipf.labColorCorrectionRegions(labView, oX, oY, oW, oH);
        ipf.logEncoding(labView);

        if ((params.colorappearance.enabled && !params.colorappearance.tonecie) || (!params.colorappearance.enabled)) {
            ipf.EPDToneMap (labView, 5, 1);
        }

        if ((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)) {
            ipf.impulsedenoise (labView);
        }

        // for all treatments Defringe, Sharpening, Contrast detail ,Microcontrast they are activated if "CIECAM" function are disabled

        if ((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)) {
            ipf.defringe (labView);
        }

        if (params.sharpenEdge.enabled) {
            ipf.MLsharpen (labView);
        }

        if (params.sharpenMicro.enabled) {
            if ((params.colorappearance.enabled && !settings->autocielab) ||  (!params.colorappearance.enabled)) {
                ipf.MLmicrocontrast (labView);    //!params.colorappearance.sharpcie
            }
        }

        if (((params.colorappearance.enabled && !settings->autocielab) || (!params.colorappearance.enabled)) && params.sharpening.enabled) {
            ipf.sharpening (labView, params.sharpening);

        }

        WaveletParams WaveParams = params.wavelet;
        WavCurve wavCLVCurve;
        WavOpacityCurveRG waOpacityCurveRG;
        WavOpacityCurveBY waOpacityCurveBY;
        WavOpacityCurveW waOpacityCurveW;
        WavOpacityCurveWL waOpacityCurveWL;

        params.wavelet.getCurves (wavCLVCurve, waOpacityCurveRG, waOpacityCurveBY, waOpacityCurveW, waOpacityCurveWL );


        // directional pyramid wavelet
        if (params.dirpyrequalizer.cbdlMethod == "aft") {
            if ((params.colorappearance.enabled && !settings->autocielab)  || !params.colorappearance.enabled) {
                ipf.dirpyrequalizer (labView, 1);    //TODO: this is the luminance tonecurve, not the RGB one
            }
        }

        bool wavcontlutili = false;

        CurveFactory::curveWavContL (wavcontlutili, params.wavelet.wavclCurve, wavclCurve,/* hist16C, dummy,*/ 1);

        if (params.wavelet.enabled) {
            ipf.ip_wavelet (labView, labView, 2, WaveParams, wavCLVCurve, waOpacityCurveRG, waOpacityCurveBY, waOpacityCurveW,  waOpacityCurveWL, wavclCurve, 1);
        }

        wavCLVCurve.Reset();

        ipf.softLight(labView);
        ipf.localContrast(labView);

        //Colorappearance and tone-mapping associated

        int f_w = 1, f_h = 1;

        if (params.colorappearance.tonecie || params.colorappearance.enabled) {
            f_w = fw;
            f_h = fh;
        }

        CieImage *cieView = new CieImage (f_w, (f_h));

        CurveFactory::curveLightBrightColor (
            params.colorappearance.curve,
            params.colorappearance.curve2,
            params.colorappearance.curve3,
            hist16, dummy,
            dummy, dummy,
            customColCurve1,
            customColCurve2,
            customColCurve3,
            1);

        if (params.colorappearance.enabled) {
            double adap;
            int imgNum = 0;
            if (imgsrc->getSensorType() == ST_BAYER) {
                imgNum = params.raw.bayersensor.imageNum;
            } else if (imgsrc->getSensorType() == ST_FUJI_XTRANS) {
                //imgNum = params.raw.xtranssensor.imageNum;
            }
            float fnum = imgsrc->getMetaData()->getFNumber (imgNum);         // F number
            float fiso = imgsrc->getMetaData()->getISOSpeed (imgNum) ;       // ISO
            float fspeed = imgsrc->getMetaData()->getShutterSpeed (imgNum) ; //speed
            float fcomp = imgsrc->getMetaData()->getExpComp (imgNum);        //compensation + -

            if (fnum < 0.3f || fiso < 5.f || fspeed < 0.00001f) {
                adap = 2000.;
            }//if no exif data or wrong
            else {
                float E_V = fcomp + log2 ((fnum * fnum) / fspeed / (fiso / 100.f));
                E_V += params.toneCurve.expcomp;// exposure compensation in tonecurve ==> direct EV
                E_V += log2 (params.raw.expos); // exposure raw white point ; log2 ==> linear to EV
                adap = powf (2.f, E_V - 3.f); //cd / m2
            }

            LUTf CAMBrightCurveJ;
            LUTf CAMBrightCurveQ;
            float CAMMean = NAN;

            float d, dj, yb;
            ipf.ciecam_02float (cieView, float (adap), 1, 2, labView, &params, customColCurve1, customColCurve2, customColCurve3, dummy, dummy, CAMBrightCurveJ, CAMBrightCurveQ, CAMMean, 5, 1, true, d, dj, yb, 1);
        }

        delete cieView;
        cieView = nullptr;




        // end tile processing...???
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if (pl) {
            pl->setProgress (0.60);
        }

        int imw, imh;
        double tmpScale = ipf.resizeScale (&params, fw, fh, imw, imh);
        bool labResize = params.resize.enabled && params.resize.method != "Nearest" && (tmpScale != 1.0 || params.prsharpening.enabled);
        LabImage *tmplab;

        // crop and convert to rgb16
        int cx = 0, cy = 0, cw = labView->W, ch = labView->H;

        if (params.crop.enabled) {
            cx = params.crop.x;
            cy = params.crop.y;
            cw = params.crop.w;
            ch = params.crop.h;

            if (labResize) { // crop lab data
                tmplab = new LabImage (cw, ch);

                for (int row = 0; row < ch; row++) {
                    for (int col = 0; col < cw; col++) {
                        tmplab->L[row][col] = labView->L[row + cy][col + cx];
                        tmplab->a[row][col] = labView->a[row + cy][col + cx];
                        tmplab->b[row][col] = labView->b[row + cy][col + cx];
                    }
                }

                delete labView;
                labView = tmplab;
                cx = 0;
                cy = 0;
            }
        }

        if (labResize) { // resize lab data
            if ((labView->W != imw || labView->H != imh) &&
                (params.resize.allowUpscaling || (labView->W >= imw && labView->H >= imh))) {
                // resize image
                tmplab = new LabImage (imw, imh);
                ipf.Lanczos (labView, tmplab, tmpScale);
                delete labView;
                labView = tmplab;
            }
            cw = labView->W;
            ch = labView->H;

            if (params.prsharpening.enabled) {
                for (int i = 0; i < ch; i++) {
                    for (int j = 0; j < cw; j++) {
                        labView->L[i][j] = labView->L[i][j] < 0.f ? 0.f : labView->L[i][j];
                    }
                }
                ipf.sharpening (labView, params.prsharpening);
            }
        }

        cmsHPROFILE jprof = nullptr;
        constexpr bool customGamma = false;
        constexpr bool useLCMS = false;
        bool bwonly = params.blackwhite.enabled && !params.colorToning.enabled && !autili && !butili ;

        ///////////// Custom output gamma has been removed, the user now has to create
        ///////////// a new output profile with the ICCProfileCreator

        // if Default gamma mode: we use the profile selected in the "Output profile" combobox;
        // gamma come from the selected profile, otherwise it comes from "Free gamma" tool

        Imagefloat* readyImg = ipf.lab2rgbOut (labView, cx, cy, cw, ch, params.icm);

        if (settings->verbose) {
            printf ("Output profile_: \"%s\"\n", params.icm.outputProfile.c_str());
        }

        delete labView;
        labView = nullptr;



        if (bwonly) { //force BW r=g=b
            if (settings->verbose) {
                printf ("Force BW\n");
            }

            for (int ccw = 0; ccw < cw; ccw++) {
                for (int cch = 0; cch < ch; cch++) {
                    readyImg->r (cch, ccw) = readyImg->g (cch, ccw);
                    readyImg->b (cch, ccw) = readyImg->g (cch, ccw);
                }
            }
        }

        if (pl) {
            pl->setProgress (0.70);
        }

        if (tmpScale != 1.0 && params.resize.method == "Nearest" &&
            (params.resize.allowUpscaling || (readyImg->getWidth() >= imw && readyImg->getHeight() >= imh))) { // resize rgb data (gamma applied)
            Imagefloat* tempImage = new Imagefloat (imw, imh);
            ipf.resize (readyImg, tempImage, tmpScale);
            delete readyImg;
            readyImg = tempImage;
        }

        switch (params.metadata.mode) {
        case MetaDataParams::TUNNEL:
            // Sending back the whole first root, which won't necessarily be the selected frame number
            // and may contain subframe depending on initial raw's hierarchy
            readyImg->setMetadata (ii->getMetaData()->getRootExifData ());
            break;
        case MetaDataParams::EDIT:
            // ask for the correct frame number, but may contain subframe depending on initial raw's hierarchy
            readyImg->setMetadata (ii->getMetaData()->getBestExifData(imgsrc, &params.raw), params.exif, params.iptc);
            break;
        default: // case MetaDataParams::STRIP
            // nothing to do
            break;
        }


        // Setting the output curve to readyImg
        if (customGamma) {
            if (!useLCMS) {
                // use corrected sRGB profile in order to apply a good TRC if present, otherwise use LCMS2 profile generated by lab2rgb16 w/ gamma
                ProfileContent pc (jprof);
                readyImg->setOutputProfile (pc.getData().c_str(), pc.getData().size());
            }
        } else {
            // use the selected output profile if present, otherwise use LCMS2 profile generate by lab2rgb16 w/ gamma

            if (params.icm.outputProfile != "" && params.icm.outputProfile != ColorManagementParams::NoICMString) {

                // if ICCStore::getInstance()->getProfile send back an object, then ICCStore::getInstance()->getContent will do too
                cmsHPROFILE jprof = ICCStore::getInstance()->getProfile (params.icm.outputProfile); //get outProfile

                if (jprof == nullptr) {
                    if (settings->verbose) {
                        printf ("\"%s\" ICC output profile not found!\n - use LCMS2 substitution\n", params.icm.outputProfile.c_str());
                    }
                } else {
                    if (settings->verbose) {
                        printf ("Using \"%s\" output profile\n", params.icm.outputProfile.c_str());
                    }

                    ProfileContent pc = ICCStore::getInstance()->getContent (params.icm.outputProfile);
                    readyImg->setOutputProfile (pc.getData().c_str(), pc.getData().size());
                }
            } else {
                // No ICM
                readyImg->setOutputProfile (nullptr, 0);
            }
        }

//    t2.set();
//    if( settings->verbose )
//           printf("Total:- %d usec\n", t2.etime(t1));

        if (!job->initialImage) {
            ii->decreaseRef ();
        }

        delete job;

        if (pl) {
            pl->setProgress (0.75);
        }

        /*  curve1.reset();curve2.reset();
            curve.reset();
            satcurve.reset();
            lhskcurve.reset();

            rCurve.reset();
            gCurve.reset();
            bCurve.reset();
            hist16.reset();
            hist16C.reset();
        */
        return readyImg;
    }

    void stage_early_resize()
    {
        procparams::ProcParams& params = job->pparams;
        ImProcFunctions &ipf = * (ipf_p.get());

        int imw, imh;
        double scale_factor = ipf.resizeScale (&params, fw, fh, imw, imh);

        if (params.crop.enabled) {
            int cx = params.crop.x;
            int cy = params.crop.y;
            int cw = params.crop.w;
            int ch = params.crop.h;
            oX = cx * scale_factor;
            oY = cy * scale_factor;

            Imagefloat *cropped = new Imagefloat(cw, ch);

            for (int row = 0; row < ch; row++) {
                for (int col = 0; col < cw; col++) {
                    cropped->r(row, col) = baseImg->r(row + cy, col + cx);
                    cropped->g(row, col) = baseImg->g(row + cy, col + cx);
                    cropped->b(row, col) = baseImg->b(row + cy, col + cx);
                }
            }

            delete baseImg;
            baseImg = cropped;
        }

        assert (params.resize.enabled);

        // resize image
        if (params.resize.allowUpscaling || (imw <= fw && imh <= fh)) {
            Imagefloat *resized = new Imagefloat(imw, imh);
            ipf.Lanczos(baseImg, resized, scale_factor);
            delete baseImg;
            baseImg = resized;
        }

//        adjust_procparams (scale_factor);
        params.resize.enabled = false;
        params.crop.enabled = false;

        oW *= scale_factor;
        oH *= scale_factor;

        fw = imw;
        fh = imh;
    }

    void adjust_procparams (double scale_factor)
    {
        procparams::ProcParams &params = job->pparams;
        procparams::ProcParams defaultparams;

        if (!params.sharpening.enabled) {
            params.sharpening = params.prsharpening;
        }
            
        ImProcFunctions &ipf = *(ipf_p.get());
        ipf.setScale(1.0 / scale_factor);

        params.wavelet.strength *= scale_factor;

        if (params.raw.xtranssensor.method == procparams::RAWParams::XTransSensor::getMethodString(procparams::RAWParams::XTransSensor::Method::THREE_PASS)) {
            params.raw.xtranssensor.method = procparams::RAWParams::XTransSensor::getMethodString(procparams::RAWParams::XTransSensor::Method::ONE_PASS);
        }

        if (params.raw.bayersensor.method == procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::PIXELSHIFT)) {
            params.raw.bayersensor.method = procparams::RAWParams::BayerSensor::getMethodString(procparams::RAWParams::BayerSensor::Method::AMAZE);
        }
    }

private:
    ProcessingJobImpl* job;
    int& errorCode;
    ProgressListener* pl;
    bool flush;

    // internal state
    std::unique_ptr<ImProcFunctions> ipf_p;
    InitialImage *ii;
    ImageSource *imgsrc;
    int fw;
    int fh;
    int oX; // parameters for the area masks 
    int oY; // in labColorCorrectionRegions
    int oW; //
    int oH; //

    int tr;
    PreviewProps pp;

    int tilesize;
    int overlap;

    ImProcFunctions::DenoiseInfoStore dnstore;

    double expcomp;
    int bright;
    int contr;
    int black;
    int hlcompr;
    int hlcomprthresh;

    ColorTemp currWB;
    Imagefloat *baseImg;
    LabImage* labView;

    LUTu hist16;

    LUTf curve1;
    LUTf curve2;
    LUTf curve;
    LUTf satcurve;
    LUTf lhskcurve;
    LUTf lumacurve;
    LUTf clcurve;
    LUTf clToningcurve;
    LUTf cl2Toningcurve;
    LUTf wavclCurve;

    LUTf rCurve;
    LUTf gCurve;
    LUTf bCurve;
    LUTu dummy;

    ToneCurve customToneCurve1, customToneCurve2;
    ColorGradientCurve ctColorCurve;
    OpacityCurve ctOpacityCurve;
    ColorAppearance customColCurve1, customColCurve2, customColCurve3 ;
    ToneCurve customToneCurvebw1;
    ToneCurve customToneCurvebw2;

    bool autili, butili;
};

} // namespace


IImagefloat* processImage (ProcessingJob* pjob, int& errorCode, ProgressListener* pl, bool flush)
{
    ImageProcessor proc (pjob, errorCode, pl, flush);
    return proc();
}

void batchProcessingThread (ProcessingJob* job, BatchProcessingListener* bpl)
{

    ProcessingJob* currentJob = job;

    while (currentJob) {
        int errorCode;
        IImagefloat* img = processImage (currentJob, errorCode, bpl, true);

        if (errorCode) {
            bpl->error (M ("MAIN_MSG_CANNOTLOAD"));
            currentJob = nullptr;
        } else {
            try {
                currentJob = bpl->imageReady (img);
            } catch (Glib::Exception& ex) {
                bpl->error (ex.what());
                currentJob = nullptr;
            }
        }
    }
}

void startBatchProcessing (ProcessingJob* job, BatchProcessingListener* bpl)
{

    if (bpl) {
        Glib::Thread::create (sigc::bind (sigc::ptr_fun (batchProcessingThread), job, bpl), 0, true, true, Glib::THREAD_PRIORITY_LOW);
    }

}

}
