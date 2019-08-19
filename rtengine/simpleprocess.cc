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
#include "metadata.h"
#undef THREAD_PRIORITY_NORMAL

namespace rtengine
{
extern const Settings* settings;

namespace {

class ImageProcessor {
public:
    ImageProcessor(ProcessingJob* pjob,int &errorCode,ProgressListener *pl,
                   bool flush):
        job(static_cast<ProcessingJobImpl*>(pjob)),
        errorCode(errorCode),
        pl(pl),
        flush(flush),
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
        labView(nullptr),
        pipeline_scale(1.0)
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
        return stage_finish(false);
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
        return stage_finish(true);
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
        } else if (imgsrc->getSensorType() == ST_FUJI_XTRANS) {
            imgsrc->setBorder(params.raw.xtranssensor.border);
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
        imgsrc->demosaic(params.raw, autoContrast, contrastThreshold);


        if (pl) {
            pl->setProgress (0.30);
        }
        pp = PreviewProps (0, 0, fw, fh, 1);

        if (pl) {
            pl->setProgress (0.40);
        }

        imgsrc->HLRecovery_Global(params.exposure);


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
        imgsrc->getImage (currWB, tr, baseImg, pp, params.exposure, params.raw);

        if (pl) {
            pl->setProgress (0.50);
        }

        if (params.exposure.enabled && params.exposure.autoexp) {
            LUTu aehist;
            int aehistcompr;
            imgsrc->getAutoExpHistogram (aehist, aehistcompr);
            params.brightContrSat.enabled = true;
            ipf.getAutoExp (aehist, aehistcompr, params.exposure.clip, params.exposure.expcomp, params.brightContrSat.brightness, params.brightContrSat.contrast, params.exposure.black, params.exposure.hlcompr, params.exposure.hlcomprthresh);
        }
        if (params.toneCurve.histmatching) {
            if (!params.toneCurve.fromHistMatching) {
                imgsrc->getAutoMatchedToneCurve(params.icm, params.toneCurve.curve);
            }
        }
        if (params.logenc.enabled && params.logenc.autocompute) {
            ipf.getAutoLog(imgsrc, params.logenc);
        }

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
            ipf.denoise(imgsrc, currWB, baseImg, dnstore, params.denoise);
        }
    }

    void stage_transform()
    {
        procparams::ProcParams &params = job->pparams;
        ImProcFunctions &ipf = *(ipf_p.get());

        imgsrc->convertColorSpace(baseImg, params.icm, currWB);

        LUTu hist16(65536);
        ipf.firstAnalysis(baseImg, params, hist16);

        ipf.dehaze(baseImg);
        ipf.dynamicRangeCompression(baseImg);

        // perform transform (excepted resizing)
        if (ipf.needsTransform()) {
            Imagefloat *trImg = nullptr;
            if (ipf.needsLuminanceOnly()) {
                trImg = baseImg;
            } else {
                trImg = new Imagefloat (fw, fh);
            }
            ipf.transform(baseImg, trImg, 0, 0, 0, 0, fw, fh, fw, fh,
                          imgsrc->getMetaData(), imgsrc->getRotateDegree(), true);
            if (trImg != baseImg) {
                delete baseImg;
                baseImg = trImg;
            }
        }
    }

    Imagefloat *stage_finish(bool is_fast)
    {
        procparams::ProcParams& params = job->pparams;
        ImProcFunctions &ipf = * (ipf_p.get());

        // RGB processing
        labView = new LabImage (fw, fh);

        DCPProfile::ApplyState as;
        DCPProfile *dcpProf = imgsrc->getDCP (params.icm, as);

        LUTu histToneCurve;

        ipf.setDCPProfile(dcpProf, as);
        ipf.rgbProc(baseImg, labView);

        // if clut was used and size of clut cache == 1 we free the memory used by the clutstore (default clut cache size = 1 for 32 bit OS)
        if ( params.filmSimulation.enabled && !params.filmSimulation.clutFilename.empty() && options.clutCacheSize == 1) {
            CLUTStore::getInstance().clearCache();
        }

        // Freeing baseImg because not used anymore
        delete baseImg;
        baseImg = nullptr;

        if (pl) {
            pl->setProgress (0.55);
        }

        ipf.sharpening(labView, params.sharpening);
        bool stop = ipf.colorCorrection(labView, oX, oY, oW, oH);
        stop = stop || ipf.guidedSmoothing(labView, oX, oY, oW, oH);
        stop = stop || ipf.contrastByDetailLevels(labView, oX, oY, oW, oH);
        if (!stop) {
            ipf.localContrast(labView);
            ipf.logEncoding(labView);
            ipf.labAdjustments(labView);
        }

        stop = stop || ipf.toneMapping(labView, oX, oY, oW, oH);
        if (!stop) {
            ipf.impulsedenoise(labView);
            ipf.defringe(labView);
            ipf.MLmicrocontrast(labView);
            //ipf.sharpening(labView, params.sharpening);
        }
        //stop = stop || ipf.contrastByDetailLevels(labView, oX, oY, oW, oH);
        if (!stop) {
            ipf.softLight(labView);
            // ipf.localContrast(labView);
            ipf.filmGrain(labView);
        }

        if (pl) {
            pl->setProgress (0.60);
        }

        int imw, imh;
        double tmpScale = ipf.resizeScale(&params, fw, fh, imw, imh);
        bool labResize = params.resize.enabled && params.resize.method != "Nearest" && (tmpScale != 1.0 || params.prsharpening.enabled);
        LabImage *tmplab;

        // crop and convert to rgb16
        int cx = 0, cy = 0, cw = labView->W, ch = labView->H;

        if (params.crop.enabled && !is_fast) {
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

        if (labResize && !is_fast) { // resize lab data
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
        }

        if (labResize || is_fast) {
            if (params.prsharpening.enabled) {
                for (int i = 0; i < ch; i++) {
                    for (int j = 0; j < cw; j++) {
                        labView->L[i][j] = labView->L[i][j] < 0.f ? 0.f : labView->L[i][j];
                    }
                }
                ipf.sharpening(labView, params.prsharpening);
            }
        }

        cmsHPROFILE jprof = nullptr;
        constexpr bool customGamma = false;
        constexpr bool useLCMS = false;
        bool bwonly = params.blackwhite.enabled;

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

        Exiv2Metadata info(imgsrc->getFileName());
        switch (params.metadata.mode) {
        case MetaDataParams::TUNNEL:
            readyImg->setMetadata(info);
            break;
        case MetaDataParams::EDIT:
            info.setExif(params.exif);
            info.setIptc(params.iptc);
            readyImg->setMetadata(info);
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

        if (!job->initialImage) {
            ii->decreaseRef ();
        }

        delete job;

        if (pl) {
            pl->setProgress (0.75);
        }

        return readyImg;
    }

    void stage_early_resize()
    {
        procparams::ProcParams& params = job->pparams;
        ImProcFunctions &ipf = * (ipf_p.get());

        int imw, imh;
        double scale_factor = ipf.resizeScale (&params, fw, fh, imw, imh);

        if (scale_factor == 1.f) {
            return;
        }

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

        ImProcFunctions &ipf = *(ipf_p.get());
        pipeline_scale = 1.0 / scale_factor;
        ipf.setScale(pipeline_scale);

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
    int oY; // in colorCorrection
    int oW; //
    int oH; //

    int tr;
    PreviewProps pp;

    int tilesize;
    int overlap;

    ImProcFunctions::DenoiseInfoStore dnstore;

    ColorTemp currWB;
    Imagefloat *baseImg;
    LabImage* labView;

    double pipeline_scale;
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
