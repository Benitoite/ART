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
#include "improccoordinator.h"
#include "curves.h"
#include "mytime.h"
#include "refreshmap.h"
#include "../rtgui/ppversion.h"
#include "colortemp.h"
#include "improcfun.h"
#include "iccstore.h"
#include <iostream>
#include <fstream>
#include <string>
#include "color.h"
#include "metadata.h"
#include "perspectivecorrection.h"
#ifdef _OPENMP
#include <omp.h>
#endif

namespace rtengine {

extern const Settings* settings;

ImProcCoordinator::ImProcCoordinator():
    orig_prev(nullptr),
    oprevi(nullptr),
    drcomp_11_dcrop_cache(nullptr),
    previmg(nullptr),
    workimg(nullptr),
    imgsrc(nullptr),
    lastAwbEqual(0.),
    ipf(&params, true),
    monitorIntent(RI_RELATIVE),
    softProof(false),
    gamutCheck(false),
    sharpMask(false),
    scale(10),
    highDetailPreprocessComputed(false),
    highDetailRawComputed(false),
    allocated(false), 

    vhist16(65536),
    histRed(256), histRedRaw(256),
    histGreen(256), histGreenRaw(256),
    histBlue(256), histBlueRaw(256),
    histLuma(256),
    histToneCurve(256),
    histLCurve(256),
    histCCurve(256),
    histLLCurve(256),
    histLCAM(256),
    histCCAM(256),
    histChroma(256),
    histLRETI(256),
    
    fw(0), fh(0), tr(0),
    fullw(1), fullh(1),
    pW(-1), pH(-1),
    plistener(nullptr),
    imageListener(nullptr),
    aeListener(nullptr),
    awbListener(nullptr),
    flatFieldAutoClipListener(nullptr),
    bayerAutoContrastListener(nullptr),
    xtransAutoContrastListener(nullptr),
    frameCountListener(nullptr),
    imageTypeListener(nullptr),
    adnListener(nullptr),
    hListener(nullptr),
    autoLogListener(nullptr),
    autoRadiusListener(nullptr),
    resultValid(false),
    lastOutputProfile("BADFOOD"),
    lastOutputIntent(RI__COUNT),
    lastOutputBPC(false),
    thread(nullptr),
    changeSinceLast(0),
    updaterRunning(false),
    destroying(false),
    highQualityComputed(false),
    customTransformIn(nullptr),
    customTransformOut(nullptr)
{
    for (int i = 0; i < 3; ++i) {
        bufs_[i] = nullptr;
    }
}

void ImProcCoordinator::assign(ImageSource* imgsrc)
{
    this->imgsrc = imgsrc;
}

ImProcCoordinator::~ImProcCoordinator()
{

    destroying = true;
    updaterThreadStart.lock();

    if (updaterRunning && thread) {
        thread->join();
    }

    mProcessing.lock();
    mProcessing.unlock();
    freeAll();

    if (drcomp_11_dcrop_cache) {
        delete drcomp_11_dcrop_cache;
        drcomp_11_dcrop_cache = nullptr;
    }

    std::vector<Crop*> toDel = crops;

    for (size_t i = 0; i < toDel.size(); i++) {
        delete toDel[i];
    }

    imgsrc->decreaseRef();

    if(customTransformIn) {
        cmsDeleteTransform(customTransformIn);
        customTransformIn = nullptr;
    }

    if(customTransformOut) {
        cmsDeleteTransform(customTransformOut);
        customTransformOut = nullptr;
    }

    updaterThreadStart.unlock();
}

DetailedCrop* ImProcCoordinator::createCrop(::EditDataProvider *editDataProvider, bool isDetailWindow)
{

    return new Crop(this, editDataProvider, isDetailWindow);
}


// todo: bitmask containing desired actions, taken from changesSinceLast
// cropCall: calling crop, used to prevent self-updates  ...doesn't seem to be used
void ImProcCoordinator::updatePreviewImage(int todo, bool panningRelatedChange)
{
    MyMutex::MyLock processingLock(mProcessing);
    int numofphases = 14;
    int readyphase = 0;

    DCPProfile *dcpProf = imgsrc->getDCP(params.icm, dcpApplyState);
    ipf.setDCPProfile(dcpProf, dcpApplyState);
    ipf.setViewport(0, 0, -1, -1);
    ipf.setOutputHistograms(&histToneCurve, &histCCurve, &histLCurve);

    if (todo == CROP && ipf.needsPCVignetting()) {
        todo |= M_LUMINANCE; //TRANSFORM;    // Change about Crop does affect TRANSFORM
    }

    bool highDetailNeeded = false;

    if (options.prevdemo == PD_Sidecar) {
        highDetailNeeded = true;    //i#2664
    } else {
        highDetailNeeded = (todo & M_HIGHQUAL);
    }

    // Check if any detail crops need high detail. If not, take a fast path short cut
    if (!highDetailNeeded) {
        for (size_t i = 0; i < crops.size(); i++)
            if (crops[i]->get_skip() == 1) {   // skip=1 -> full resolution
                highDetailNeeded = true;
                break;
            }
    }

    ipf.setPipetteBuffer(nullptr);
    bool stop = false;
                    
    if (((todo & ALL) == ALL) || (todo & M_MONITOR) || panningRelatedChange || (highDetailNeeded && options.prevdemo != PD_Sidecar)) {
        if (todo == CROP && ipf.needsPCVignetting()) {
            todo |= TRANSFORM;    // Change about Crop does affect TRANSFORM
        }
    
        RAWParams rp = params.raw;
        if (!highDetailNeeded) {
            // if below 100% magnification, take a fast path
            if (rp.bayersensor.method != RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::NONE) && rp.bayersensor.method != RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::MONO)) {
                rp.bayersensor.method = RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::FAST);
            }
    
            //bayerrp.all_enhance = false;
    
            if (rp.xtranssensor.method != RAWParams::XTransSensor::getMethodString(RAWParams::XTransSensor::Method::NONE) && rp.xtranssensor.method != RAWParams::XTransSensor::getMethodString(RAWParams::XTransSensor::Method::MONO)) {
                rp.xtranssensor.method = RAWParams::XTransSensor::getMethodString(RAWParams::XTransSensor::Method::FAST);
            }
    
            rp.bayersensor.ccSteps = 0;
            rp.xtranssensor.ccSteps = 0;
            //rp.deadPixelFilter = rp.hotPixelFilter = false;
        }
    
        progress("Applying white balance, color correction & sRGB conversion...", 100 * readyphase / numofphases);
    
        if (frameCountListener) {
            frameCountListener->FrameCountChanged(imgsrc->getFrameCount(), params.raw.bayersensor.imageNum);
        }
    
        // raw auto CA is bypassed if no high detail is needed, so we have to compute it when high detail is needed
        if ((todo & M_PREPROC) || (!highDetailPreprocessComputed && highDetailNeeded)) {
            imgsrc->setCurrentFrame(params.raw.bayersensor.imageNum);
    
            imgsrc->preprocess(rp, params.lensProf, params.coarse);
            if (flatFieldAutoClipListener && rp.ff_AutoClipControl) {
                flatFieldAutoClipListener->flatFieldAutoClipValueChanged(imgsrc->getFlatFieldAutoClipValue());
            }
            imgsrc->getRAWHistogram(histRedRaw, histGreenRaw, histBlueRaw);
    
            highDetailPreprocessComputed = highDetailNeeded;

            if (todo & M_RAW) {
                imgsrc->filmNegativeProcess(params.filmNegative);
            }                
        }
    
        /*
        Demosaic is kicked off only when
        Detail considerations:
            accurate detail is not displayed yet needed based on preview specifics (driven via highDetailNeeded flag)
        OR
        HLR considerations:
            Color HLR alters rgb output of demosaic, so re-demosaic is needed when Color HLR is being turned off;
            if HLR is enabled and changing method *from* Color to any other method
            OR HLR gets disabled when Color method was selected
        */
        // If high detail (=100%) is newly selected, do a demosaic update, since the last was just with FAST
    
        if (imageTypeListener) {
            imageTypeListener->imageTypeChanged(imgsrc->isRAW(), imgsrc->getSensorType() == ST_BAYER, imgsrc->getSensorType() == ST_FUJI_XTRANS, imgsrc->isMono());
        }

        const bool hrcolor = params.exposure.enabled && params.exposure.hrmode == procparams::ExposureParams::HR_COLOR;
        
        if ((todo & M_RAW) || (!highDetailRawComputed && highDetailNeeded) || (!hrcolor && imgsrc->isRGBSourceModified())) {
            if (settings->verbose) {
                if (imgsrc->getSensorType() == ST_BAYER) {
                    printf("Demosaic Bayer image n.%d using method: %s\n", rp.bayersensor.imageNum + 1, rp.bayersensor.method.c_str());
                } else if (imgsrc->getSensorType() == ST_FUJI_XTRANS) {
                    printf("Demosaic X-Trans image with using method: %s\n", rp.xtranssensor.method.c_str());
                }
            }
            if(imgsrc->getSensorType() == ST_BAYER) {
                if(params.raw.bayersensor.method != RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::PIXELSHIFT)) {
                    imgsrc->setBorder(params.raw.bayersensor.border);
                } else {
                    imgsrc->setBorder(std::max(params.raw.bayersensor.border, 2));
                }
            } else if (imgsrc->getSensorType() == ST_FUJI_XTRANS) {
                imgsrc->setBorder(params.raw.xtranssensor.border);
            }
            bool autoContrast = imgsrc->getSensorType() == ST_BAYER ? params.raw.bayersensor.dualDemosaicAutoContrast : params.raw.xtranssensor.dualDemosaicAutoContrast;
            double contrastThreshold = imgsrc->getSensorType() == ST_BAYER ? params.raw.bayersensor.dualDemosaicContrast : params.raw.xtranssensor.dualDemosaicContrast;
            imgsrc->demosaic(rp, autoContrast, contrastThreshold); //enabled demosaic

            if (imgsrc->getSensorType() == ST_BAYER && bayerAutoContrastListener && autoContrast) {
                bayerAutoContrastListener->autoContrastChanged(autoContrast ? contrastThreshold : -1.0);
            }
            if (imgsrc->getSensorType() == ST_FUJI_XTRANS && xtransAutoContrastListener && autoContrast) {
                xtransAutoContrastListener->autoContrastChanged(autoContrast ? contrastThreshold : -1.0);
            }
    
            // if a demosaic happened we should also call getimage later, so we need to set the M_INIT flag
            todo |= M_INIT;
    
            if (highDetailNeeded) {
                highDetailRawComputed = true;
            } else {
                highDetailRawComputed = false;
            }
        }   
    
        int tr = getCoarseBitMask(params.coarse);    
        imgsrc->getFullSize(fw, fh, tr);
        PreviewProps pp(0, 0, fw, fh, scale);
        ipf.setScale(scale);

        if (todo & (M_INIT | M_LINDENOISE | M_HDR)) {
            MyMutex::MyLock initLock(minit);  // Also used in crop window
    
            imgsrc->HLRecovery_Global(params.exposure);   // this handles Color HLRecovery
    
    
            if (settings->verbose) {
                printf("Applying white balance, color correction & sRBG conversion...\n");
            }
    
            currWB = ColorTemp(params.wb.temperature, params.wb.green, params.wb.equal, params.wb.method);
    
            if (!params.wb.enabled) {
                currWB = ColorTemp();
            } else if (params.wb.method == "Camera") {
                currWB = imgsrc->getWB();
            } else if (params.wb.method == "Auto") {
                if (lastAwbEqual != params.wb.equal) {
                    double rm, gm, bm;
                    imgsrc->getAutoWBMultipliers(rm, gm, bm);
    
                    if (rm != -1.) {
                        autoWB.update(rm, gm, bm, params.wb.equal);
                        lastAwbEqual = params.wb.equal;
                    } else {
                        lastAwbEqual = -1.;
                        autoWB.useDefaults(params.wb.equal);
                    }
                }
    
                currWB = autoWB;
            }
    
            if (params.wb.enabled) {
                params.wb.temperature = currWB.getTemp();
                params.wb.green = currWB.getGreen();
            }
    
            if (params.wb.method == "Auto" && awbListener && params.wb.enabled) {
                awbListener->WBChanged(params.wb.temperature, params.wb.green);
            }
    
            // int tr = getCoarseBitMask(params.coarse);    
            // imgsrc->getFullSize(fw, fh, tr);
            // PreviewProps pp(0, 0, fw, fh, scale);
            // ipf.setScale(scale);
            // Will (re)allocate the preview's buffers
            setScale(scale);
            imgsrc->getImage(currWB, tr, orig_prev, pp, params.exposure, params.raw);
            denoiseInfoStore.valid = false;
            imgsrc->convertColorSpace(orig_prev, params.icm, currWB);
    
            ipf.firstAnalysis(orig_prev, params, vhist16);
        }

        orig_prev->assignColorSpace(params.icm.workingProfile);
        readyphase++;
    
        if ((todo & M_HDR) && (params.fattal.enabled || params.dehaze.enabled)) {
            if (drcomp_11_dcrop_cache) {
                delete drcomp_11_dcrop_cache;
                drcomp_11_dcrop_cache = nullptr;
            }
    
            stop = ipf.process(ImProcFunctions::Pipeline::NAVIGATOR, ImProcFunctions::Stage::STAGE_0, orig_prev);
    
            if (oprevi != orig_prev) {
                delete oprevi;
            }
        }
    
        oprevi = orig_prev;
    
        progress("Rotate / Distortion...", 100 * readyphase / numofphases);
        // Remove transformation if unneeded
        if (ipf.needsTransform()) {
            assert(oprevi);
            Imagefloat *op = oprevi;
            oprevi = new Imagefloat(pW, pH, op);
            ipf.transform(op, oprevi, 0, 0, 0, 0, pW, pH, fw, fh,
                          imgsrc->getMetaData(), imgsrc->getRotateDegree(), false);
        }
    
        readyphase++;
        progress("Preparing shadow/highlight map...", 100 * readyphase / numofphases);
    
        readyphase++;

        if ((todo & M_INIT) && autoRadiusListener) {
            if (!imgsrc->getDeconvAutoRadius(nullptr)) {
                autoRadiusListener->autoDeconvRadiusChanged(-1);
            } else {
                autoRadiusListener->autoDeconvRadiusChanged(params.sharpening.deconvradius);
            }
        }
        
        if (todo & M_AUTOEXP) {
            if (params.toneCurve.histmatching) {
                if (!params.toneCurve.fromHistMatching) {
                    imgsrc->getAutoMatchedToneCurve(params.icm, params.toneCurve.curve);
                }
                params.toneCurve.fromHistMatching = true;
    
                if (aeListener) {
                    aeListener->autoMatchedToneCurveChanged(params.toneCurve.curveMode, params.toneCurve.curve);
                }
            }

            if (params.logenc.enabled && params.logenc.autocompute) {
                ipf.getAutoLog(imgsrc, params.logenc);
                if (autoLogListener) {
                    autoLogListener->logEncodingChanged(params.logenc);
                }
            }

            if (params.sharpening.deconvAutoRadius) {
                float r = -1.f;
                if (imgsrc->getDeconvAutoRadius(&r)) {
                    params.sharpening.deconvradius = r;
                } else {
                    r = -1.f;
                }
                if (autoRadiusListener) {
                    autoRadiusListener->autoDeconvRadiusChanged(r);
                }
            }
        }
    
        progress("Exposure curve & CIELAB conversion...", 100 * readyphase / numofphases);
    
        if ((todo & M_RGBCURVE) || (todo & M_CROP)) {
            // if it's just crop we just need the histogram, no image updates
            if (todo & M_RGBCURVE) {
                //initialize rrm bbm ggm different from zero to avoid black screen in some cases
                oprevi->copyTo(bufs_[0]);
                stop = stop || ipf.process(ImProcFunctions::Pipeline::NAVIGATOR, ImProcFunctions::Stage::STAGE_1, bufs_[0]);
                
            }
    
            // compute L channel histogram
            int x1, y1, x2, y2;
            params.crop.mapToResized(pW, pH, scale, x1, x2,  y1, y2);
        }
    
        readyphase++;
    
        if (todo & M_LUMACURVE) {
            bufs_[0]->copyTo(bufs_[1]);
            stop = stop || ipf.process(ImProcFunctions::Pipeline::NAVIGATOR, ImProcFunctions::Stage::STAGE_2, bufs_[1]);
        }

        if (todo & (M_LUMINANCE | M_COLOR)) {
            bufs_[1]->copyTo(bufs_[2]);
            stop = stop || ipf.process(ImProcFunctions::Pipeline::NAVIGATOR, ImProcFunctions::Stage::STAGE_3, bufs_[2]);
        }
    
        // Update the monitor color transform if necessary
        if ((todo & M_MONITOR) || (lastOutputProfile != params.icm.outputProfile) || lastOutputIntent != params.icm.outputIntent || lastOutputBPC != params.icm.outputBPC) {
            lastOutputProfile = params.icm.outputProfile;
            lastOutputIntent = params.icm.outputIntent;
            lastOutputBPC = params.icm.outputBPC;
            ipf.updateColorProfiles(monitorProfile, monitorIntent, softProof, gamutCheck);
        }
    }

    // process crop, if needed
    for (size_t i = 0; i < crops.size(); i++)
        if (crops[i]->hasListener() && (panningRelatedChange || (highDetailNeeded && options.prevdemo != PD_Sidecar) || (todo & (M_MONITOR | M_RGBCURVE | M_LUMACURVE)) || crops[i]->get_skip() == 1)) {
            crops[i]->update(todo);     // may call ourselves
        }

    if (panningRelatedChange || (todo & M_MONITOR)) {
        progress("Conversion to RGB...", 100 * readyphase / numofphases);

        if ((todo != CROP && todo != MINUPDATE) || (todo & M_MONITOR)) {
            MyMutex::MyLock prevImgLock(previmg->getMutex());

            try {
                // Computing the preview image, i.e. converting from WCS->Monitor color space (soft-proofing disabled) or WCS->Printer profile->Monitor color space (soft-proofing enabled)
                ipf.lab2monitorRgb(bufs_[2], previmg);

                // Computing the internal image for analysis, i.e. conversion from WCS->Output profile
                delete workimg;
                workimg = ipf.lab2rgb(bufs_[2], 0, 0, pW, pH, params.icm);
            } catch (char * str) {
                progress("Error converting file...", 0);
                return;
            }
        }

        if (!resultValid) {
            resultValid = true;

            if (imageListener) {
                imageListener->setImage(previmg, scale, params.crop);
            }
        }

        if (imageListener)
            // TODO: The WB tool should be advertised too in order to get the AutoWB's temp and green values
        {
            imageListener->imageReady(params.crop);
        }

        readyphase++;

        if (hListener) {
            updateLRGBHistograms();
            hListener->histogramChanged(histRed, histGreen, histBlue, histLuma, histToneCurve, histLCurve, histCCurve, /*histCLurve, histLLCurve,*/ histLCAM, histCCAM, histRedRaw, histGreenRaw, histBlueRaw, histChroma, histLRETI);
        }
    }
    if (orig_prev != oprevi) {
        delete oprevi;
        oprevi = nullptr;
    }

    
}


void ImProcCoordinator::freeAll()
{

    if (allocated) {
        if (orig_prev != oprevi) {
            delete oprevi;
        }

        oprevi    = nullptr;
        delete orig_prev;
        orig_prev = nullptr;
        for (int i = 3; i > 0; --i) {
            if (bufs_[i-1]) {
                delete bufs_[i-1];
                bufs_[i-1] = nullptr;
            }
        }

        if (imageListener) {
            imageListener->delImage(previmg);
        } else {
            delete previmg;
        }

        delete workimg;

    }

    allocated = false;
}

/** @brief Handles image buffer (re)allocation and trigger sizeChanged of SizeListener[s]
 * If the scale change, this method will free all buffers and reallocate ones of the new size.
 * It will then tell to the SizeListener that size has changed (sizeChanged)
 *
 * @param prevscale New Preview's scale.
 */
void ImProcCoordinator::setScale(int prevscale)
{

    tr = getCoarseBitMask(params.coarse);

    int nW, nH;
    imgsrc->getFullSize(fw, fh, tr);

    prevscale++;

    do {
        prevscale--;
        PreviewProps pp(0, 0, fw, fh, prevscale);
        imgsrc->getSize(pp, nW, nH);
    } while (nH < 400 && prevscale > 1 && (nW * nH < 1000000));  // sctually hardcoded values, perhaps a better choice is possible

    if (nW != pW || nH != pH) {

        freeAll();

        pW = nW;
        pH = nH;

        orig_prev = new Imagefloat(pW, pH);
        oprevi = orig_prev;
        for (int i = 0; i < 3; ++i) {
            bufs_[i] = new Imagefloat(pW, pH);
        }
        previmg = new Image8(pW, pH);
        workimg = new Image8(pW, pH);

        allocated = true;
    }

    scale = prevscale;
    resultValid = false;
    fullw = fw;
    fullh = fh;

    orig_prev->assignColorSpace(params.icm.workingProfile);
    if (oprevi && oprevi != orig_prev) {
        oprevi->assignColorSpace(params.icm.workingProfile);
    }
    for (int i = 0; i < 3; ++i) {
        bufs_[i]->assignColorSpace(params.icm.workingProfile);
    }
    
    if (!sizeListeners.empty())
        for (size_t i = 0; i < sizeListeners.size(); i++) {
            sizeListeners[i]->sizeChanged(fullw, fullh, fw, fh);
        }
}


void ImProcCoordinator::updateLRGBHistograms()
{

    int x1, y1, x2, y2;
    params.crop.mapToResized(pW, pH, scale, x1, x2, y1, y2);

#ifdef _OPENMP
    #pragma omp parallel sections
#endif
    {
#ifdef _OPENMP
        #pragma omp section
#endif
        {
            histLuma.clear();
            histChroma.clear();

            for (int i = y1; i < y2; i++)
                for (int j = x1; j < x2; j++)
                {
                    float L, a, b;
                    bufs_[2]->getLab(i, j, L, a, b);
                    histChroma[(int)(sqrtf(SQR(a) + SQR(b)) / 188.f)]++;      //188 = 48000/256
                    histLuma[(int)(L / 128.f)]++;
                }
        }
// #ifdef _OPENMP
//         #pragma omp section
// #endif
//         {
//             histLuma.clear();

//             for (int i = y1; i < y2; i++)
//                 for (int j = x1; j < x2; j++)
//                 {
//                     histLuma[(int)(lab->L[i][j] / 128.f)]++;
//                 }
//         }
#ifdef _OPENMP
        #pragma omp section
#endif
        {
            histRed.clear();
            histGreen.clear();
            histBlue.clear();

            for (int i = y1; i < y2; i++)
            {
                int ofs = (i * pW + x1) * 3;

                for (int j = x1; j < x2; j++) {
                    int r = workimg->data[ofs++];
                    int g = workimg->data[ofs++];
                    int b = workimg->data[ofs++];

                    histRed[r]++;
                    histGreen[g]++;
                    histBlue[b]++;
                }
            }
        }
    }

}

void ImProcCoordinator::progress(Glib::ustring str, int pr)
{

    /*  if (plistener) {
        plistener->setProgressStr (str);
        plistener->setProgress ((double)pr / 100.0);
      }*/
}

bool ImProcCoordinator::getAutoWB(double& temp, double& green, double equal)
{

    if (imgsrc) {
        if (lastAwbEqual != equal) {
// Issue 2500            MyMutex::MyLock lock(minit);  // Also used in crop window
            double rm, gm, bm;
            imgsrc->getAutoWBMultipliers(rm, gm, bm);

            if (rm != -1) {
                autoWB.update(rm, gm, bm, equal);
                lastAwbEqual = equal;
            } else {
                lastAwbEqual = -1.;
                autoWB.useDefaults(equal);
            }
        }

        temp = autoWB.getTemp();
        green = autoWB.getGreen();
        return true;
    } else {
        //temp = autoWB.getTemp();
        temp = -1.0;
        green = -1.0;
        return false;
    }
}

void ImProcCoordinator::getCamWB(double& temp, double& green)
{

    if (imgsrc) {
        temp = imgsrc->getWB().getTemp();
        green = imgsrc->getWB().getGreen();
    }
}

void ImProcCoordinator::getSpotWB(int x, int y, int rect, double& temp, double& tgreen)
{

    ColorTemp ret;

    {
        MyMutex::MyLock lock(mProcessing);
        std::vector<Coord2D> points, red, green, blue;

        for (int i = y - rect; i <= y + rect; i++)
            for (int j = x - rect; j <= x + rect; j++) {
                points.push_back(Coord2D(j, i));
            }

        ipf.transCoord(fw, fh, points, red, green, blue);

        int tr = getCoarseBitMask(params.coarse);

        ret = imgsrc->getSpotWB(red, green, blue, tr, params.wb.equal);
        currWB = ColorTemp(params.wb.temperature, params.wb.green, params.wb.equal, params.wb.method);
        //double rr,gg,bb;
        //currWB.getMultipliers(rr,gg,bb);

    } // end of mutex lockong

    if (ret.getTemp() > 0) {
        temp = ret.getTemp();
        tgreen = ret.getGreen();
    } else {
        temp = currWB.getTemp();
        tgreen = currWB.getGreen();
    }
}

void ImProcCoordinator::getAutoCrop(double ratio, int &x, int &y, int &w, int &h)
{

    MyMutex::MyLock lock(mProcessing);

    LensCorrection *pLCPMap = nullptr;

    if (params.lensProf.useLcp() && imgsrc->getMetaData()->getFocalLen() > 0) {
        const std::shared_ptr<LCPProfile> pLCPProf = LCPStore::getInstance()->getProfile(params.lensProf.lcpFile);

        if (pLCPProf) pLCPMap = new LCPMapper(pLCPProf, imgsrc->getMetaData()->getFocalLen(), imgsrc->getMetaData()->getFocalLen35mm(), imgsrc->getMetaData()->getFocusDist(),
                                                  0, false, params.lensProf.useDist, fullw, fullh, params.coarse, imgsrc->getRotateDegree());
    }

    double fillscale = ipf.getTransformAutoFill(fullw, fullh, pLCPMap);

    if (ratio > 0) {
        w = fullw * fillscale;
        h = w / ratio;

        if (h > fullh * fillscale) {
            h = fullh * fillscale;
            w = h * ratio;
        }
    } else {
        w = fullw * fillscale;
        h = fullh * fillscale;
    }

    x = (fullw - w) / 2;
    y = (fullh - h) / 2;

    if (params.perspective.enabled && !params.commonTrans.autofill) {
        int xx, yy, ww, hh;
        PerspectiveCorrection::autocrop(w, h, ratio > 0, params.perspective, imgsrc->getMetaData(), xx, yy, ww, hh);
        x += xx;
        y += yy;
        w = ww;
        h = hh;
    }
}


bool ImProcCoordinator::getFilmNegativeExponents(int xA, int yA, int xB, int yB, std::array<float, 3>& newExps)
{
    MyMutex::MyLock lock(mProcessing);

    const auto xlate =
        [this](int x, int y) -> Coord2D
        {
            const std::vector<Coord2D> points = {Coord2D(x, y)};

            std::vector<Coord2D> red;
            std::vector<Coord2D> green;
            std::vector<Coord2D> blue;
            ipf.transCoord(fw, fh, points, red, green, blue);

            return green[0];
        };

    const int tr = getCoarseBitMask(params.coarse);

    const Coord2D p1 = xlate(xA, yA);
    const Coord2D p2 = xlate(xB, yB);

    return imgsrc->getFilmNegativeExponents(p1, p2, tr, params.filmNegative, newExps);
}


void ImProcCoordinator::setMonitorProfile(const Glib::ustring& profile, RenderingIntent intent)
{
    monitorProfile = profile;
    monitorIntent = intent;
}

void ImProcCoordinator::getMonitorProfile(Glib::ustring& profile, RenderingIntent& intent) const
{
    profile = monitorProfile;
    intent = monitorIntent;
}

void ImProcCoordinator::setSoftProofing(bool softProof, bool gamutCheck)
{
    this->softProof = softProof;
    this->gamutCheck = gamutCheck;
}

void ImProcCoordinator::getSoftProofing(bool &softProof, bool &gamutCheck)
{
    softProof = this->softProof;
    gamutCheck = this->gamutCheck;
}

void ImProcCoordinator::setSharpMask (bool sharpMask)
{
    this->sharpMask = sharpMask;
}

void ImProcCoordinator::saveInputICCReference(const Glib::ustring& fname, bool apply_wb)
{

    MyMutex::MyLock lock(mProcessing);

    int fW, fH;

    int tr = getCoarseBitMask(params.coarse);

    imgsrc->getFullSize(fW, fH, tr);
    PreviewProps pp(0, 0, fW, fH, 1);
    ProcParams ppar = params;
    ppar.exposure.hrmode = procparams::ExposureParams::HR_OFF;
    ppar.icm.inputProfile = "(none)";
    Imagefloat* im = new Imagefloat(fW, fH);
    im->assignColorSpace(ppar.icm.workingProfile);
    imgsrc->preprocess(ppar.raw, ppar.lensProf, ppar.coarse);
    double dummy = 0.0;
    imgsrc->demosaic(ppar.raw, false, dummy);
    ColorTemp currWB = ColorTemp(params.wb.temperature, params.wb.green, params.wb.equal, params.wb.method);

    if (params.wb.method == "Camera") {
        currWB = imgsrc->getWB();
    } else if (params.wb.method == "Auto") {
        if (lastAwbEqual != params.wb.equal) {
            double rm, gm, bm;
            imgsrc->getAutoWBMultipliers(rm, gm, bm);

            if (rm != -1.) {
                autoWB.update(rm, gm, bm, params.wb.equal);
                lastAwbEqual = params.wb.equal;
            } else {
                lastAwbEqual = -1.;
                autoWB.useDefaults(params.wb.equal);
            }
        }

        currWB = autoWB;
    }

    if (!apply_wb) {
        currWB = ColorTemp(); // = no white balance
    }

    imgsrc->getImage(currWB, tr, im, pp, ppar.exposure, ppar.raw);
    ImProcFunctions ipf(&ppar, true);

    if (ipf.needsTransform()) {
        Imagefloat* trImg = new Imagefloat(fW, fH, im);
        ipf.transform(im, trImg, 0, 0, 0, 0, fW, fH, fW, fH,
                      imgsrc->getMetaData(), imgsrc->getRotateDegree(), true);
        delete im;
        im = trImg;
    }

    if (params.crop.enabled) {
        Imagefloat *tmpim = new Imagefloat(params.crop.w, params.crop.h, im);
        int cx = params.crop.x;
        int cy = params.crop.y;
        int cw = params.crop.w;
        int ch = params.crop.h;
#ifdef _OPENMP
        #pragma omp parallel for
#endif

        for (int i = cy; i < cy + ch; i++) {
            for (int j = cx; j < cx + cw; j++) {
                tmpim->r(i - cy, j - cx) = im->r(i, j);
                tmpim->g(i - cy, j - cx) = im->g(i, j);
                tmpim->b(i - cy, j - cx) = im->b(i, j);
            }
        }

        delete im;
        im = tmpim;
    }

    // image may contain out of range samples, clip them to avoid wrap-arounds
#ifdef _OPENMP
    #pragma omp parallel for
#endif

    for (int i = 0; i < im->getHeight(); i++) {
        for (int j = 0; j < im->getWidth(); j++) {
            im->r(i, j) = CLIP(im->r(i, j));
            im->g(i, j) = CLIP(im->g(i, j));
            im->b(i, j) = CLIP(im->b(i, j));
        }
    }

    int imw, imh;
    double tmpScale = ipf.resizeScale(&params, fW, fH, imw, imh);

    if (tmpScale != 1.0) {
        Imagefloat* tempImage = new Imagefloat(imw, imh, im);
        ipf.resize(im, tempImage, tmpScale);
        delete im;
        im = tempImage;
    }

    // im->setMetadata(imgsrc->getMetaData()->getRootExifData());
    im->setMetadata(Exiv2Metadata(imgsrc->getFileName(), false));

    im->saveTIFF(fname, 16, false, true);
    delete im;

    if (plistener) {
        plistener->setProgressState(false);
    }

    //im->saveJPEG (fname, 85);
}

void ImProcCoordinator::stopProcessing()
{

    updaterThreadStart.lock();

    if (updaterRunning && thread) {
        changeSinceLast = 0;
        thread->join();
    }

    updaterThreadStart.unlock();
}

void ImProcCoordinator::startProcessing()
{

#undef THREAD_PRIORITY_NORMAL

    if (!destroying) {
        if (!updaterRunning) {
            updaterThreadStart.lock();
            thread = nullptr;
            updaterRunning = true;
            updaterThreadStart.unlock();

            //batchThread->yield(); //the running batch should wait other threads to avoid conflict

            thread = Glib::Thread::create(sigc::mem_fun(*this, &ImProcCoordinator::process), 0, true, true, Glib::THREAD_PRIORITY_NORMAL);

        }
    }
}

void ImProcCoordinator::startProcessing(int changeCode)
{
    paramsUpdateMutex.lock();
    changeSinceLast |= changeCode;
    paramsUpdateMutex.unlock();

    startProcessing();
}

void ImProcCoordinator::process()
{
    if (plistener) {
        plistener->setProgressState(true);
        ipf.setProgressListener(plistener, crops.size() + 1);
    }

    paramsUpdateMutex.lock();

    while (changeSinceLast) {
        const bool panningRelatedChange = true;
        params = nextParams;
        int change = changeSinceLast;
        changeSinceLast = 0;
        paramsUpdateMutex.unlock();

        // M_VOID means no update, and is a bit higher that the rest
        if (change & (M_VOID - 1)) {
            updatePreviewImage(change, panningRelatedChange);
        }

        paramsUpdateMutex.lock();
    }

    paramsUpdateMutex.unlock();
    updaterRunning = false;

    if (plistener) {
        plistener->setProgressState(false);
    }
}

ProcParams* ImProcCoordinator::beginUpdateParams()
{
    paramsUpdateMutex.lock();

    return &nextParams;
}

void ImProcCoordinator::endUpdateParams(ProcEvent change)
{
    int action = RefreshMapper::getInstance()->getAction(change);
    endUpdateParams(action);
}

void ImProcCoordinator::endUpdateParams(int changeFlags)
{
    changeSinceLast |= changeFlags;

    paramsUpdateMutex.unlock();
    startProcessing();
}

bool ImProcCoordinator::getHighQualComputed()
{
    // this function may only be called from detail windows
    if (!highQualityComputed) {
        if (options.prevdemo == PD_Sidecar) {
            // we already have high quality preview
            setHighQualComputed();
        } else {
            for (size_t i = 0; i < crops.size() - 1; ++i) { // -1, because last entry is the freshly created detail window
                if (crops[i]->get_skip() == 1) {   // there is at least one crop with skip == 1 => we already have high quality preview
                    setHighQualComputed();
                    break;
                }
            }
        }
    }

    return highQualityComputed;
}

void ImProcCoordinator::setHighQualComputed()
{
    highQualityComputed = true;
}


bool ImProcCoordinator::getDeltaELCH(EditUniqueID id, int x, int y, float &L, float &C, float &H)
{
    int change = ipf.setDeltaEData(id, x, y);
    if (!change) {
        return false;
    }
    startProcessing(change);

    bool ret = false;
    updaterThreadStart.lock();
    if (updaterRunning && thread) {
        thread->join();
        if (ipf.deltaE.ok) {
            ret = true;
            L = ipf.deltaE.L;
            C = ipf.deltaE.C;
            H = ipf.deltaE.H;
        }
    }
    ipf.setDeltaEData(EUID_None, -1, -1);
    updaterThreadStart.unlock();

    return ret;
}

} // namespace rtengine
