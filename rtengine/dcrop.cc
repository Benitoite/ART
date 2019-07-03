/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
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
#include "dcrop.h"
#include "curves.h"
#include "mytime.h"
#include "refreshmap.h"
#include "rt_math.h"

namespace
{

// "ceil" rounding
template<typename T>
constexpr T skips(T a, T b)
{
    return a / b + static_cast<bool>(a % b);
}

}

namespace rtengine
{

extern const Settings* settings;

Crop::Crop(ImProcCoordinator* parent, EditDataProvider *editDataProvider, bool isDetailWindow)
    : PipetteBuffer(editDataProvider), origCrop(nullptr), laboCrop(nullptr), labnCrop(nullptr),
      cropImg (nullptr), transCrop (nullptr), 
      updating(false), newUpdatePending(false), skip(10),
      cropx(0), cropy(0), cropw(-1), croph(-1),
      trafx(0), trafy(0), trafw(-1), trafh(-1),
      rqcropx(0), rqcropy(0), rqcropw(-1), rqcroph(-1),
      borderRequested(32), upperBorder(0), leftBorder(0),
      cropAllocated(false),
      cropImageListener(nullptr), parent(parent), isDetailWindow(isDetailWindow)
{
    parent->crops.push_back(this);
}

Crop::~Crop()
{

    MyMutex::MyLock cropLock(cropMutex);

    std::vector<Crop*>::iterator i = std::find(parent->crops.begin(), parent->crops.end(), this);

    if (i != parent->crops.end()) {
        parent->crops.erase(i);
    }

    MyMutex::MyLock processingLock(parent->mProcessing);
    freeAll();
}

void Crop::destroy()
{
    MyMutex::MyLock lock(cropMutex);
    MyMutex::MyLock processingLock(parent->mProcessing);
    freeAll();
}

void Crop::setListener(DetailedCropListener* il)
{
    // We can make reads in the IF, because the mProcessing lock is only needed for change
    if (cropImageListener != il) {
        MyMutex::MyLock lock(cropMutex);
        cropImageListener = il;
    }
}

EditUniqueID Crop::getCurrEditID()
{
    EditSubscriber *subscriber = PipetteBuffer::dataProvider ? PipetteBuffer::dataProvider->getCurrSubscriber() : nullptr;
    return subscriber ? subscriber->getEditID() : EUID_None;
}

/*
 * Delete the edit image buffer if there's no subscriber anymore.
 * If allocation has to be done, it is deferred to Crop::update
 */
void Crop::setEditSubscriber(EditSubscriber* newSubscriber)
{
    MyMutex::MyLock lock(cropMutex);

    // At this point, editCrop.dataProvider->currSubscriber is the old subscriber
    EditSubscriber *oldSubscriber = PipetteBuffer::dataProvider ? PipetteBuffer::dataProvider->getCurrSubscriber() : nullptr;

    if (newSubscriber == nullptr || (oldSubscriber != nullptr && oldSubscriber->getPipetteBufferType() != newSubscriber->getPipetteBufferType())) {
        if (PipetteBuffer::imgFloatBuffer != nullptr) {
            delete PipetteBuffer::imgFloatBuffer;
            PipetteBuffer::imgFloatBuffer = nullptr;
        }

        if (PipetteBuffer::LabBuffer != nullptr) {
            delete PipetteBuffer::LabBuffer;
            PipetteBuffer::LabBuffer = nullptr;
        }

        if (PipetteBuffer::singlePlaneBuffer.getWidth() != -1) {
            PipetteBuffer::singlePlaneBuffer.flushData();
        }
    }

    // If oldSubscriber == NULL && newSubscriber != NULL && newSubscriber->getEditingType() == ET_PIPETTE-> the image will be allocated when necessary
}

bool Crop::hasListener()
{
    MyMutex::MyLock cropLock(cropMutex);
    return cropImageListener;
}

void Crop::update(int todo)
{
    MyMutex::MyLock cropLock(cropMutex);

    ProcParams& params = parent->params;
//       CropGUIListener* cropgl;

    // No need to update todo here, since it has already been changed in ImprocCoordinator::updatePreviewImage,
    // and Crop::update ask to do ALL anyway

    // give possibility to the listener to modify crop window (as the full image dimensions are already known at this point)
    int wx, wy, ww, wh, ws;
    const bool overrideWindow = cropImageListener;

    if (overrideWindow) {
        cropImageListener->getWindow(wx, wy, ww, wh, ws);
    }

    // re-allocate sub-images and arrays if their dimensions changed
    bool needsinitupdate = false;

    if (!overrideWindow) {
        needsinitupdate = setCropSizes(rqcropx, rqcropy, rqcropw, rqcroph, skip, true);
    } else {
        needsinitupdate = setCropSizes(wx, wy, ww, wh, ws, true);     // this set skip=ws
    }

    // it something has been reallocated, all processing steps have to be performed
    if (needsinitupdate || (todo & M_HIGHQUAL)) {
        todo = ALL;
    }

    // Tells to the ImProcFunctions' tool what is the preview scale, which may lead to some simplifications
    parent->ipf.setScale(skip);
    parent->ipf.setPipetteBuffer(this);

    Imagefloat* baseCrop = origCrop;

    bool needstransform  = parent->ipf.needsTransform();
    bool show_denoise = params.denoise.enabled && (skip == 1 || options.denoiseZoomedOut);

    if (todo & (M_INIT | M_LINDENOISE | M_HDR)) {
        MyMutex::MyLock lock(parent->minit);  // Also used in improccoord

        int tr = getCoarseBitMask(params.coarse);

        if (!needsinitupdate) {
            setCropSizes(rqcropx, rqcropy, rqcropw, rqcroph, skip, true);
        }

        PreviewProps pp(trafx, trafy, trafw * skip, trafh * skip, skip);
        parent->imgsrc->getImage(parent->currWB, tr, origCrop, pp, params.toneCurve, params.raw);

        if (show_denoise) {
            parent->ipf.denoiseComputeParams(parent->imgsrc, parent->currWB, parent->denoiseInfoStore, params.denoise);
        }

        if ((!isDetailWindow) && parent->adnListener && show_denoise) {
            parent->adnListener->chromaChanged(params.denoise.chrominance, params.denoise.chrominanceRedGreen, params.denoise.chrominanceBlueYellow);
        }

        if ((todo & M_LINDENOISE) && show_denoise) {
            parent->ipf.denoise(parent->imgsrc, parent->currWB, origCrop, parent->denoiseInfoStore, params.denoise);

            if (parent->adnListener && params.denoise.chrominanceMethod == DenoiseParams::ChrominanceMethod::AUTOMATIC) {
                parent->adnListener->chromaChanged(params.denoise.chrominance, params.denoise.chrominanceRedGreen, params.denoise.chrominanceBlueYellow);
            }                
        }

        parent->imgsrc->convertColorSpace(origCrop, params.icm, parent->currWB);
    }

    // has to be called after setCropSizes! Tools prior to this point can't handle the Edit mechanism, but that shouldn't be a problem.
    createBuffer(cropw, croph);

    std::unique_ptr<Imagefloat> drCompCrop;

    if ((todo & M_HDR) && (params.fattal.enabled || params.dehaze.enabled)) {
        Imagefloat *f = origCrop;
        int fw = skips(parent->fw, skip);
        int fh = skips(parent->fh, skip);
        bool need_cropping = false;
        bool need_drcomp = true;

        if (trafx || trafy || trafw != fw || trafh != fh) {
            need_cropping = true;

            // fattal needs to work on the full image. So here we get the full
            // image from imgsrc, and replace the denoised crop in case
            if (!params.denoise.enabled && skip == 1 && parent->drcomp_11_dcrop_cache) {
                f = parent->drcomp_11_dcrop_cache;
                need_drcomp = false;
            } else {
                f = new Imagefloat(fw, fh);
                drCompCrop.reset(f);
                PreviewProps pp(0, 0, parent->fw, parent->fh, skip);
                int tr = getCoarseBitMask(params.coarse);
                parent->imgsrc->getImage(parent->currWB, tr, f, pp, params.toneCurve, params.raw);
                parent->imgsrc->convertColorSpace(f, params.icm, parent->currWB);

                if (params.denoise.enabled) {
                    // copy the denoised crop
                    int oy = trafy / skip;
                    int ox = trafx / skip;
#ifdef _OPENMP
                    #pragma omp parallel for
#endif

                    for (int y = 0; y < baseCrop->getHeight(); ++y) {
                        int dy = oy + y;

                        for (int x = 0; x < baseCrop->getWidth(); ++x) {
                            int dx = ox + x;
                            f->r(dy, dx) = baseCrop->r(y, x);
                            f->g(dy, dx) = baseCrop->g(y, x);
                            f->b(dy, dx) = baseCrop->b(y, x);
                        }
                    }
                } else if (skip == 1) {
                    parent->drcomp_11_dcrop_cache = f; // cache this globally
                    drCompCrop.release();
                }
            }
        }

        if (need_drcomp) {
            parent->ipf.dehaze(f);
            parent->ipf.dynamicRangeCompression(f);
        }

        // crop back to the size expected by the rest of the pipeline
        if (need_cropping) {
            Imagefloat *c = origCrop;

            int oy = trafy / skip;
            int ox = trafx / skip;
#ifdef _OPENMP
            #pragma omp parallel for
#endif

            for (int y = 0; y < trafh; ++y) {
                int cy = y + oy;

                for (int x = 0; x < trafw; ++x) {
                    int cx = x + ox;
                    c->r(y, x) = f->r(cy, cx);
                    c->g(y, x) = f->g(cy, cx);
                    c->b(y, x) = f->b(cy, cx);
                }
            }

            baseCrop = c;
        } else {
            baseCrop = f;
        }
    }

    // transform
    if (needstransform) {
        if (!transCrop) {
            transCrop = new Imagefloat(cropw, croph);
        }

        if (needstransform)
            parent->ipf.transform(baseCrop, transCrop, cropx / skip, cropy / skip, trafx / skip, trafy / skip, skips(parent->fw, skip), skips(parent->fh, skip), parent->getFullWidth(), parent->getFullHeight(),
                                  parent->imgsrc->getMetaData(),
                                  parent->imgsrc->getRotateDegree(), false);
        else {
            baseCrop->copyData(transCrop);
        }

        if (transCrop) {
            baseCrop = transCrop;
        }
    } else {
        if (transCrop) {
            delete transCrop;
        }

        transCrop = nullptr;
    }

    if (todo & M_RGBCURVE) {
        Imagefloat *workingCrop = baseCrop;
        parent->ipf.rgbProc(workingCrop, laboCrop);
        
        if (workingCrop != baseCrop) {
            delete workingCrop;
        }
    }

    int offset_x = cropx / skip;
    int offset_y = cropy / skip;
    int full_width = parent->getFullWidth() / skip;
    int full_height = parent->getFullHeight() / skip;
    
    // apply luminance operations
    if (todo & M_LUMACURVE) {
        parent->ipf.labColorCorrectionRegions(laboCrop, offset_x, offset_y, full_width, full_height);
        parent->ipf.guidedSmoothing(laboCrop, offset_x, offset_y, full_width, full_height);
        parent->ipf.logEncoding(laboCrop);
        
        // bool utili = parent->utili;
        // bool autili = parent->autili;
        // bool butili = parent->butili;
        // bool ccutili = parent->ccutili;
        // bool clcutili = parent->clcutili;
        // bool cclutili = parent->cclutili;

        // LUTu dummy;
        // parent->ipf.chromiLuminanceCurve(1, laboCrop, laboCrop, parent->chroma_acurve, parent->chroma_bcurve, parent->satcurve, parent->lhskcurve,  parent->clcurve, parent->lumacurve, utili, autili, butili, ccutili, cclutili, clcutili, dummy, dummy);
        parent->ipf.labAdjustments(laboCrop);
    }
    
    if (todo & (M_LUMINANCE | M_COLOR)) {
        labnCrop->CopyFrom(laboCrop);

        parent->ipf.toneMapping(labnCrop, offset_x, offset_y, full_width, full_height);

        //parent->ipf.EPDToneMap(labnCrop, 5, 1);    //Go with much fewer than normal iterates for fast redisplay.
        // for all treatments Defringe, Sharpening, Contrast detail , Microcontrast they are activated if "CIECAM" function are disabled
        if (skip == 1) {
            parent->ipf.impulsedenoise(labnCrop);
            parent->ipf.defringe(labnCrop);
            parent->ipf.MLmicrocontrast (labnCrop);
            parent->ipf.sharpening (labnCrop, params.sharpening, parent->sharpMask);
        }

        parent->ipf.contrastByDetailLevels(labnCrop, offset_x, offset_y, full_width, full_height); 

        parent->ipf.softLight(labnCrop);
        parent->ipf.localContrast(labnCrop);
        parent->ipf.filmGrain(labnCrop, cropx / skip, cropy / skip, parent->getFullWidth() / skip, parent->getFullHeight() / skip);
        
    }

    // all pipette buffer processing should be finished now
    PipetteBuffer::setReady();

    // Computing the preview image, i.e. converting from lab->Monitor color space (soft-proofing disabled) or lab->Output profile->Monitor color space (soft-proofing enabled)
    parent->ipf.lab2monitorRgb(labnCrop, cropImg);

    if (cropImageListener) {
        // Computing the internal image for analysis, i.e. conversion from lab->Output profile (rtSettings.HistogramWorking disabled) or lab->WCS (rtSettings.HistogramWorking enabled)

        // internal image in output color space for analysis
        Image8 *cropImgtrue = parent->ipf.lab2rgb(labnCrop, 0, 0, cropw, croph, params.icm);

        int finalW = rqcropw;

        if (cropImg->getWidth() - leftBorder < finalW) {
            finalW = cropImg->getWidth() - leftBorder;
        }

        int finalH = rqcroph;

        if (cropImg->getHeight() - upperBorder < finalH) {
            finalH = cropImg->getHeight() - upperBorder;
        }

        Image8* final = new Image8(finalW, finalH);
        Image8* finaltrue = new Image8(finalW, finalH);

        for (int i = 0; i < finalH; i++) {
            memcpy(final->data + 3 * i * finalW, cropImg->data + 3 * (i + upperBorder)*cropw + 3 * leftBorder, 3 * finalW);
            memcpy(finaltrue->data + 3 * i * finalW, cropImgtrue->data + 3 * (i + upperBorder)*cropw + 3 * leftBorder, 3 * finalW);
        }

        cropImageListener->setDetailedCrop(final, finaltrue, params.icm, params.crop, rqcropx, rqcropy, rqcropw, rqcroph, skip);
        delete final;
        delete finaltrue;
        delete cropImgtrue;
    }
}

void Crop::freeAll()
{

    if (cropAllocated) {
        if (origCrop) {
            delete    origCrop;
            origCrop = nullptr;
        }

        if (transCrop) {
            delete    transCrop;
            transCrop = nullptr;
        }

        if (laboCrop) {
            delete    laboCrop;
            laboCrop = nullptr;
        }

        if (labnCrop) {
            delete    labnCrop;
            labnCrop = nullptr;
        }

        if (cropImg) {
            delete    cropImg;
            cropImg = nullptr;
        }

        PipetteBuffer::flush();
    }

    cropAllocated = false;
}


namespace {

bool check_need_larger_crop_for_transform(int fw, int fh, int x, int y, int w, int h, const ProcParams &params, double &adjust)
{
    if (x == 0 && y == 0 && w == fw && h == fh) {
        return false;
    }

    if (params.perspective.enabled) {
        adjust = 1; // TODO -- ask PerspectiveCorrection the right value
        return true;
    } else if (params.lensProf.useDist && (params.lensProf.useLensfun() || params.lensProf.useLcp())) {
        adjust = 0.15;
        return true;
    }

    return false;
}

} // namespace

/** @brief Handles crop's image buffer reallocation and trigger sizeChanged of SizeListener[s]
 * If the scale changes, this method will free all buffers and reallocate ones of the new size.
 * It will then tell to the SizeListener that size has changed (sizeChanged)
 */
bool Crop::setCropSizes(int rcx, int rcy, int rcw, int rch, int skip, bool internal)
{

    if (!internal) {
        cropMutex.lock();
    }

    bool changed = false;

    rqcropx = rcx;
    rqcropy = rcy;
    rqcropw = rcw;
    rqcroph = rch;

    // store and set requested crop size
    int rqx1 = LIM(rqcropx, 0, parent->fullw - 1);
    int rqy1 = LIM(rqcropy, 0, parent->fullh - 1);
    int rqx2 = rqx1 + rqcropw - 1;
    int rqy2 = rqy1 + rqcroph - 1;
    rqx2 = LIM(rqx2, 0, parent->fullw - 1);
    rqy2 = LIM(rqy2, 0, parent->fullh - 1);

    this->skip = skip;

    // add border, if possible
    int bx1 = rqx1 - skip * borderRequested;
    int by1 = rqy1 - skip * borderRequested;
    int bx2 = rqx2 + skip * borderRequested;
    int by2 = rqy2 + skip * borderRequested;
    // clip it to fit into image area
    bx1 = LIM(bx1, 0, parent->fullw - 1);
    by1 = LIM(by1, 0, parent->fullh - 1);
    bx2 = LIM(bx2, 0, parent->fullw - 1);
    by2 = LIM(by2, 0, parent->fullh - 1);
    int bw = bx2 - bx1 + 1;
    int bh = by2 - by1 + 1;

    // determine which part of the source image is required to compute the crop rectangle
    int orx, ory, orw, orh;
    orx = bx1;
    ory = by1;
    orw = bw;
    orh = bh;

    parent->ipf.transCoord(parent->fw, parent->fh, bx1, by1, bw, bh, orx, ory, orw, orh);

    double adjust = 0.f;
    if (check_need_larger_crop_for_transform(parent->fw, parent->fh, orx, ory, orw, orh, parent->params, adjust)) {
        // TODO - "adjust" is an estimate of the max distortion relative to the image size. It is hardcoded to be 15% for lens distortion correction, and 100% for perspective (because we don't know better yet -- ideally this should be calculated)
        int dW = int (double (parent->fw) * adjust / (2 * skip));
        int dH = int (double (parent->fh) * adjust / (2 * skip));
        int x1 = orx - dW;
        int x2 = orx + orw + dW;
        int y1 = ory - dH;
        int y2 = ory + orh + dH;

        if (x1 < 0) {
            x2 += -x1;
            x1 = 0;
        }

        if (x2 > parent->fw) {
            x1 -= x2 - parent->fw;
            x2 = parent->fw;
        }

        if (y1 < 0) {
            y2 += -y1;
            y1 = 0;
        }

        if (y2 > parent->fh) {
            y1 -= y2 - parent->fh;
            y2 = parent->fh;
        }

        orx = max(x1, 0);
        ory = max(y1, 0);
        orw = min(x2 - x1, parent->fw - orx);
        orh = min(y2 - y1, parent->fh - ory);
    }

    leftBorder  = skips(rqx1 - bx1, skip);
    upperBorder = skips(rqy1 - by1, skip);

    PreviewProps cp(orx, ory, orw, orh, skip);
    int orW, orH;
    parent->imgsrc->getSize(cp, orW, orH);

    trafx = orx;
    trafy = ory;

    int cw = skips(bw, skip);
    int ch = skips(bh, skip);

    EditType editType = ET_PIPETTE;

    if (const auto editProvider = PipetteBuffer::getDataProvider()) {
        if (const auto editSubscriber = editProvider->getCurrSubscriber()) {
            editType = editSubscriber->getEditingType();
        }
    }

    if (cw != cropw || ch != croph || orW != trafw || orH != trafh) {

        cropw = cw;
        croph = ch;
        trafw = orW;
        trafh = orH;

        if (!origCrop) {
            origCrop = new Imagefloat;
        }

        origCrop->allocate(trafw, trafh);  // Resizing the buffer (optimization)

        // if transCrop doesn't exist yet, it'll be created where necessary
        if (transCrop) {
            transCrop->allocate(cropw, croph);
        }

        if (laboCrop) {
            delete laboCrop;    // laboCrop can't be resized
        }

        laboCrop = new LabImage(cropw, croph);

        if (labnCrop) {
            delete labnCrop;    // labnCrop can't be resized
        }

        labnCrop = new LabImage(cropw, croph);

        if (!cropImg) {
            cropImg = new Image8;
        }

        cropImg->allocate(cropw, croph);  // Resizing the buffer (optimization)

        if (editType == ET_PIPETTE) {
            PipetteBuffer::resize(cropw, croph);
        } else if (PipetteBuffer::bufferCreated()) {
            PipetteBuffer::flush();
        }

        cropAllocated = true;

        changed = true;
    }

    cropx = bx1;
    cropy = by1;

    if (!internal) {
        cropMutex.unlock();
    }

    return changed;
}

/** @brief Look out if a new thread has to be started to process the update
  *
  * @return If true, a new updating thread has to be created. If false, the current updating thread will be used
  */
bool Crop::tryUpdate()
{
    bool needsNewThread = true;

    if (updating) {
        // tells to the updater thread that a new update is pending
        newUpdatePending = true;
        // no need for a new thread, the current one will do the job
        needsNewThread = false;
    } else
        // the crop is now being updated ...well, when fullUpdate will be called
    {
        updating = true;
    }

    return needsNewThread;
}

/* @brief Handles Crop updating in its own thread
 *
 * This method will cycle updates as long as Crop::newUpdatePending will be true. During the processing,
 * intermediary update will be automatically flushed by Crop::tryUpdate.
 *
 * This method is called when the visible part of the crop has changed (resize, zoom, etc..), so it needs a full update
 */
void Crop::fullUpdate()
{

    parent->updaterThreadStart.lock();

    if (parent->updaterRunning && parent->thread) {
        // Do NOT reset changes here, since in a long chain of events it will lead to chroma_scale not being updated,
        // causing Color::lab2rgb to return a black image on some opens
        //parent->changeSinceLast = 0;
        parent->thread->join();
    }

    if (parent->plistener) {
        parent->plistener->setProgressState(true);
    }

    // If there are more update request, the following WHILE will collect it
    newUpdatePending = true;

    while (newUpdatePending) {
        newUpdatePending = false;
        update(ALL);
    }

    updating = false;  // end of crop update

    if (parent->plistener) {
        parent->plistener->setProgressState(false);
    }

    parent->updaterThreadStart.unlock();
}

int Crop::get_skip()
{
    MyMutex::MyLock lock(cropMutex);
    return skip;
}

int Crop::getLeftBorder()
{
    MyMutex::MyLock lock(cropMutex);
    return leftBorder;
}

int Crop::getUpperBorder()
{
    MyMutex::MyLock lock(cropMutex);
    return upperBorder;
}

}
