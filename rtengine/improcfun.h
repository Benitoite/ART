/* -*- C++ -*-
 *  
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
#ifndef _IMPROCFUN_H_
#define _IMPROCFUN_H_

#include "imagefloat.h"
#include "image16.h"
#include "image8.h"
#include "procparams.h"
#include "coord2d.h"
#include "color.h"
#include "labimage.h"
#include "LUT.h"
#include "lcp.h"
#include "dcp.h"
#include "curves.h"
#include "cplx_wavelet_dec.h"
#include "pipettebuffer.h"
#include "gamutwarning.h"

namespace rtengine {

using namespace procparams;

struct ImProcData {
    const ProcParams *params;
    double scale;
    bool multiThread;

    explicit ImProcData(const ProcParams *p=nullptr, double s=1.0, bool m=true):
        params(p), scale(s), multiThread(m) {}
};


class ImProcFunctions {
public:
    //----------------------------------------------------------------------
    // constructor/destructor and initialization/state manipulation
    //----------------------------------------------------------------------
    ImProcFunctions(const ProcParams* iparams, bool imultiThread=true);
    ~ImProcFunctions();
    
    void setScale(double iscale);

    void updateColorProfiles(const Glib::ustring& monitorProfile, RenderingIntent monitorIntent, bool softProof, bool gamutCheck);

    void setDCPProfile(DCPProfile *dcp, const DCPProfile::ApplyState &as)
    {
        dcpProf = dcp;
        dcpApplyState = &as;
    }

    void setPipetteBuffer(PipetteBuffer *pb)
    {
        pipetteBuffer = pb;
    }
    //----------------------------------------------------------------------

    //----------------------------------------------------------------------
    // pipeline management
    //----------------------------------------------------------------------
    enum class Stage {
        STAGE_0,
        STAGE_1,
        STAGE_2,
        STAGE_3
    };
    enum class Pipeline {
        THUMBNAIL,
        NAVIGATOR,
        PREVIEW,
        OUTPUT
    };
    bool process(Pipeline pipeline, Stage stage, Imagefloat *img);

    void setViewport(int ox, int oy, int fw, int fh);
    void setOutputHistograms(LUTu *histToneCurve, LUTu *histCCurve, LUTu *histLCurve);
    void setShowSharpeningMask(bool yes);
    //----------------------------------------------------------------------
    
    //----------------------------------------------------------------------
    // image processing operations
    //----------------------------------------------------------------------
    void firstAnalysis(const Imagefloat* const working, const ProcParams &params, LUTu & vhist16);

    void labAdjustments(Imagefloat *rgb);
    bool sharpening(Imagefloat *rgb, const SharpeningParams &sharpenParam, bool showMask=false);
    void transform(Imagefloat* original, Imagefloat* transformed, int cx, int cy, int sx, int sy, int oW, int oH, int fW, int fH, const FramesMetaData *metadata, int rawRotationDeg, bool highQuality);    
    void resize(Imagefloat* src, Imagefloat* dst, float dScale);
    void Lanczos(const LabImage* src, LabImage* dst, float scale);
    void Lanczos(Imagefloat *src, Imagefloat *dst, float scale);
    void impulsedenoise(Imagefloat *rgb);   //Emil's impulse denoise
    bool textureBoost(Imagefloat *rgb);

    struct DenoiseInfoStore {
        DenoiseInfoStore () : chM (0), max_r{}, max_b{}, ch_M{}, valid (false)  {}
        float chM;
        float max_r[9];
        float max_b[9];
        float ch_M[9];
        bool valid;
    };
    void denoiseComputeParams(ImageSource *imgsrc, const ColorTemp &currWB, DenoiseInfoStore &store, procparams::DenoiseParams &dnparams);
    void denoise(ImageSource *imgsrc, const ColorTemp &currWB, Imagefloat *img, const DenoiseInfoStore &store, const procparams::DenoiseParams &dnparams);
    
    void defringe(Imagefloat *rgb);
    void dehaze(Imagefloat *rgb);
    void dynamicRangeCompression(Imagefloat *rgb);
    bool localContrast(Imagefloat *rgb);
    void toneEqualizer(Imagefloat *rgb);
    void softLight(Imagefloat *rgb);
    bool colorCorrection(Imagefloat *rgb);
    void logEncoding(Imagefloat *rgb);
    // bool contrastByDetailLevels(Imagefloat *rgb);
    void filmGrain(Imagefloat *rgb);
    bool guidedSmoothing(Imagefloat *rgb);
    void hslEqualizer(Imagefloat *rgb);
    void channelMixer(Imagefloat *rgb);
    void exposure(Imagefloat *rgb);
    void rgbCurves(Imagefloat *rgb);
    void blackAndWhite(Imagefloat *rgb);
    void toneCurve(Imagefloat *img);
    void saturationVibrance(Imagefloat *img);
    void filmSimulation(Imagefloat *img);
    void creativeGradients(Imagefloat *img);
    //----------------------------------------------------------------------

    //----------------------------------------------------------------------
    // Lab/RGB conversion
    //----------------------------------------------------------------------
    void lab2monitorRgb(Imagefloat *img, Image8* image);
    
    Image8 *lab2rgb(Imagefloat *img, int cx, int cy, int cw, int ch, const procparams::ColorManagementParams &icm, bool consider_histogram_settings = true);

    Imagefloat *lab2rgbOut(Imagefloat *img, int cx, int cy, int cw, int ch, const procparams::ColorManagementParams &icm);

    void rgb2lab(Imagefloat &src, LabImage &dst, const Glib::ustring &workingSpace);
    void rgb2lab(Imagefloat &src, LabImage &dst) { rgb2lab(src, dst, params->icm.workingProfile); }
    
    void lab2rgb(const LabImage &src, Imagefloat &dst, const Glib::ustring &workingSpace);    
    void lab2rgb(const LabImage &src, Imagefloat &dst) { lab2rgb(src, dst, params->icm.workingProfile); }
    //----------------------------------------------------------------------

    //----------------------------------------------------------------------
    // misc helper functions for image processing ops
    //----------------------------------------------------------------------
    bool needsLuminanceOnly();
    bool needsTransform();
    bool needsPCVignetting();
    
    float resizeScale(const ProcParams* params, int fw, int fh, int &imw, int &imh);
    
    void getAutoLog(ImageSource *imgsrc, procparams::LogEncodingParams &params);
    
    static void getAutoExp(const LUTu & histogram, int histcompr, double clip, double& expcomp, int& bright, int& contr, int& black, int& hlcompr, int& hlcomprthresh);
    
    static double getAutoDistor(const Glib::ustring& fname, int thumb_size);
    bool transCoord(int W, int H, int x, int y, int w, int h, int& xv, int& yv, int& wv, int& hv, double ascaleDef = -1, const LensCorrection *pLCPMap = nullptr);
    bool transCoord(int W, int H, const std::vector<Coord2D> &src, std::vector<Coord2D> &red,  std::vector<Coord2D> &green, std::vector<Coord2D> &blue, double ascaleDef = -1, const LensCorrection *pLCPMap = nullptr);
    double getTransformAutoFill(int oW, int oH, const LensCorrection *pLCPMap = nullptr);
    //----------------------------------------------------------------------

    class DeltaEData {
    public:
        bool ok;
        float L;
        float C;
        float H;
        double x;
        double y;

        DeltaEData():
            ok(false), L(0), C(0), H(0), x(-1), y(-1) {}
    };
    DeltaEData deltaE;
    int setDeltaEData(EditUniqueID id, double x, double y);

private:
    cmsHTRANSFORM monitorTransform;
    std::unique_ptr<GamutWarning> gamutWarning;

    const ProcParams* params;
    double scale;
    bool multiThread;

    DCPProfile *dcpProf;
    const DCPProfile::ApplyState *dcpApplyState;

    PipetteBuffer *pipetteBuffer;

    double lumimul[3];

    int offset_x;
    int offset_y;
    int full_width;
    int full_height;

    LUTu *histToneCurve;
    LUTu *histCCurve;
    LUTu *histLCurve;

    bool show_sharpening_mask;
    
private:
    void transformLuminanceOnly(Imagefloat* original, Imagefloat* transformed, int cx, int cy, int oW, int oH, int fW, int fH, bool creative);
    void transformGeneral(bool highQuality, Imagefloat *original, Imagefloat *transformed, int cx, int cy, int sx, int sy, int oW, int oH, int fW, int fH, const LensCorrection *pLCPMap);
    void transformLCPCAOnly(Imagefloat *original, Imagefloat *transformed, int cx, int cy, const LensCorrection *pLCPMap);

    bool needsCA();
    bool needsDistortion();
    bool needsRotation();
    bool needsPerspective();
    bool needsGradient();
    bool needsVignetting();
    bool needsLCP();
    bool needsLensfun();
};


} // namespace rtengine

#endif
