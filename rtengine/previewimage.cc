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

#include "previewimage.h"
#include "iimage.h"
#include "utils.h"
#include "iimage.h"
#include "rtthumbnail.h"
#include "rawimagesource.h"
#include "stdimagesource.h"
#include "iccstore.h"


namespace rtengine {
extern const Settings *settings;
} // namespace rtengine


using namespace rtengine;
using namespace procparams;

PreviewImage::PreviewImage(const Glib::ustring &fname, const Glib::ustring &ext, int width, int height, bool enable_cms)
{
    auto lext = ext.lowercase();

    if (lext == "jpg" || lext == "jpeg" || lext == "png" || lext == "tif" || lext == "tiff") {
        img_.reset(load_img(fname, width, height));
    } else if (settings->thumbnail_inspector_mode == Settings::ThumbnailInspectorMode::RAW) {
        img_.reset(load_raw(fname, width, height));
        if (settings->thumbnail_inspector_raw_curve == Settings::ThumbnailInspectorRawCurve::RAW_CLIPPING) {
            enable_cms = false;
        }
    } else {
        img_.reset(load_raw_preview(fname, width, height));
    }

    if (img_) {
        try {
            previewImage = Cairo::ImageSurface::create(Cairo::FORMAT_RGB24, img_->getWidth(), img_->getHeight());
            previewImage->flush();
            render(enable_cms);
        } catch (std::exception &exc) {
            if (settings->verbose) {
                std::cout << "ERROR in creating PreviewImage: " << exc.what() << std::endl;
            }
            previewImage.clear();
        }
    }    
}


void PreviewImage::render(bool enable_cms)
{
    if (img_) {
        cmsHTRANSFORM xform = nullptr;
        if (enable_cms) {
            cmsHPROFILE mprof = ICCStore::getInstance()->getProfile(ICCStore::getInstance()->getDefaultMonitorProfileName());
            cmsHPROFILE iprof = ICCStore::getInstance()->getsRGBProfile();
            if (mprof) {
                xform = cmsCreateTransform(iprof, TYPE_RGB_8, mprof, TYPE_RGB_8, settings->monitorIntent, cmsFLAGS_NOCACHE | (settings->monitorBPC ? cmsFLAGS_BLACKPOINTCOMPENSATION : 0));
            }
        }
        const unsigned char *data = img_->data;
        int w = img_->getWidth();
        int h = img_->getHeight();

#ifdef _OPENMP
#       pragma omp parallel
#endif
        {
            std::vector<unsigned char> line(w * 3);
            unsigned char *buf = &line[0];
            
#ifdef _OPENMP
#           pragma omp for
#endif
            for (unsigned int i = 0; i < (unsigned int)(h); ++i) {
                const unsigned char *src = data + i * w * 3;
                unsigned char *dst = previewImage->get_data() + i * w * 4;

                if (xform) {
                    cmsDoTransform(xform, src, buf, w);
                    src = buf;
                }
                
                for (unsigned int j = 0; j < (unsigned int)(w); ++j) {
                    unsigned char r = *(src++);
                    unsigned char g = *(src++);
                    unsigned char b = *(src++);

                    poke255_uc(dst, r, g, b);
                }
            }
        }
        previewImage->mark_dirty();

        if (xform) {
            cmsDeleteTransform(xform);
        }
    }
}


Cairo::RefPtr<Cairo::ImageSurface> PreviewImage::getImage()
{
    return previewImage;
}


Image8 *PreviewImage::load_img(const Glib::ustring &fname, int w, int h)
{
    StdImageSource imgSrc;
    if (imgSrc.load (fname)) {
        return nullptr;
    }

    ImageIO *img = imgSrc.getImageIO();

    if (w < 0) {
        w = img->getWidth();
        h = img->getHeight();
    } else {
        std::cout << "BOUNDING BOX: " << w << "x" << h << std::endl;
        // (w, h) is a bounding box
        double sw = std::max(double(img->getWidth()) / w, 1.0);
        double sh = std::max(double(img->getHeight()) / h, 1.0);
        if (sw > sh) {
            h = img->getHeight() / sw;
            w = img->getWidth() / sw;
        } else {
            w = img->getWidth() / sh;
            h = img->getHeight() / sh;
        }
    }

    Image8 *ret = new Image8(w, h);

    if (img->getType() == sImage8) {
        static_cast<Image8 *>(img)->resizeImgTo(w, h, TI_Bilinear, ret);
    } else if (img->getType() == sImage16) {
        static_cast<Image16 *>(img)->resizeImgTo(w, h, TI_Bilinear, ret);
    } else if (img->getType() == sImagefloat) {
        static_cast<Imagefloat *>(img)->resizeImgTo(w, h, TI_Bilinear, ret);
    } else {
        delete ret;
        ret = nullptr;
    }
    
    return ret;
}


Image8 *PreviewImage::load_raw_preview(const Glib::ustring &fname, int w, int h)
{
    RawImage ri(fname);
    unsigned int imageNum = 0;
    int r = ri.loadRaw(false, imageNum, false);

    if (r) {
        return nullptr;
    }

    int err = 1;

    // See if it is something we support
    if (!ri.checkThumbOk()) {
        return nullptr;
    }

    Image8 *img = new Image8();
    // No sample format detection occurred earlier, so we set them here,
    // as they are mandatory for the setScanline method
    img->setSampleFormat(IIOSF_UNSIGNED_CHAR);
    img->setSampleArrangement(IIOSA_CHUNKY);

    const char *data ((const char*)fdata (ri.get_thumbOffset(), ri.get_file()));

    if ((unsigned char)data[1] == 0xd8) {
        err = img->loadJPEGFromMemory(data, ri.get_thumbLength());
    } else if (ri.is_ppmThumb()) {
        err = img->loadPPMFromMemory(data, ri.get_thumbWidth(), ri.get_thumbHeight(), ri.get_thumbSwap(), ri.get_thumbBPS());
    }

    // did we succeed?
    if (err) {
        delete img;
        return nullptr;
    }

    if (w > 0 && h > 0) {
        double fw = img->getWidth();
        double fh = img->getHeight();
        if ((ri.get_rotateDegree() == 90 || ri.get_rotateDegree() == 270) && ri.thumbNeedsRotation()) {
            std::swap(w, h);
        }
        double sw = std::max(fw / w, 1.0);
        double sh = std::max(fh / h, 1.0);
        if (sw > sh) {
            h = fh / sw;
            w = fw / sw;
        } else {
            w = fw / sh;
            h = fh / sh;
        }
        if ((w != fw || h != fh) && w <= fw && h <= fh) {
            Image8 *img2 = new Image8(w, h);
            img->resizeImgTo(w, h, TI_Bilinear, img2);
            delete img;
            img = img2;
        }
    }

    if (ri.get_rotateDegree() > 0 && ri.thumbNeedsRotation()) {
        img->rotate(ri.get_rotateDegree());
    }

    return img;
}


namespace {

class PreviewRawImageSource: public RawImageSource {
public:
    PreviewRawImageSource(int bw, int bh):
        RawImageSource(),
        bbox_W_(bw),
        bbox_H_(bh)
    {
    }

    void rescale()
    {
        if (bbox_W_ < 0 || bbox_H_ < 0) {
            return;
        }
        if (fuji || d1x) {
            return;
        }
        if (numFrames != 1) {
            return;
        }

        // (w, h) is a bounding box
        double sw = std::max(double(W) / bbox_W_, 1.0);
        double sh = std::max(double(H) / bbox_H_, 1.0);
        int skip = std::max(sw, sh);
        
        if (skip <= 1) {
            return;
        }

        if (ri->getSensorType() == ST_BAYER) {
            skip /= 2;
            if (skip & 1) {
                --skip;
            }
            if (skip <= 1) {
                return;
            }

            if (settings->verbose) {
                std::cout << "SKIP: " << skip << ", FROM: " << W << "x" << H
                          << " to " << (W/skip) << "x" << (H/skip) << std::endl;
            }
            
            array2D<float> tmp;
            tmp(W, H, static_cast<float *>(rawData), 0);
            W /= skip;
            H /= skip;
            rawData.free();
            rawData(W, H);
            
#ifdef _OPENMP
#           pragma omp parallel for
#endif
            for (int y = 0; y < H; ++y) {
                int yy = y * skip + int(y & 1);
                for (int x = 0; x < W; ++x) {
                    int xx = x * skip + int(x & 1);
                    rawData[y][x] = tmp[yy][xx];
                }
            }
            flushRGB();
            red(W, H);
            green(W, H);
            blue(W, H);
        } else if (ri->getSensorType() == ST_FUJI_XTRANS) {
            return; // TODO
        } else {
            int c = ri->get_colors();
            if (c > 3) {
                return;
            }

            if (settings->verbose) {
                std::cout << "SKIP: " << skip << ", FROM: " << W << "x" << H
                          << " to " << (W/skip) << "x" << (H/skip) << std::endl;
            }
            
            array2D<float> tmp;
            tmp(W * c, H, static_cast<float *>(rawData), 0);
            W /= skip;
            H /= skip;
            rawData.free();
            rawData(W * c, H);
            
#ifdef _OPENMP
#           pragma omp parallel for
#endif
            for (int y = 0; y < H; ++y) {
                int yy = y * skip;
                for (int x = 0; x < W; ++x) {
                    int xx = x * skip;
                    for (int j = 0; j < c; ++j) {
                        rawData[y][c * x + j] = tmp[yy][c * xx + j];
                    }
                }
            }
            flushRGB();
            red(W, H);
            green(W, H);
            blue(W, H);
        }
    }

    bool mark_clipped()
    {
        array2D<float> *channels[3] = { &red, &green, &blue };

        constexpr float hl = 65534.f;
        constexpr float sh = hl / 2.f;
        
        if ((ri->getSensorType() == ST_BAYER || ri->getSensorType() == ST_FUJI_XTRANS) && ri->get_colors() == 3) {
            const bool bayer = ri->getSensorType() == ST_BAYER;
#ifdef _OPENMP
#           pragma omp parallel for
#endif
            for (int y = 0; y < H; ++y) {
                for (int x = 0; x < W; ++x) {
                    int c = bayer ? ri->FC(y, x) : ri->XTRANSFC(y, x);
                    if (bayer && c == 3) {
                        c = 1;
                    }
                    if (rawData[y][x] >= clmax[c]) {
                        for (int i = 0; i < 3; ++i) {
                            (*channels[i])[y][x] = 0.f;
                        }
                        (*channels[c])[y][x] = hl;
                    } else if (rawData[y][x] <= 0.f) {
                        for (int i = 0; i < 3; ++i) {
                            (*channels[i])[y][x] = 0.f;
                        }
                        (*channels[c])[y][x] = sh;
                    }     
                }
            }
            return true;
        } else if (ri->get_colors() == 1) {
            const float chmax = 65535.f;
            
#ifdef _OPENMP
#           pragma omp parallel for
#endif
            for (int y = 0; y < H; ++y) {
                for (int xx = 0; xx < W; ++xx) {
                    int x = 3 * xx;
                    bool ch = true;
                    bool cl = true;
                    for (int i = 0; i < 3; ++i) {
                        if ((*channels[i])[y][x] < chmax) {
                            ch = false;
                        }
                        if ((*channels[i])[y][x] > 0.f) {
                            cl = false;
                        }
                    }
                    if (ch) {
                        (*channels[0])[y][x] = (*channels[1])[y][x] = hl;
                        (*channels[2])[y][x] = 0.f;
                    } else if (cl) {
                        (*channels[0])[y][x] = (*channels[1])[y][x] = sh;
                        (*channels[2])[y][x] = 0.f;
                    }
                }
            }
            return true;
        }
        return false;
    }


    void mark_clipped_rgb(Imagefloat *img)
    {
        const int iw = img->getWidth();
        const int ih = img->getHeight();
        
#ifdef _OPENMP
#       pragma omp parallel for
#endif
        for (int y = 0; y < ih; ++y) {
            for (int x = 0; x < iw; ++x) {
                float &r = img->r(y, x);
                float &g = img->g(y, x);
                float &b = img->b(y, x);
                if (r >= MAXVALF && b >= MAXVALF && b >= MAXVALF) {
                    r = g = 65534.f;
                    b = 0.f;
                } else if (r <= 0.f && g <= 0.f && b <= 0.f) {
                    r = g = 65534.f / 2.f;
                    b = 0.f;
                } else {
                    r = g = b = Color::rgbLuminance(r, g, b);
                }
            }
        }
    }

private:
    int bbox_W_;
    int bbox_H_;
};

} // namespace


Image8 *PreviewImage::load_raw(const Glib::ustring &fname, int w, int h)
{
    PreviewRawImageSource src(w, h);
    int err = src.load(fname, true);
    if (err) {
        return nullptr;
    }

    const bool show_clip = settings->thumbnail_inspector_raw_curve == Settings::ThumbnailInspectorRawCurve::RAW_CLIPPING;

    ProcParams neutral;
    neutral.raw.bayersensor.method = RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::FAST);
    neutral.raw.xtranssensor.method = RAWParams::XTransSensor::getMethodString(RAWParams::XTransSensor::Method::FAST);
    neutral.icm.inputProfile = "(camera)";
    neutral.icm.workingProfile = options.rtSettings.srgb;

    if (show_clip) {
        neutral.raw.bayersensor.method = RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::MONO);
        neutral.raw.xtranssensor.method = RAWParams::XTransSensor::getMethodString(RAWParams::XTransSensor::Method::MONO);
    }

    src.preprocess(neutral.raw, neutral.lensProf, neutral.coarse, false);
    double thresholdDummy = 0.f;

    src.rescale();
    src.demosaic(neutral.raw, false, thresholdDummy);

    bool marked = show_clip && src.mark_clipped();
    // if (show_clip) {
    //     src.mark_clipped();
    // }

    int fw, fh;
    src.getFullSize(fw, fh);
    double scale = 1.0;

    if (w < 0) {
        w = fw;
        h = fh;
        scale = 1.0;
    } else {
        // (w, h) is a bounding box
        double sw = std::max(double(fw) / w, 1.0);
        double sh = std::max(double(fh) / h, 1.0);
        if (sw > sh) {
            scale = sw;
        } else {
            scale = sh;
        }
        w = fw / scale;
        h = fh / scale;
    }
        
    PreviewProps pp(0, 0, fw, fh, int(scale));//1);
    int iw = fw / int(scale);
    int ih = fh / int(scale);

    Imagefloat tmp(iw, ih);
    ColorTemp wb = src.getWB();
    if (show_clip) {
        wb = ColorTemp();
    }
    src.getImage(wb, TR_NONE, &tmp, pp, neutral.exposure, neutral.raw);
    if (!show_clip) {
        src.convertColorSpace(&tmp, neutral.icm, wb);
    } else if (!marked) {
        src.mark_clipped_rgb(&tmp);
    }

    LUTi gamma(65536);
    const bool apply_curve = settings->thumbnail_inspector_raw_curve != Settings::ThumbnailInspectorRawCurve::LINEAR && !show_clip;
    if (apply_curve) {
        static const std::vector<double> shadowcurve = {
            DCT_CatumullRom,
            0, 0,
            0.1, 0.4,
            0.70, 0.75,
            1, 1
        };
        DiagonalCurve curve(settings->thumbnail_inspector_raw_curve == Settings::ThumbnailInspectorRawCurve::FILM ? curves::filmcurve_def : shadowcurve);
        for (int i = 0; i < 65536; ++i) {
            float x = Color::gamma_srgbclipped(i) / 65535.f;
            float y = curve.getVal(x) * 255.f;
            gamma[i] = y;
        }
    } else {
        for (int i = 0; i < 65536; ++i) {
            gamma[i] = int(Color::gamma_srgbclipped(i)) * 255 / 65535;
        }
    }

    Image8 *img = new Image8(iw, ih);
    const int maxval = MAXVALF;
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (int y = 0; y < ih; ++y) {
        for (int x = 0; x < iw; ++x) {
            int r = tmp.r(y, x);
            int g = tmp.g(y, x);
            int b = tmp.b(y, x);
            // avoid magenta highlights
            if (r >= maxval && b >= maxval) {
                int v = CLIP((r + g + b) / 3) * 255 / 65535;
                img->r(y, x) = img->g(y, x) = img->b(y, x) = v;
            } else {
                img->r(y, x) = gamma[r];
                img->g(y, x) = gamma[g];
                img->b(y, x) = gamma[b];
            }
        }
    }

    if (w != iw || h != ih) {
        Image8 *rimg = new Image8(w, h);
        img->resizeImgTo(w, h, TI_Bilinear, rimg);
        delete img;
        img = rimg;
    }

    return img;
}
