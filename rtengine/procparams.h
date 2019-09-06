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
#pragma once

#include <cmath>
#include <cstdio>
#include <type_traits>
#include <vector>

#include <glibmm.h>
#include <lcms2.h>

#include "noncopyable.h"
#include "../rtgui/paramsedited.h"

class ParamsEdited;

namespace rtengine {

enum RenderingIntent {
    RI_PERCEPTUAL = INTENT_PERCEPTUAL,
    RI_RELATIVE = INTENT_RELATIVE_COLORIMETRIC,
    RI_SATURATION = INTENT_SATURATION,
    RI_ABSOLUTE = INTENT_ABSOLUTE_COLORIMETRIC,
    RI__COUNT
};


namespace procparams {

class KeyFile {
public:
    explicit KeyFile(const Glib::ustring &prefix=""): prefix_(prefix) {}
    
    bool has_group(const Glib::ustring &grp) const;
    bool has_key(const Glib::ustring &grp, const Glib::ustring &key) const;
    Glib::ArrayHandle<Glib::ustring> get_keys(const Glib::ustring &grp) const;
    
    Glib::ustring get_string(const Glib::ustring &grp, const Glib::ustring &key) const;
    int get_integer(const Glib::ustring &grp, const Glib::ustring &key) const;
    double get_double(const Glib::ustring &grp, const Glib::ustring &key) const;
    bool get_boolean(const Glib::ustring &grp, const Glib::ustring &key) const;
    Glib::ArrayHandle<Glib::ustring> get_string_list(const Glib::ustring &grp, const Glib::ustring &key) const;
    Glib::ArrayHandle<int> get_integer_list(const Glib::ustring &grp, const Glib::ustring &key) const;
    Glib::ArrayHandle<double> get_double_list(const Glib::ustring &grp, const Glib::ustring &key) const;
    
    void set_string(const Glib::ustring &grp, const Glib::ustring &key, const Glib::ustring& string);
    void set_boolean(const Glib::ustring &grp, const Glib::ustring &key, bool value);
    void set_integer(const Glib::ustring &grp, const Glib::ustring &key, int value);
    void set_double(const Glib::ustring &grp, const Glib::ustring &key, double value);
    void set_string_list(const Glib::ustring &grp, const Glib::ustring &key, const Glib::ArrayHandle<Glib::ustring> &list);
    void set_integer_list(const Glib::ustring &grp, const Glib::ustring &key, const Glib::ArrayHandle<int> &list);
    void set_double_list(const Glib::ustring &grp, const Glib::ustring &key, const Glib::ArrayHandle<double> &list);

    bool load_from_file(const Glib::ustring &fn);
    Glib::ustring to_data();

    Glib::ustring get_prefix() const { return prefix_; }
    void set_prefix(const Glib::ustring &prefix) { prefix_ = prefix; }
    
private:
    Glib::ustring GRP(const Glib::ustring &g) const { return prefix_ + g; }
    
    Glib::ustring prefix_;
    Glib::KeyFile kf_;
};


struct AreaMask {
    struct Shape {
        enum Mode {
            ADD,
            SUBTRACT,
            INTERSECT
        };
        double x; // [-100,100], with 0 as center of the image
        double y; // [-100,100]
        double width; // [0,200], with 100 as image width
        double height; // [0,200]
        double angle; // in degrees
        double roundness; // [0,100] (0 = rectangle, 100 = ellipse)
        Mode mode;
        Shape();
        bool operator==(const Shape &other) const;
        bool operator!=(const Shape &other) const;
    };
    bool inverted;
    double feather; // [0,100]
    std::vector<double> contrast; // curve
    std::vector<Shape> shapes;
        
    AreaMask();
    bool operator==(const AreaMask &other) const;
    bool operator!=(const AreaMask &other) const;
    bool isTrivial() const;
};


class LabCorrectionMask {
public:
    std::vector<double> hueMask;
    std::vector<double> chromaticityMask;
    std::vector<double> lightnessMask;
    double maskBlur;
    bool areaEnabled;
    AreaMask areaMask;

    LabCorrectionMask();
    bool operator==(const LabCorrectionMask &other) const;
    bool operator!=(const LabCorrectionMask &other) const;

    bool load(const KeyFile &keyfile, const Glib::ustring &group_name,
              const Glib::ustring &prefix, const Glib::ustring &suffix);
    void save(KeyFile &keyfile, const Glib::ustring &group_name,
              const Glib::ustring &prefix, const Glib::ustring &suffix) const;
};





template<typename T>
class Threshold final
{
public:
    Threshold(T _bottom, T _top, bool _start_at_one) :
        Threshold(_bottom, _top, 0, 0, _start_at_one, false)
    {
    }

    Threshold(T _bottom_left, T _top_left, T _bottom_right, T _top_right, bool _start_at_one) :
        Threshold(_bottom_left, _top_left, _bottom_right, _top_right, _start_at_one, true)
    {
    }

    template<typename U = T>
    typename std::enable_if<std::is_floating_point<U>::value, bool>::type operator ==(const Threshold<U>& rhs) const
    {
        if (is_double) {
            return
                std::fabs(bottom_left - rhs.bottom_left) < 1e-10
                && std::fabs(top_left - rhs.top_left) < 1e-10
                && std::fabs(bottom_right - rhs.bottom_right) < 1e-10
                && std::fabs(top_right - rhs.top_right) < 1e-10;
        } else {
            return
                std::fabs(bottom_left - rhs.bottom_left) < 1e-10
                && std::fabs(top_left - rhs.top_left) < 1e-10;
        }
    }

    template<typename U = T>
    typename std::enable_if<std::is_integral<U>::value, bool>::type operator ==(const Threshold<U>& rhs) const
    {
        if (is_double) {
            return
                bottom_left == rhs.bottom_left
                && top_left == rhs.top_left
                && bottom_right == rhs.bottom_right
                && top_right == rhs.top_right;
        } else {
            return
                bottom_left == rhs.bottom_left
                && top_left == rhs.top_left;
        }
    }

    T getBottom() const
    {
        return bottom_left;
    }

    T getTop() const
    {
        return top_left;
    }

    T getBottomLeft() const
    {
        return bottom_left;
    }

    T getTopLeft() const
    {
        return top_left;
    }

    T getBottomRight() const
    {
        return bottom_right;
    }

    T getTopRight() const
    {
        return top_right;
    }

    void setValues(T bottom, T top)
    {
        bottom_left = bottom;
        top_left = top;
    }

    void setValues(T bottom_left, T top_left, T bottom_right, T top_right)
    {
        this->bottom_left = bottom_left;
        this->top_left = top_left;
        this->bottom_right = bottom_right;
        this->top_right = top_right;
    }

    bool isDouble() const
    {
        return is_double;
    }

    std::vector<T> toVector() const
    {
        if (is_double) {
            return {
                bottom_left,
                top_left,
                bottom_right,
                top_right
            };
        } else {
            return {
                bottom_left,
                top_left
            };
        }
    }

    // RT: Type of the returned value
    // RV: Type of the value on the X axis
    // RV2: Type of the maximum value on the Y axis
    template <typename RT, typename RV, typename RV2>
    RT multiply(RV x, RV2 y_max) const
    {
        const double val = x;

        if (init_eql) {
            if (is_double) {
                if (val == static_cast<double>(bottom_right) && static_cast<double>(bottom_right) == static_cast<double>(top_right)) {
                    // This handles the special case where the 2 right values are the same, then bottom one is sent back,
                    // useful if one wants to keep the bottom value even beyond the x max bound
                    return 0;
                }

                if (val >= static_cast<double>(top_right)) {
                    return y_max;
                }

                if (val > static_cast<double>(bottom_right)) {
                    return static_cast<double>(y_max * (val - static_cast<double>(bottom_right)) / (static_cast<double>(top_right) - static_cast<double>(bottom_right)));
                }
            }

            if (val >= static_cast<double>(bottom_left)) {
                return 0;
            }

            if (val > static_cast<double>(top_left)) {
                return static_cast<double>(y_max * (1. - (val - static_cast<double>(bottom_left)) / (static_cast<double>(top_left) - static_cast<double>(bottom_left))));
            }

            return y_max;
        } else {
            if (is_double) {
                if (val == static_cast<double>(bottom_right) && static_cast<double>(bottom_right) == static_cast<double>(top_right)) {
                    // This handles the special case where the 2 right values are the same, then top one is sent back,
                    // useful if one wants to keep the top value even beyond the x max bound
                    return y_max;
                }

                if (val >= static_cast<double>(bottom_right)) {
                    return 0;
                }

                if (val > static_cast<double>(top_right)) {
                    return static_cast<double>(y_max * (1.0 - (val - static_cast<double>(top_right)) / (static_cast<double>(bottom_right) - static_cast<double>(top_right))));
                }
            }

            if (val >= static_cast<double>(top_left)) {
                return y_max;
            }

            if (val > static_cast<double>(bottom_left)) {
                return static_cast<double>(y_max * (val - static_cast<double>(bottom_left)) / (static_cast<double>(top_left) - static_cast<double>(bottom_left)));
            }

            return 0;
        }
    }

private:
    Threshold(T _bottom_left, T _top_left, T _bottom_right, T _top_right, bool _start_at_one, bool _is_double) :
        bottom_left(_bottom_left),
        top_left(_top_left),
        bottom_right(_bottom_right),
        top_right(_top_right),
        init_eql(_start_at_one),
        is_double(_is_double)
    {
    }

    T bottom_left;
    T top_left;
    T bottom_right;
    T top_right;
    bool init_eql;
    bool is_double;
};


struct ExposureParams {
    bool enabled;
    bool autoexp;
    double clip;
    enum HighlightReconstruction {
        HR_OFF,
        HR_BLEND,
        HR_COLOR
    };
    HighlightReconstruction hrmode;
    double expcomp;
    int black;
    int shcompr;
    int hlcompr;        // Highlight Recovery's compression
    int hlcomprthresh;  // Highlight Recovery's threshold
    bool clampOOG; // clamp out of gamut colours

    ExposureParams();

    bool operator==(const ExposureParams &other) const;
    bool operator!=(const ExposureParams &other) const;
};


struct BrightnessContrastSaturationParams {
    bool enabled;
    int brightness;
    int contrast;
    int saturation;
    int vibrance;

    BrightnessContrastSaturationParams();

    bool operator ==(const BrightnessContrastSaturationParams &other) const;
    bool operator !=(const BrightnessContrastSaturationParams &other) const;
};

/**
  * Parameters of the tone curve
  */
struct ToneCurveParams {
    bool enabled;
    
    enum class TcMode {
        STD,               // Standard modes, the curve is applied on all component individually
        WEIGHTEDSTD,       // Weighted standard mode
        FILMLIKE,          // Film-like mode, as defined in Adobe's reference code
        SATANDVALBLENDING, // Modify the Saturation and Value channel
        LUMINANCE,         // Modify the Luminance channel with coefficients from Rec 709's
        PERCEPTUAL         // Keep color appearance constant using perceptual modeling
    };

    std::vector<double> curve;
    std::vector<double> curve2;
    TcMode curveMode;
    TcMode curveMode2;
    bool histmatching; // histogram matching
    bool fromHistMatching;
    std::vector<double> saturation;

    ToneCurveParams();

    bool operator==(const ToneCurveParams &other) const;
    bool operator!=(const ToneCurveParams &other) const;
};


/**
  * Parameters of the luminance curve
  */
struct LCurveParams {
    bool enabled;
    std::vector<double>   lcurve;
    std::vector<double>   acurve;
    std::vector<double>   bcurve;
    std::vector<double>   cccurve;
    std::vector<double>   chcurve;
    std::vector<double>   lhcurve;
    std::vector<double>   hhcurve;
    std::vector<double>   lccurve;
    std::vector<double>   clcurve;
    int     brightness;
    int     contrast;
    int     chromaticity;
    bool    avoidcolorshift;
    double  rstprotection;
    bool    lcredsk;

    LCurveParams();

    bool operator ==(const LCurveParams& other) const;
    bool operator !=(const LCurveParams& other) const;
};


/**
 * Parameters for local contrast
 */
struct LocalContrastParams {
    bool enabled;
    enum Mode {
        USM,
        WAVELETS
    };
    Mode mode;
    int radius;
    double amount;
    double darkness;
    double lightness;
    double contrast;
    std::vector<double> curve;

    LocalContrastParams();

    bool operator==(const LocalContrastParams &other) const;
    bool operator!=(const LocalContrastParams &other) const;
};


/**
  * Parameters of the RGB curves
  */
struct RGBCurvesParams {
    bool enabled;
    std::vector<double>   rcurve;
    std::vector<double>   gcurve;
    std::vector<double>   bcurve;

    RGBCurvesParams();

    bool operator ==(const RGBCurvesParams& other) const;
    bool operator !=(const RGBCurvesParams& other) const;
};


/**
  * Parameters of the sharpening
  */
struct SharpeningParams {
    bool           enabled;
    double         contrast;
    double         blurradius;
    double         radius;
    int            amount;
    Threshold<int> threshold;
    bool           edgesonly;
    double         edges_radius;
    int            edges_tolerance;
    bool           halocontrol;
    int            halocontrol_amount;
    Glib::ustring  method;
    int            deconvamount;
    double         deconvradius;
    int            deconviter;
    int            deconvdamping;

    SharpeningParams();

    bool operator ==(const SharpeningParams& other) const;
    bool operator !=(const SharpeningParams& other) const;
};


/**
  * Parameters of the white balance adjustments
  */
struct WBEntry {
    enum class Type {
        CAMERA,
        AUTO,
        DAYLIGHT,
        CLOUDY,
        SHADE,
        WATER,
        TUNGSTEN,
        FLUORESCENT,
        LAMP,
        FLASH,
        LED,
        // CUSTOM one must remain the last one!
        CUSTOM
    };

    Glib::ustring ppLabel;
    Type type;
    Glib::ustring GUILabel;
    int temperature;
    double green;
    double equal;
    double tempBias;
};

struct WBParams {
    bool enabled;
    Glib::ustring   method;
    int             temperature;
    double          green;
    double          equal;
    double          tempBias;

    WBParams();

    bool operator ==(const WBParams& other) const;
    bool operator !=(const WBParams& other) const;

    static const std::vector<WBEntry>& getWbEntries();
};


/**
 * Parameters of defringing
 */
struct DefringeParams {
    bool    enabled;
    double  radius;
    int     threshold;
    std::vector<double> huecurve;

    DefringeParams();

    bool operator ==(const DefringeParams& other) const;
    bool operator !=(const DefringeParams& other) const;
};

/**
  * Parameters of impulse denoising
  */
struct ImpulseDenoiseParams {
    bool    enabled;
    int     thresh;

    ImpulseDenoiseParams();

    bool operator ==(const ImpulseDenoiseParams& other) const;
    bool operator !=(const ImpulseDenoiseParams& other) const;
};

/**
 * Parameters of the directional pyramid denoising
 */
struct DenoiseParams {
    enum class ColorSpace {
        LAB,
        RGB
    };
    
    enum class ChrominanceMethod {
        MANUAL,
        AUTOMATIC
    };

    enum class MedianType {
        TYPE_3X3_SOFT,
        TYPE_3X3_STRONG,
        TYPE_5X5_SOFT,
        TYPE_5X5_STRONG,
        TYPE_7X7,
        TYPE_9X9
    };

    enum class MedianMethod {
        LUMINANCE,
        CHROMINANCE,
        LAB_WEIGHTED,
        LAB,
        RGB
    };

    enum class SmoothingMethod {
        MEDIAN,
        GUIDED
    };
    
    bool enabled;

    ColorSpace colorSpace;

    bool aggressive;
    double gamma;

    double luminance;
    double luminanceDetail;
    int luminanceDetailThreshold;

    ChrominanceMethod chrominanceMethod;
    double chrominanceAutoFactor;
    double chrominance;
    double chrominanceRedGreen;
    double chrominanceBlueYellow;

    bool smoothingEnabled;
    SmoothingMethod smoothingMethod;

    MedianType medianType;
    MedianMethod medianMethod;
    int medianIterations;

    int guidedLumaRadius;
    int guidedChromaRadius;
    int guidedLumaStrength;
    int guidedChromaStrength;

    DenoiseParams();

    bool operator ==(const DenoiseParams& other) const;
    bool operator !=(const DenoiseParams& other) const;
};

// EPD related parameters.
struct TextureBoostParams {
    struct Region {
        double strength;
        double edgeStopping;
        double scale;

        Region();
        bool operator==(const Region &other) const;
        bool operator!=(const Region &other) const;
    };

    bool   enabled;
    std::vector<Region> regions;
    std::vector<LabCorrectionMask> labmasks;
    int showMask;

    TextureBoostParams();

    bool operator ==(const TextureBoostParams& other) const;
    bool operator !=(const TextureBoostParams& other) const;
};


struct LogEncodingParams {
    bool enabled;
    bool autocompute;
    bool autogray;
    double sourceGray;
    double targetGray;
    double blackEv;
    double whiteEv;
    int detail;

    LogEncodingParams();

    bool operator==(const LogEncodingParams &other) const;
    bool operator !=(const LogEncodingParams &other) const;
};


struct FattalToneMappingParams {
    bool enabled;
    int threshold;
    int amount;

    FattalToneMappingParams();

    bool operator ==(const FattalToneMappingParams& other) const;
    bool operator !=(const FattalToneMappingParams& other) const;
};

/**
  * Parameters of the shadow/highlight enhancement
  */
struct SHParams {
    bool    enabled;
    int     highlights;
    int     htonalwidth;
    int     shadows;
    int     stonalwidth;
    int     radius;
    bool    lab;

    SHParams();

    bool operator ==(const SHParams& other) const;
    bool operator !=(const SHParams& other) const;
};


struct ToneEqualizerParams {
    bool enabled;
    std::array<int, 5> bands;
    int detail;

    ToneEqualizerParams();

    bool operator==(const ToneEqualizerParams &other) const;
    bool operator!=(const ToneEqualizerParams &other) const;
};

/**
  * Parameters of the cropping
  */
struct CropParams {
    bool    enabled;
    int     x;
    int     y;
    int     w;
    int     h;
    bool    fixratio;
    Glib::ustring   ratio;
    Glib::ustring   orientation;
    Glib::ustring   guide;

    CropParams();

    bool operator ==(const CropParams& other) const;
    bool operator !=(const CropParams& other) const;

    void mapToResized(int resizedWidth, int resizedHeight, int scale, int& x1, int& x2, int& y1, int& y2) const;
};

/**
  * Parameters of the coarse transformations like 90 deg rotations and h/v flipping
  */
struct CoarseTransformParams {
    int     rotate;
    bool    hflip;
    bool    vflip;

    CoarseTransformParams();

    bool operator ==(const CoarseTransformParams& other) const;
    bool operator !=(const CoarseTransformParams& other) const;
};

/**
  * Common transformation parameters
  */
struct CommonTransformParams {
    bool autofill;

    CommonTransformParams();

    bool operator ==(const CommonTransformParams& other) const;
    bool operator !=(const CommonTransformParams& other) const;
};

/**
  * Parameters of the rotation
  */
struct RotateParams {
    bool enabled;
    double  degree;

    RotateParams();

    bool operator ==(const RotateParams& other) const;
    bool operator !=(const RotateParams& other) const;
};

/**
  * Parameters of the distortion correction
  */
struct DistortionParams {
    bool enabled;
    double  amount;

    DistortionParams();

    bool operator ==(const DistortionParams& other) const;
    bool operator !=(const DistortionParams& other) const;
};

// Lens profile correction parameters
struct LensProfParams {
    enum class LcMode {
        NONE,               // No lens correction
        LENSFUNAUTOMATCH,   // Lens correction using auto matched lensfun database entry
        LENSFUNMANUAL,      // Lens correction using manually selected lensfun database entry
        LCP                 // Lens correction using lcp file
    };

    LcMode lcMode;
    Glib::ustring lcpFile;
    bool useDist, useVign, useCA;
    Glib::ustring lfCameraMake;
    Glib::ustring lfCameraModel;
    Glib::ustring lfLens;

    LensProfParams();

    bool operator ==(const LensProfParams& other) const;
    bool operator !=(const LensProfParams& other) const;

    bool useLensfun() const;
    bool lfAutoMatch() const;
    bool useLcp() const;
    bool lfManual() const;

    const std::vector<const char*>& getMethodStrings() const;
    Glib::ustring getMethodString(LcMode mode) const;
    LcMode getMethodNumber(const Glib::ustring& mode) const;
};


/**
  * Parameters of the perspective correction
  */
struct PerspectiveParams {
    bool enabled;
    double  horizontal;
    double  vertical;
    double angle;
    double shear;
    double flength;
    double cropfactor;
    double aspect;

    PerspectiveParams();

    bool operator ==(const PerspectiveParams& other) const;
    bool operator !=(const PerspectiveParams& other) const;
};

/**
  * Parameters of the gradient filter
  */
struct GradientParams {
    bool   enabled;
    double degree;
    int    feather;
    double strength;
    int    centerX;
    int    centerY;

    GradientParams();

    bool operator ==(const GradientParams& other) const;
    bool operator !=(const GradientParams& other) const;
};

/**
  * Parameters of the post-crop vignette filter
  */
struct PCVignetteParams {
    bool   enabled;
    double strength;
    int    feather;
    int    roundness;

    PCVignetteParams();

    bool operator ==(const PCVignetteParams& other) const;
    bool operator !=(const PCVignetteParams& other) const;
};

/**
  * Parameters of the vignetting correction
  */
struct VignettingParams {
    bool enabled;
    int  amount;
    int  radius;
    int  strength;
    int  centerX;
    int  centerY;

    VignettingParams();

    bool operator ==(const VignettingParams& other) const;
    bool operator !=(const VignettingParams& other) const;
};

/**
  * Parameters of the color mixer
  */
struct ChannelMixerParams {
    bool enabled;
    int red[3];
    int green[3];
    int blue[3];

    ChannelMixerParams();

    bool operator ==(const ChannelMixerParams& other) const;
    bool operator !=(const ChannelMixerParams& other) const;
};

struct BlackWhiteParams {
    bool enabled;

    Glib::ustring filter;
    Glib::ustring setting;
    int mixerRed;
    int mixerGreen;
    int mixerBlue;
    int gammaRed;
    int gammaGreen;
    int gammaBlue;

    BlackWhiteParams();

    bool operator ==(const BlackWhiteParams& other) const;
    bool operator !=(const BlackWhiteParams& other) const;
};

struct HSLEqualizerParams {
    bool enabled;
    std::vector<double> hCurve;
    std::vector<double> sCurve;
    std::vector<double> lCurve;
    int smoothing;

    HSLEqualizerParams();

    bool operator==(const HSLEqualizerParams &other) const;
    bool operator!=(const HSLEqualizerParams &other) const;
};

/**
  * Parameters of the c/a correction
  */
struct CACorrParams {
    bool enabled;
    double red;
    double blue;

    CACorrParams();

    bool operator ==(const CACorrParams& other) const;
    bool operator !=(const CACorrParams& other) const;
};

/**
  * Parameters of the resizing
  */
struct ResizeParams {
    bool enabled;
    double scale;
    Glib::ustring appliesTo;
    Glib::ustring method;
    int dataspec;
    int width;
    int height;
    bool allowUpscaling;

    ResizeParams();

    bool operator ==(const ResizeParams& other) const;
    bool operator !=(const ResizeParams& other) const;
};

/**
  * Parameters of the color spaces used during the processing
  */
struct ColorManagementParams {
    Glib::ustring inputProfile;
    bool toneCurve;
    bool applyLookTable;
    bool applyBaselineExposureOffset;
    bool applyHueSatMap;
    int dcpIlluminant;

    Glib::ustring workingProfile;

    Glib::ustring outputProfile;
    RenderingIntent outputIntent;
    bool outputBPC;

    static const Glib::ustring NoICMString;

    ColorManagementParams();

    bool operator ==(const ColorManagementParams& other) const;
    bool operator !=(const ColorManagementParams& other) const;
};


/**
  * Parameters for metadata handling
  */
struct MetaDataParams {
    enum Mode {
        TUNNEL,
        EDIT,
        STRIP
    };
    Mode mode;

    MetaDataParams();

    bool operator ==(const MetaDataParams &other) const;
    bool operator !=(const MetaDataParams &other) const;
};


/**
  * Typedef for representing a key/value for the exif metadata information
  */
typedef std::map<Glib::ustring, Glib::ustring> ExifPairs;

/**
  * The IPTC key/value pairs
  */
typedef std::map<Glib::ustring, std::vector<Glib::ustring>> IPTCPairs;

/**
* Directional pyramid equalizer params
*/
struct DirPyrEqualizerParams {
    bool enabled;
    // bool gamutlab;
    struct Levels {
        double mult[6];
        double threshold;
        Levels();
        bool operator==(const Levels &other) const;
        bool operator!=(const Levels &other) const;
    };
    std::vector<Levels> levels;
    std::vector<LabCorrectionMask> labmasks;
    int showMask;
    // double skinprotect;
    // Threshold<int> hueskin;
    // Glib::ustring cbdlMethod;

    DirPyrEqualizerParams();

    bool operator ==(const DirPyrEqualizerParams& other) const;
    bool operator !=(const DirPyrEqualizerParams& other) const;
};


/**
 *  Film simualtion params
 */
struct FilmSimulationParams {
    bool enabled;
    Glib::ustring clutFilename;
    int strength;

    FilmSimulationParams();

    bool operator ==(const FilmSimulationParams& other) const;
    bool operator !=(const FilmSimulationParams& other) const;
};


struct SoftLightParams {
    bool enabled;
    int strength;

    SoftLightParams();

    bool operator==(const SoftLightParams &other) const;
    bool operator!=(const SoftLightParams &other) const;
};


struct DehazeParams {
    bool enabled;
    int strength;
    bool showDepthMap;
    int depth;
    bool luminance;

    DehazeParams();

    bool operator==(const DehazeParams &other) const;
    bool operator!=(const DehazeParams &other) const;
};


struct GrainParams {
    bool enabled;
    int iso;
    int strength;
    int scale;

    GrainParams();

    bool operator==(const GrainParams &other) const;
    bool operator!=(const GrainParams &other) const;
};


struct GuidedSmoothingParams {
    struct Region {
        enum class Channel {
            LUMINANCE,
            CHROMINANCE,
            RGB
        };
        Channel channel;
        int radius;
        int epsilon;

        Region();
        bool operator==(const Region &other) const;
        bool operator!=(const Region &other) const;
    };
    bool enabled;
    std::vector<Region> regions;
    std::vector<LabCorrectionMask> labmasks;
    int showMask;

    GuidedSmoothingParams();

    bool operator==(const GuidedSmoothingParams &other) const;
    bool operator!=(const GuidedSmoothingParams &other) const;
};


struct ColorCorrectionParams {
    struct Region {
        enum { CHAN_ALL = -1, CHAN_R, CHAN_G, CHAN_B };
        double a;
        double b;
        double saturation;
        double slope;
        double offset;
        double power;
        double pivot;
        int channel;

        Region();
        bool operator==(const Region &other) const;
        bool operator!=(const Region &other) const;
    };

    bool enabled;
    std::vector<Region> regions;
    std::vector<LabCorrectionMask> labmasks;
    int showMask;

    ColorCorrectionParams();
    bool operator==(const ColorCorrectionParams &other) const;
    bool operator!=(const ColorCorrectionParams &other) const;
};


/**
  * Parameters for RAW demosaicing, common to all sensor type
  */
struct RAWParams {
    /**
     * Parameters for RAW demosaicing specific to Bayer sensors
     */
    struct BayerSensor {
        enum class Method {
            AMAZE,
            AMAZEVNG4,
            RCD,
            RCDVNG4,
            DCB,
            DCBVNG4,
            LMMSE,
            IGV,
            AHD,
            EAHD,
            HPHD,
            VNG4,
            FAST,
            MONO,
            PIXELSHIFT,
            NONE
        };

        enum class PSMotionCorrectionMethod {
            OFF,
            AUTO,
            CUSTOM
        };

        enum class PSDemosaicMethod {
            AMAZE,
            AMAZEVNG4,
            LMMSE
        };

        Glib::ustring method;
        int border;
        int imageNum;
        int ccSteps;
        double black0;
        double black1;
        double black2;
        double black3;
        bool twogreen;
        int linenoise;
        enum class LineNoiseDirection {
            HORIZONTAL = 1,
            VERTICAL,
            BOTH,
            PDAF_LINES = 5
        };
        LineNoiseDirection linenoiseDirection;
        int greenthresh;
        int dcb_iterations;
        int lmmse_iterations;
        bool dualDemosaicAutoContrast;
        double dualDemosaicContrast;
        PSMotionCorrectionMethod pixelShiftMotionCorrectionMethod;
        double pixelShiftEperIso;
        double pixelShiftSigma;
        bool pixelShiftShowMotion;
        bool pixelShiftShowMotionMaskOnly;
        bool pixelShiftHoleFill;
        bool pixelShiftMedian;
        bool pixelShiftGreen;
        bool pixelShiftBlur;
        double pixelShiftSmoothFactor;
        bool pixelShiftEqualBright;
        bool pixelShiftEqualBrightChannel;
        bool pixelShiftNonGreenCross;
        Glib::ustring pixelShiftDemosaicMethod;
        bool dcb_enhance;
        bool pdafLinesFilter;

        // some enable flags
        bool enable_black;
        bool enable_preproc;

        BayerSensor();

        bool operator ==(const BayerSensor& other) const;
        bool operator !=(const BayerSensor& other) const;

        void setPixelShiftDefaults();

        static const std::vector<const char*>& getMethodStrings();
        static Glib::ustring getMethodString(Method method);

        static const std::vector<const char*>& getPSDemosaicMethodStrings();
        static Glib::ustring getPSDemosaicMethodString(PSDemosaicMethod method);
    };

    /**
     * Parameters for RAW demosaicing specific to X-Trans sensors
     */
    struct XTransSensor {
        enum class Method {
            FOUR_PASS,
            THREE_PASS,
            TWO_PASS,
            ONE_PASS,
            FAST,
            MONO,
            NONE
        };

        Glib::ustring method;
        bool dualDemosaicAutoContrast;
        double dualDemosaicContrast;
        int border;
        int ccSteps;
        double blackred;
        double blackgreen;
        double blackblue;

        bool enable_black;

        XTransSensor();

        bool operator ==(const XTransSensor& other) const;
        bool operator !=(const XTransSensor& other) const;

        static const std::vector<const char*>& getMethodStrings();
        static Glib::ustring getMethodString(Method method);
    };

    BayerSensor bayersensor;         ///< RAW parameters for Bayer sensors
    XTransSensor xtranssensor;       ///< RAW parameters for X-Trans sensors

    enum class FlatFieldBlurType {
        AREA,
        V,
        H,
        VH,
    };

    Glib::ustring dark_frame;
    bool df_autoselect;

    Glib::ustring ff_file;
    bool ff_AutoSelect;
    int ff_BlurRadius;
    Glib::ustring ff_BlurType;
    bool ff_AutoClipControl;
    int ff_clipControl;

    bool ca_autocorrect;
    bool ca_avoidcolourshift;
    int caautoiterations;
    double cared;
    double cablue;

    // exposure before interpolation
    double expos;

    bool hotPixelFilter;
    bool deadPixelFilter;
    int hotdeadpix_thresh;

    // some enable flags
    bool enable_darkframe;
    bool enable_flatfield;
    bool enable_ca;
    bool enable_hotdeadpix;
    bool enable_whitepoint;

    RAWParams();

    bool operator ==(const RAWParams& other) const;
    bool operator !=(const RAWParams& other) const;

    static const std::vector<const char*>& getFlatFieldBlurTypeStrings();
    static Glib::ustring getFlatFieldBlurTypeString(FlatFieldBlurType type);
};


/**
  * This class holds all the processing parameters applied on the images
  */
class ProcParams {
public:
    ExposureParams          exposure;
    BrightnessContrastSaturationParams brightContrSat;
    ToneCurveParams         toneCurve;       ///< Tone curve parameters
    LCurveParams            labCurve;        ///< CIELAB luminance curve parameters
    LocalContrastParams     localContrast;   ////< Local contrast parameters
    RGBCurvesParams         rgbCurves;       ///< RGB curves parameters
    SharpeningParams        sharpening;      ///< Sharpening parameters
    SharpeningParams        prsharpening;    ///< Sharpening parameters for post resize sharpening
    WBParams                wb;              ///< White balance parameters
    DefringeParams          defringe;        ///< Defringing parameters
    ImpulseDenoiseParams    impulseDenoise;  ///< Impulse denoising parameters
    DenoiseParams           denoise;   ///< Directional Pyramid denoising parameters
    TextureBoostParams      textureBoost;  ///< Edge Preserving Decomposition parameters
    FattalToneMappingParams fattal;          ///< Dynamic Range Compression
    LogEncodingParams       logenc;
    SHParams                sh;              ///< Shadow/highlight enhancement parameters
    ToneEqualizerParams     toneEqualizer;
    CropParams              crop;            ///< Crop parameters
    CoarseTransformParams   coarse;          ///< Coarse transformation (90, 180, 270 deg rotation, h/v flipping) parameters
    CommonTransformParams   commonTrans;     ///< Common transformation parameters (autofill)
    RotateParams            rotate;          ///< Rotation parameters
    DistortionParams        distortion;      ///< Lens distortion correction parameters
    LensProfParams          lensProf;        ///< Lens correction profile parameters
    PerspectiveParams       perspective;     ///< Perspective correction parameters
    GradientParams          gradient;        ///< Gradient filter parameters
    PCVignetteParams        pcvignette;      ///< Post-crop vignette filter parameters
    CACorrParams            cacorrection;    ///< Lens c/a correction parameters
    VignettingParams        vignetting;      ///< Lens vignetting correction parameters
    ChannelMixerParams      chmixer;         ///< Channel mixer parameters
    BlackWhiteParams        blackwhite;      ///< Black&  White parameters
    HSLEqualizerParams      hsl;
    ResizeParams            resize;          ///< Resize parameters
    ColorManagementParams   icm;             ///< profiles/color spaces used during the image processing
    RAWParams               raw;             ///< RAW parameters before demosaicing
    DirPyrEqualizerParams   dirpyrequalizer; ///< directional pyramid wavelet parameters
    FilmSimulationParams    filmSimulation;  ///< film simulation parameters
    SoftLightParams         softlight;       ///< softlight parameters
    DehazeParams            dehaze;          ///< dehaze parameters
    GrainParams             grain;
    GuidedSmoothingParams   smoothing;
    ColorCorrectionParams   colorcorrection;
    int                     rank;            ///< Custom image quality ranking
    int                     colorlabel;      ///< Custom color label
    bool                    inTrash;         ///< Marks deleted image
    Glib::ustring           appVersion;      ///< Version of the application that generated the parameters
    int                     ppVersion;       ///< Version of the PP file from which the parameters have been read

    MetaDataParams          metadata;        ///< Metadata parameters
    ExifPairs               exif;            ///< List of modifications appplied on the exif tags of the input image
    IPTCPairs               iptc;            ///< The IPTC tags and values to be saved to the output image

    /**
      * The constructor only sets the hand-wired defaults.
      */
    ProcParams();
    /**
      * Sets the hand-wired defaults parameters.
      */
    void setDefaults();
    /**
      * Saves the parameters to possibly two files. This is a performance improvement if a function has to
      * save the same file in two different location, i.e. the cache and the image's directory
      * @param fname   the name of the first file (can be an empty string)
      * @param fname2  the name of the second file (can be an empty string) (optional)
      * @param fnameAbsolute set to false if embedded filenames (if any, darkframe/flatfield) should be stored as relative
      * filenames if they are inside the same directory or in a sub-directory to fname's directory.
      * @param pedited pointer to a ParamsEdited object (optional) to store which values has to be saved
      * @return Error code (=0 if all supplied filenames where created correctly)
      */
    int save(const Glib::ustring& fname, const Glib::ustring& fname2 = Glib::ustring(), bool fnameAbsolute = true, const ParamsEdited *pedited = nullptr);
    /**
      * Loads the parameters from a file.
      * @param fname the name of the file
      * @params pedited pointer to a ParamsEdited object (optional) to store which values has been loaded
      * @return Error code (=0 if no error)
      */
    int load(const Glib::ustring& fname, const ParamsEdited *pedited=nullptr);

    int load(const KeyFile &keyFile, const ParamsEdited *pedited=nullptr,
             bool resetOnError=true, const Glib::ustring &fname="");
    int save(KeyFile &keyFile, const ParamsEdited *pedited=nullptr,
             const Glib::ustring &fname="", bool fnameAbsolute=true) const;

    /** Creates a new instance of ProcParams.
      * @return a pointer to the new ProcParams instance. */
    static ProcParams* create();

    /** Destroys an instance of ProcParams.
      * @param pp a pointer to the ProcParams instance to destroy. */
    static void destroy(ProcParams* pp);

    static void init();
    static void cleanup();

    bool operator ==(const ProcParams& other) const;
    bool operator !=(const ProcParams& other) const;

private:
    /** Write the ProcParams's text in the file of the given name.
    * @param fname the name of the file
    * @param content the text to write
    * @return Error code (=0 if no error)
    * */
    int write(const Glib::ustring& fname, const Glib::ustring& content) const;

    int load(bool load_general,
             const KeyFile &keyFile, const ParamsEdited *pedited,
             bool resetOnError, const Glib::ustring &fname);
    int save(bool save_general,
             KeyFile &keyFile, const ParamsEdited *pedited,
             const Glib::ustring &fname, bool fnameAbsolute) const;

    friend class ProcParamsWithSnapshots;
};


class ProcParamsWithSnapshots {
public:
    int load(const Glib::ustring &fname);
    int save(const Glib::ustring &fname, const Glib::ustring &fname2=Glib::ustring(), bool fnameAbsolute=true);

    ProcParams master;
    std::vector<std::pair<Glib::ustring, ProcParams>> snapshots;
};


class PartialProfile {
public:
    virtual ~PartialProfile() = default;
    virtual bool applyTo(ProcParams &pp) const = 0;
};


class FullPartialProfile: public PartialProfile {
public:
    explicit FullPartialProfile(const ProcParams &pp);
    bool applyTo(ProcParams &pp) const override;

private:
    ProcParams pp_;
};


class FilePartialProfile: public PartialProfile {
public:
    FilePartialProfile(): fname_(""), full_(true) {}
    explicit FilePartialProfile(const Glib::ustring &fname, bool full=false);
    bool applyTo(ProcParams &pp) const override;

private:
    Glib::ustring fname_;
    bool full_;
};


class PEditedPartialProfile: public PartialProfile {
public:
    PEditedPartialProfile(const Glib::ustring &fname, const ParamsEdited &pe);
    PEditedPartialProfile(const ProcParams &pp, const ParamsEdited &pe);
    bool applyTo(ProcParams &pp) const override;

private:
    Glib::ustring fname_;
    ProcParams pp_;
    ParamsEdited pe_;
};


}} // namespace rtengine::procparams

