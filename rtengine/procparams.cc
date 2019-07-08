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

#include <map>

#include <locale.h>

#include <glib/gstdio.h>

#include "curves.h"
#include "procparams.h"

#include "../rtgui/multilangmgr.h"
#include "../rtgui/options.h"
#include "../rtgui/paramsedited.h"
#include "../rtgui/ppversion.h"
#include "../rtgui/version.h"

using namespace std;

namespace
{

Glib::ustring expandRelativePath(const Glib::ustring &procparams_fname, const Glib::ustring &prefix, Glib::ustring embedded_fname)
{
    if (embedded_fname == "" || !Glib::path_is_absolute(procparams_fname)) {
        return embedded_fname;
    }

    if (prefix != "") {
        if (embedded_fname.length() < prefix.length() || embedded_fname.substr(0, prefix.length()) != prefix) {
            return embedded_fname;
        }

        embedded_fname = embedded_fname.substr(prefix.length());
    }

    if (Glib::path_is_absolute(embedded_fname)) {
        return prefix + embedded_fname;
    }

    Glib::ustring absPath = prefix + Glib::path_get_dirname(procparams_fname) + G_DIR_SEPARATOR_S + embedded_fname;
    return absPath;
}

Glib::ustring relativePathIfInside(const Glib::ustring &procparams_fname, bool fnameAbsolute, Glib::ustring embedded_fname)
{
    if (fnameAbsolute || embedded_fname == "" || !Glib::path_is_absolute(procparams_fname)) {
        return embedded_fname;
    }

    Glib::ustring prefix = "";

    if (embedded_fname.length() > 5 && embedded_fname.substr(0, 5) == "file:") {
        embedded_fname = embedded_fname.substr(5);
        prefix = "file:";
    }

    if (!Glib::path_is_absolute(embedded_fname)) {
        return prefix + embedded_fname;
    }

    Glib::ustring dir1 = Glib::path_get_dirname(procparams_fname) + G_DIR_SEPARATOR_S;
    Glib::ustring dir2 = Glib::path_get_dirname(embedded_fname) + G_DIR_SEPARATOR_S;

    if (dir2.substr(0, dir1.length()) != dir1) {
        // it's in a different directory, ie not inside
        return prefix + embedded_fname;
    }

    return prefix + embedded_fname.substr(dir1.length());
}

void getFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    int& value
)
{
    value = keyfile.get_integer(group_name, key);
}

void getFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    double& value
)
{
    value = keyfile.get_double(group_name, key);
}

void getFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    bool& value
)
{
    try {
        value = keyfile.get_boolean(group_name, key);
    } catch (Glib::KeyFileError &e) {
        int v = keyfile.get_integer(group_name, key);
        value = v;
    }
}

void getFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    Glib::ustring& value
)
{
    value = keyfile.get_string(group_name, key);
}

void getFromKeyfile(
    const Glib::KeyFile& keyfile,
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    std::vector<double>& value
)
{
    value = keyfile.get_double_list(group_name, key);
    rtengine::sanitizeCurve(value);
}

template<typename T>
bool assignFromKeyfile(const Glib::KeyFile& keyfile, const Glib::ustring& group_name, const Glib::ustring& key, T &value)
{
    if (keyfile.has_key(group_name, key)) {
        getFromKeyfile(keyfile, group_name, key, value);

        return true;
    }

    return false;
}

template<typename T, typename = typename std::enable_if<std::is_enum<T>::value>::type>
bool assignFromKeyfile(const Glib::KeyFile& keyfile, const Glib::ustring& group_name, const Glib::ustring& key, const std::map<std::string, T>& mapping, T& value)
{
    if (keyfile.has_key(group_name, key)) {
        Glib::ustring v;
        getFromKeyfile(keyfile, group_name, key, v);

        const typename std::map<std::string, T>::const_iterator m = mapping.find(v);

        if (m != mapping.end()) {
            value = m->second;
        } else {
            return false;
        }

        return true;
    }

    return false;
}

void putToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    int value,
    Glib::KeyFile& keyfile
)
{
    keyfile.set_integer(group_name, key, value);
}

void putToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    double value,
    Glib::KeyFile& keyfile
)
{
    keyfile.set_double(group_name, key, value);
}

void putToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    bool value,
    Glib::KeyFile& keyfile
)
{
    keyfile.set_boolean(group_name, key, value);
}

void putToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    const Glib::ustring& value,
    Glib::KeyFile& keyfile
)
{
    keyfile.set_string(group_name, key, value);
}

void putToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    const std::vector<int>& value,
    Glib::KeyFile& keyfile
)
{
    const Glib::ArrayHandle<int> list = value;
    keyfile.set_integer_list(group_name, key, list);
}

void putToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    const std::vector<double>& value,
    Glib::KeyFile& keyfile
)
{
    const Glib::ArrayHandle<double> list = value;
    keyfile.set_double_list(group_name, key, list);
}

template<typename T>
bool saveToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    const T& value,
    Glib::KeyFile& keyfile
)
{
    putToKeyfile(group_name, key, value, keyfile);
    return true;
}

template<typename T, typename = typename std::enable_if<std::is_enum<T>::value>::type>
bool saveToKeyfile(
    const Glib::ustring& group_name,
    const Glib::ustring& key,
    const std::map<T, const char*>& mapping,
    const T& value,
    Glib::KeyFile& keyfile
)
{
    const typename std::map<T, const char*>::const_iterator m = mapping.find(value);

    if (m != mapping.end()) {
        keyfile.set_string(group_name, key, m->second);
        return true;
    }

    return false;
}


const std::map<Glib::ustring, Glib::ustring> exif_keys = {
    {"Copyright", "Exif.Image.Copyright"},
    {"Artist", "Exif.Image.Artist"},
    {"ImageDescription", "Exif.Image.ImageDescription"},
    {"Exif.UserComment", "Exif.Photo.UserComment"}
};

const std::map<Glib::ustring, Glib::ustring> iptc_keys = {
    {"Title", "Iptc.Application2.ObjectName"},
    {"Category", "Iptc.Application2.Category"},
    {"SupplementalCategories", "Iptc.Application2.SuppCategory"},
    {"Keywords", "Iptc.Application2.Keywords"},
    {"Instructions", "Iptc.Application2.SpecialInstructions"},
    {"DateCreated", "Iptc.Application2.DateCreated"},
    {"Creator", "Iptc.Application2.Byline"},
    {"CreatorJobTitle", "Iptc.Application2.BylineTitle"},
    {"City", "Iptc.Application2.City"},
    {"Province", "Iptc.Application2.ProvinceState"},
    {"Country", "Iptc.Application2.CountryName"},
    {"TransReference", "Iptc.Application2.TransmissionReference"},
    {"Headline", "Iptc.Application2.Headline"},
    {"Credit", "Iptc.Application2.Credit"},
    {"Source", "Iptc.Application2.Source"},
    {"Copyright", "Iptc.Application2.Copyright"},
    {"Caption", "Iptc.Application2.Caption"},
    {"CaptionWriter", "Iptc.Application2.Writer"}
};

} // namespace

namespace rtengine { namespace procparams {

AreaMask::Shape::Shape():
    x(0),
    y(0),
    width(100),
    height(100),
    angle(0),
    roundness(0)
{
}


bool AreaMask::Shape::operator==(const Shape &other) const
{
    return
        x == other.x
        && y == other.y
        && width == other.width
        && height == other.height
        && angle == other.angle
        && roundness == other.roundness;
}


bool AreaMask::Shape::operator!=(const Shape &other) const
{
    return !(*this == other);
}


AreaMask::AreaMask():
    inverted(false),
    feather(0),
    contrast{DCT_Linear},
    shapes{Shape()}
{
}


bool AreaMask::operator==(const AreaMask &other) const
{
    return inverted == other.inverted
        && feather == other.feather
        && contrast == other.contrast
        && shapes == other.shapes;
}


bool AreaMask::operator!=(const AreaMask &other) const
{
    return !(*this == other);
}


bool AreaMask::isTrivial() const
{
    return (*this == AreaMask());
}


LabCorrectionMask::LabCorrectionMask():
    hueMask{
        FCT_MinMaxCPoints,
            0.166666667,
            1.,
            0.35,
            0.35,
            0.8287775246,
            1.,
            0.35,
            0.35
    },
    chromaticityMask{
        FCT_MinMaxCPoints,
            0.,
            1.,
            0.35,
            0.35,
            1.,
            1.,
            0.35,
            0.35
            },
    lightnessMask{
        FCT_MinMaxCPoints,
            0.,
            1.,
            0.35,
            0.35,
            1.,
            1.,
            0.35,
            0.35
            },
    maskBlur(0),
    areaEnabled(false),
    areaMask()
{
}


bool LabCorrectionMask::operator==(const LabCorrectionMask &other) const
{
    return hueMask == other.hueMask
        && chromaticityMask == other.chromaticityMask
        && lightnessMask == other.lightnessMask
        && maskBlur == other.maskBlur
        && areaEnabled == other.areaEnabled
        && areaMask == other.areaMask;
}


bool LabCorrectionMask::operator!=(const LabCorrectionMask &other) const
{
    return !(*this == other);
}


bool LabCorrectionMask::load(const Glib::KeyFile &keyfile, const Glib::ustring &group_name, const Glib::ustring &prefix, const Glib::ustring &suffix)
{
    bool ret = false;
    ret |= assignFromKeyfile(keyfile, group_name, prefix + "HueMask" + suffix, hueMask);
    ret |= assignFromKeyfile(keyfile, group_name, prefix + "ChromaticityMask" + suffix, chromaticityMask);
    ret |= assignFromKeyfile(keyfile, group_name, prefix + "LightnessMask" + suffix, lightnessMask);
    ret |= assignFromKeyfile(keyfile, group_name, prefix + "MaskBlur" + suffix, maskBlur);
    ret |= assignFromKeyfile(keyfile, group_name, prefix + "AreaMaskEnabled" + suffix, areaEnabled);
    ret |= assignFromKeyfile(keyfile, group_name, prefix + "AreaMaskInverted" + suffix, areaMask.inverted);
    ret |= assignFromKeyfile(keyfile, group_name, prefix + "AreaMaskFeather" + suffix, areaMask.feather);
    ret |= assignFromKeyfile(keyfile, group_name, prefix + "AreaMaskContrast" + suffix, areaMask.contrast);
    if (areaMask.contrast.empty() || areaMask.contrast[0] < DCT_Linear || areaMask.contrast[0] >= DCT_Unchanged) {
        areaMask.contrast = {DCT_Linear};
    }
    std::vector<AreaMask::Shape> s;
    for (int i = 0; ; ++i) {
        AreaMask::Shape a;
        bool found = false;
        std::string n = i ? std::string("_") + std::to_string(i) + "_" : "";
        found |= assignFromKeyfile(keyfile, group_name, prefix + "AreaMask" + n + "X" + suffix, a.x);
        found |= assignFromKeyfile(keyfile, group_name, prefix + "AreaMask" + n + "Y" + suffix, a.y);
        found |= assignFromKeyfile(keyfile, group_name, prefix + "AreaMask" + n + "Width" + suffix, a.width);
        found |= assignFromKeyfile(keyfile, group_name, prefix + "AreaMask" + n + "Height" + suffix, a.height);
        found |= assignFromKeyfile(keyfile, group_name, prefix + "AreaMask" + n + "Angle" + suffix, a.angle);
        found |= assignFromKeyfile(keyfile, group_name, prefix + "AreaMask" + n + "Roundness" + suffix, a.roundness);
        if (found) {
            s.emplace_back(a);
            ret = true;
        } else {
            break;
        }
    }
    if (!s.empty()) {
        areaMask.shapes = std::move(s);
    }
    return ret;
}


void LabCorrectionMask::save(Glib::KeyFile &keyfile, const Glib::ustring &group_name, const Glib::ustring &prefix, const Glib::ustring &suffix) const
{
    putToKeyfile(group_name, prefix + "HueMask" + suffix, hueMask, keyfile);
    putToKeyfile(group_name, prefix + "ChromaticityMask" + suffix, chromaticityMask, keyfile);
    putToKeyfile(group_name, prefix + "LightnessMask" + suffix, lightnessMask, keyfile);
    putToKeyfile(group_name, prefix + "MaskBlur" + suffix, maskBlur, keyfile);
    putToKeyfile(group_name, prefix + "AreaMaskEnabled" + suffix, areaEnabled, keyfile);
    putToKeyfile(group_name, prefix + "AreaMaskInverted" + suffix, areaMask.inverted, keyfile);
    putToKeyfile(group_name, prefix + "AreaMaskFeather" + suffix, areaMask.feather, keyfile);
    putToKeyfile(group_name, prefix + "AreaMaskContrast" + suffix, areaMask.contrast, keyfile);
    for (size_t i = 0; i < areaMask.shapes.size(); ++i) {
        auto &a = areaMask.shapes[i];
        std::string n = i ? std::string("_") + std::to_string(i) + "_" : "";
        putToKeyfile(group_name, prefix + "AreaMask" + n + "X" + suffix, a.x, keyfile);
        putToKeyfile(group_name, prefix + "AreaMask" + n + "Y" + suffix, a.y, keyfile);
        putToKeyfile(group_name, prefix + "AreaMask" + n + "Width" + suffix, a.width, keyfile);
        putToKeyfile(group_name, prefix + "AreaMask" + n + "Height" + suffix, a.height, keyfile);
        putToKeyfile(group_name, prefix + "AreaMask" + n + "Angle" + suffix, a.angle, keyfile);
        putToKeyfile(group_name, prefix + "AreaMask" + n + "Roundness" + suffix, a.roundness, keyfile);
    }
}


ExposureParams::ExposureParams():
    enabled(true),
    autoexp(false),
    clip(0.02),
    hrenabled(false),
    method("Blend"),
    expcomp(0),
    black(0),
    shcompr(50),
    hlcompr(0),
    hlcomprthresh(0),
    clampOOG(true)
{
}


bool ExposureParams::operator==(const ExposureParams &other) const
{
    return enabled == other.enabled
        && autoexp == other.autoexp
        && clip == other.clip
        && hrenabled == other.hrenabled
        && method == other.method
        && expcomp == other.expcomp
        && black == other.black
        && shcompr == other.shcompr
        && hlcompr == other.hlcompr
        && hlcomprthresh == other.hlcomprthresh
        && clampOOG == other.clampOOG;    
}


bool ExposureParams::operator!=(const ExposureParams &other) const
{
    return !(*this == other);
}


BrightnessContrastSaturationParams::BrightnessContrastSaturationParams():
    enabled(false),
    brightness(0),
    contrast(0),
    saturation(0)
{
}


bool BrightnessContrastSaturationParams::operator==(const BrightnessContrastSaturationParams &other) const
{
    return enabled == other.enabled
        && brightness == other.brightness
        && contrast == other.contrast
        && saturation == other.saturation;
}


bool BrightnessContrastSaturationParams::operator!=(const BrightnessContrastSaturationParams &other) const
{
    return !(*this == other);
}


ToneCurveParams::ToneCurveParams():
    enabled(false),
    curve{
        DCT_Linear
    },
    curve2{
        DCT_Linear
    },
    curveMode(ToneCurveParams::TcMode::STD),
    curveMode2(ToneCurveParams::TcMode::STD),
    histmatching(false),
    fromHistMatching(false)
{
}


bool ToneCurveParams::operator ==(const ToneCurveParams& other) const
{
    return enabled == other.enabled
        && curve == other.curve
        && curve2 == other.curve2
        && curveMode == other.curveMode
        && curveMode2 == other.curveMode2
        && histmatching == other.histmatching
        && fromHistMatching == other.fromHistMatching;
}


bool ToneCurveParams::operator !=(const ToneCurveParams& other) const
{
    return !(*this == other);
}


LCurveParams::LCurveParams() :
    enabled(false),
    lcurve{
        DCT_Linear
    },
    acurve{
        DCT_Linear
    },
    bcurve{
        DCT_Linear
    },
    cccurve{
        DCT_Linear
    },
    chcurve{
        FCT_Linear
    },
    lhcurve{
        FCT_Linear
    },
    hhcurve{
        FCT_Linear
    },
    lccurve{
        DCT_Linear
    },
    clcurve{
        DCT_Linear
    },
    brightness(0),
    contrast(0),
    chromaticity(0),
    avoidcolorshift(false),
    rstprotection(0),
    lcredsk(true)
{
}

bool LCurveParams::operator ==(const LCurveParams& other) const
{
    return
        enabled == other.enabled
        && lcurve == other.lcurve
        && acurve == other.acurve
        && bcurve == other.bcurve
        && cccurve == other.cccurve
        && chcurve == other.chcurve
        && lhcurve == other.lhcurve
        && hhcurve == other.hhcurve
        && lccurve == other.lccurve
        && clcurve == other.clcurve
        && brightness == other.brightness
        && contrast == other.contrast
        && chromaticity == other.chromaticity
        && avoidcolorshift == other.avoidcolorshift
        && rstprotection == other.rstprotection
        && lcredsk == other.lcredsk;
}

bool LCurveParams::operator !=(const LCurveParams& other) const
{
    return !(*this == other);
}

RGBCurvesParams::RGBCurvesParams() :
    enabled(false),
    lumamode(false),
    rcurve{
        DCT_Linear
    },
    gcurve{
        DCT_Linear
    },
    bcurve{
        DCT_Linear
    }
{
}

bool RGBCurvesParams::operator ==(const RGBCurvesParams& other) const
{
    return
        enabled == other.enabled
        && lumamode == other.lumamode
        && rcurve == other.rcurve
        && gcurve == other.gcurve
        && bcurve == other.bcurve;
}

bool RGBCurvesParams::operator !=(const RGBCurvesParams& other) const
{
    return !(*this == other);
}


LocalContrastParams::LocalContrastParams():
    enabled(false),
    mode(USM),
    radius(80),
    amount(0.2),
    darkness(1.0),
    lightness(1.0),
    contrast(0),
    curve{
        static_cast<double>(FCT_MinMaxCPoints),
        0.0,
        0.50,
        0.35,
        0.35,
        1.00,
        0.50,
        0.35,
        0.35
    }
{
}


bool LocalContrastParams::operator==(const LocalContrastParams &other) const
{
    return
        enabled == other.enabled
        && mode == other.mode
        && radius == other.radius
        && amount == other.amount
        && darkness == other.darkness
        && lightness == other.lightness
        && contrast == other.contrast
        && curve == other.curve;
}


bool LocalContrastParams::operator!=(const LocalContrastParams &other) const
{
    return !(*this == other);
}


SharpeningParams::SharpeningParams() :
    enabled(false),
    contrast(20.0),
    blurradius(0.2),
    radius(0.5),
    amount(200),
    threshold(20, 80, 2000, 1200, false),
    edgesonly(false),
    edges_radius(1.9),
    edges_tolerance(1800),
    halocontrol(false),
    halocontrol_amount(85),
    method("usm"),
    deconvamount(100),
    deconvradius(0.75),
    deconviter(30),
    deconvdamping(0)
{
}

bool SharpeningParams::operator ==(const SharpeningParams& other) const
{
    return
        enabled == other.enabled
        && contrast == other.contrast
        && blurradius == other.blurradius
        && radius == other.radius
        && amount == other.amount
        && threshold == other.threshold
        && edgesonly == other.edgesonly
        && edges_radius == other.edges_radius
        && edges_tolerance == other.edges_tolerance
        && halocontrol == other.halocontrol
        && halocontrol_amount == other.halocontrol_amount
        && method == other.method
        && deconvamount == other.deconvamount
        && deconvradius == other.deconvradius
        && deconviter == other.deconviter
        && deconvdamping == other.deconvdamping;
}

bool SharpeningParams::operator !=(const SharpeningParams& other) const
{
    return !(*this == other);
}


SharpenMicroParams::SharpenMicroParams() :
    enabled(false),
    matrix(false),
    amount(20.0),
    contrast(20.0),
    uniformity(5)
{
}

bool SharpenMicroParams::operator ==(const SharpenMicroParams& other) const
{
    return
        enabled == other.enabled
        && matrix == other.matrix
        && amount == other.amount
        && contrast == other.contrast
        && uniformity == other.uniformity;
}

bool SharpenMicroParams::operator !=(const SharpenMicroParams& other) const
{
    return !(*this == other);
}

WBParams::WBParams() :
    enabled(true),
    method("Camera"),
    temperature(6504),
    green(1.0),
    equal(1.0),
    tempBias(0.0)
{
}

bool WBParams::operator ==(const WBParams& other) const
{
    return
        enabled == other.enabled
        && method == other.method
        && temperature == other.temperature
        && green == other.green
        && equal == other.equal
        && tempBias == other.tempBias;
}

bool WBParams::operator !=(const WBParams& other) const
{
    return !(*this == other);
}

const std::vector<WBEntry>& WBParams::getWbEntries()
{
    static const std::vector<WBEntry> wb_entries = {
        {"Camera",               WBEntry::Type::CAMERA,      M("TP_WBALANCE_CAMERA"),         0, 1.f,   1.f,   0.f},
        {"Auto",                 WBEntry::Type::AUTO,        M("TP_WBALANCE_AUTO"),           0, 1.f,   1.f,   0.f},
        {"Daylight",             WBEntry::Type::DAYLIGHT,    M("TP_WBALANCE_DAYLIGHT"),    5300, 1.f,   1.f,   0.f},
        {"Cloudy",               WBEntry::Type::CLOUDY,      M("TP_WBALANCE_CLOUDY"),      6200, 1.f,   1.f,   0.f},
        {"Shade",                WBEntry::Type::SHADE,       M("TP_WBALANCE_SHADE"),       7600, 1.f,   1.f,   0.f},
        {"Water 1",              WBEntry::Type::WATER,       M("TP_WBALANCE_WATER1"),     35000, 0.3f,  1.1f,  0.f},
        {"Water 2",              WBEntry::Type::WATER,       M("TP_WBALANCE_WATER2"),     48000, 0.63f, 1.38f, 0.f},
        {"Tungsten",             WBEntry::Type::TUNGSTEN,    M("TP_WBALANCE_TUNGSTEN"),    2856, 1.f,   1.f,   0.f},
        {"Fluo F1",              WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO1"),       6430, 1.f,   1.f,   0.f},
        {"Fluo F2",              WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO2"),       4230, 1.f,   1.f,   0.f},
        {"Fluo F3",              WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO3"),       3450, 1.f,   1.f,   0.f},
        {"Fluo F4",              WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO4"),       2940, 1.f,   1.f,   0.f},
        {"Fluo F5",              WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO5"),       6350, 1.f,   1.f,   0.f},
        {"Fluo F6",              WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO6"),       4150, 1.f,   1.f,   0.f},
        {"Fluo F7",              WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO7"),       6500, 1.f,   1.f,   0.f},
        {"Fluo F8",              WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO8"),       5020, 1.f,   1.f,   0.f},
        {"Fluo F9",              WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO9"),       4330, 1.f,   1.f,   0.f},
        {"Fluo F10",             WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO10"),      5300, 1.f,   1.f,   0.f},
        {"Fluo F11",             WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO11"),      4000, 1.f,   1.f,   0.f},
        {"Fluo F12",             WBEntry::Type::FLUORESCENT, M("TP_WBALANCE_FLUO12"),      3000, 1.f,   1.f,   0.f},
        {"HMI Lamp",             WBEntry::Type::LAMP,        M("TP_WBALANCE_HMI"),         4800, 1.f,   1.f,   0.f},
        {"GTI Lamp",             WBEntry::Type::LAMP,        M("TP_WBALANCE_GTI"),         5000, 1.f,   1.f,   0.f},
        {"JudgeIII Lamp",        WBEntry::Type::LAMP,        M("TP_WBALANCE_JUDGEIII"),    5100, 1.f,   1.f,   0.f},
        {"Solux Lamp 3500K",     WBEntry::Type::LAMP,        M("TP_WBALANCE_SOLUX35"),     3480, 1.f,   1.f,   0.f},
        {"Solux Lamp 4100K",     WBEntry::Type::LAMP,        M("TP_WBALANCE_SOLUX41"),     3930, 1.f,   1.f,   0.f},
        {"Solux Lamp 4700K",     WBEntry::Type::LAMP,        M("TP_WBALANCE_SOLUX47"),     4700, 1.f,   1.f,   0.f},
        {"NG Solux Lamp 4700K",  WBEntry::Type::LAMP,        M("TP_WBALANCE_SOLUX47_NG"),  4480, 1.f,   1.f,   0.f},
        {"LED LSI Lumelex 2040", WBEntry::Type::LED,         M("TP_WBALANCE_LED_LSI"),     2970, 1.f,   1.f,   0.f},
        {"LED CRS SP12 WWMR16",  WBEntry::Type::LED,         M("TP_WBALANCE_LED_CRS"),     3050, 1.f,   1.f,   0.f},
        {"Flash 5500K",          WBEntry::Type::FLASH,       M("TP_WBALANCE_FLASH55"),     5500, 1.f,   1.f,   0.f},
        {"Flash 6000K",          WBEntry::Type::FLASH,       M("TP_WBALANCE_FLASH60"),     6000, 1.f,   1.f,   0.f},
        {"Flash 6500K",          WBEntry::Type::FLASH,       M("TP_WBALANCE_FLASH65"),     6500, 1.f,   1.f,   0.f},
        // Should remain the last one
        {"Custom",               WBEntry::Type::CUSTOM,      M("TP_WBALANCE_CUSTOM"),        0, 1.f,   1.f,   0.f}
    };

    return wb_entries;
}


DefringeParams::DefringeParams() :
    enabled(false),
    radius(2.0),
    threshold(13),
    huecurve{
        FCT_MinMaxCPoints,
        0.166666667,
        0.,
        0.35,
        0.35,
        0.347,
        0.,
        0.35,
        0.35,
        0.513667426,
        0,
        0.35,
        0.35,
        0.668944571,
        0.,
        0.35,
        0.35,
        0.8287775246,
        0.97835991,
        0.35,
        0.35,
        0.9908883827,
        0.,
        0.35,
        0.35
    }
{
}

bool DefringeParams::operator ==(const DefringeParams& other) const
{
    return
        enabled == other.enabled
        && radius == other.radius
        && threshold == other.threshold
        && huecurve == other.huecurve;
}

bool DefringeParams::operator !=(const DefringeParams& other) const
{
    return !(*this == other);
}

ImpulseDenoiseParams::ImpulseDenoiseParams() :
    enabled(false),
    thresh(50)
{
}

bool ImpulseDenoiseParams::operator ==(const ImpulseDenoiseParams& other) const
{
    return
        enabled == other.enabled
        && thresh == other.thresh;
}

bool ImpulseDenoiseParams::operator !=(const ImpulseDenoiseParams& other) const
{
    return !(*this == other);
}

DenoiseParams::DenoiseParams() :
    enabled(false),
    colorSpace(ColorSpace::LAB),
    aggressive(false),
    gamma(1.7),
    luminance(0),
    luminanceDetail(0),
    chrominanceMethod(ChrominanceMethod::AUTOMATIC),
    chrominanceAutoFactor(1),
    chrominance(15),
    chrominanceRedGreen(0),
    chrominanceBlueYellow(0),
    smoothingEnabled(false),
    smoothingMethod(SmoothingMethod::MEDIAN),
    medianType(MedianType::TYPE_3X3_SOFT),
    medianMethod(MedianMethod::CHROMINANCE),
    medianIterations(1),
    guidedLumaRadius(2),
    guidedChromaRadius(4),
    guidedLumaStrength(0),
    guidedChromaStrength(100)
{
}


bool DenoiseParams::operator ==(const DenoiseParams& other) const
{
    return
        enabled == other.enabled
        && colorSpace == other.colorSpace
        && aggressive == other.aggressive
        && gamma == other.gamma
        && luminance == other.luminance
        && luminanceDetail == other.luminanceDetail
        && chrominanceMethod == other.chrominanceMethod
        && chrominanceAutoFactor == other.chrominanceAutoFactor
        && chrominance == other.chrominance
        && chrominanceRedGreen == other.chrominanceRedGreen
        && chrominanceBlueYellow == other.chrominanceBlueYellow
        && smoothingEnabled == other.smoothingEnabled
        && smoothingMethod == other.smoothingMethod
        && medianType == other.medianType
        && medianMethod == other.medianMethod
        && medianIterations == other.medianIterations
        && guidedLumaRadius == other.guidedLumaRadius
        && guidedChromaRadius == other.guidedChromaRadius
        && guidedLumaStrength == other.guidedLumaStrength
        && guidedChromaStrength == other.guidedChromaStrength;
}


bool DenoiseParams::operator !=(const DenoiseParams& other) const
{
    return !(*this == other);
}


EPDParams::Region::Region():
    strength(0.5),
    gamma(1.0),
    edgeStopping(1.4),
    scale(1.0),
    reweightingIterates(0)
{
}


bool EPDParams::Region::operator==(const Region &other) const
{
    return strength == other.strength
        && gamma == other.gamma
        && edgeStopping == other.edgeStopping
        && scale == other.scale
        && reweightingIterates == other.reweightingIterates;
}


bool EPDParams::Region::operator!=(const Region &other) const
{
    return !(*this == other);
}


EPDParams::EPDParams() :
    enabled(false),
    regions{Region()},
    labmasks{LabCorrectionMask()},
    showMask(-1)
{
}

bool EPDParams::operator ==(const EPDParams& other) const
{
    return
        enabled == other.enabled
        && regions == other.regions
        && labmasks == other.labmasks
        && showMask == other.showMask;
}

bool EPDParams::operator !=(const EPDParams& other) const
{
    return !(*this == other);
}


LogEncodingParams::LogEncodingParams():
    enabled(false),
    autocompute(true),
    autogray(true),
    sourceGray(18.0),
    targetGray(18.0),
    blackEv(-5.0),
    whiteEv(10.0),
    detail(1)
{
}

bool LogEncodingParams::operator ==(const LogEncodingParams& other) const
{
    return
        enabled == other.enabled
        && autocompute == other.autocompute
        && autogray == other.autogray
        && sourceGray == other.sourceGray
        && blackEv == other.blackEv
        && whiteEv == other.whiteEv
        && targetGray == other.targetGray
        && detail == other.detail;
}

bool LogEncodingParams::operator !=(const LogEncodingParams& other) const
{
    return !(*this == other);
}


FattalToneMappingParams::FattalToneMappingParams() :
    enabled(false),
    threshold(30),
    amount(20),
    anchor(50)
{
}

bool FattalToneMappingParams::operator ==(const FattalToneMappingParams& other) const
{
    return
        enabled == other.enabled
        && threshold == other.threshold
        && amount == other.amount
        && anchor == other.anchor;
}

bool FattalToneMappingParams::operator !=(const FattalToneMappingParams& other) const
{
    return !(*this == other);
}

SHParams::SHParams() :
    enabled(false),
    highlights(0),
    htonalwidth(70),
    shadows(0),
    stonalwidth(30),
    radius(40),
    lab(false)
{
}

bool SHParams::operator ==(const SHParams& other) const
{
    return
        enabled == other.enabled
        && highlights == other.highlights
        && htonalwidth == other.htonalwidth
        && shadows == other.shadows
        && stonalwidth == other.stonalwidth
        && radius == other.radius
        && lab == other.lab;
}

bool SHParams::operator !=(const SHParams& other) const
{
    return !(*this == other);
}


ToneEqualizerParams::ToneEqualizerParams():
    enabled(false),
    bands{0,0,0,0,0},
    detail(0)
{
}


bool ToneEqualizerParams::operator ==(const ToneEqualizerParams& other) const
{
    return
        enabled == other.enabled
        && bands == other.bands
        && detail == other.detail;
}


bool ToneEqualizerParams::operator !=(const ToneEqualizerParams& other) const
{
    return !(*this == other);
}


CropParams::CropParams() :
    enabled(false),
    x(-1),
    y(-1),
    w(15000),
    h(15000),
    fixratio(true),
    ratio("As Image"),
    orientation("As Image"),
    guide("Frame")
{
}

bool CropParams::operator ==(const CropParams& other) const
{
    return
        enabled == other.enabled
        && x == other.x
        && y == other.y
        && w == other.w
        && h == other.h
        && fixratio == other.fixratio
        && ratio == other.ratio
        && orientation == other.orientation
        && guide == other.guide;
}

bool CropParams::operator !=(const CropParams& other) const
{
    return !(*this == other);
}

void CropParams::mapToResized(int resizedWidth, int resizedHeight, int scale, int& x1, int& x2, int& y1, int& y2) const
{
    x1 = 0, x2 = resizedWidth, y1 = 0, y2 = resizedHeight;

    if (enabled) {
        x1 = min(resizedWidth - 1, max(0, x / scale));
        y1 = min(resizedHeight - 1, max(0, y / scale));
        x2 = min(resizedWidth, max(0, (x + w) / scale));
        y2 = min(resizedHeight, max(0, (y + h) / scale));
    }
}

CoarseTransformParams::CoarseTransformParams() :
    rotate(0),
    hflip(false),
    vflip(false)
{
}

bool CoarseTransformParams::operator ==(const CoarseTransformParams& other) const
{
    return
        rotate == other.rotate
        && hflip == other.hflip
        && vflip == other.vflip;
}

bool CoarseTransformParams::operator !=(const CoarseTransformParams& other) const
{
    return !(*this == other);
}

CommonTransformParams::CommonTransformParams() :
    autofill(true)
{
}

bool CommonTransformParams::operator ==(const CommonTransformParams& other) const
{
    return autofill == other.autofill;
}

bool CommonTransformParams::operator !=(const CommonTransformParams& other) const
{
    return !(*this == other);
}

RotateParams::RotateParams() :
    enabled(false),
    degree(0.0)
{
}

bool RotateParams::operator ==(const RotateParams& other) const
{
    return enabled == other.enabled && degree == other.degree;
}

bool RotateParams::operator !=(const RotateParams& other) const
{
    return !(*this == other);
}

DistortionParams::DistortionParams() :
    enabled(false),
    amount(0.0)
{
}

bool DistortionParams::operator ==(const DistortionParams& other) const
{
    return enabled == other.enabled && amount == other.amount;
}

bool DistortionParams::operator !=(const DistortionParams& other) const
{
    return !(*this == other);
}

LensProfParams::LensProfParams() :
    lcMode(LcMode::NONE),
    useDist(true),
    useVign(true),
    useCA(false)
{
}

bool LensProfParams::operator ==(const LensProfParams& other) const
{
    return
        lcMode == other.lcMode
        && lcpFile == other.lcpFile
        && useCA == other.useCA
        && lfCameraMake == other.lfCameraMake
        && lfCameraModel == other.lfCameraModel
        && lfLens == other.lfLens;
}

bool LensProfParams::operator !=(const LensProfParams& other) const
{
    return !(*this == other);
}

bool LensProfParams::useLensfun() const
{
    return lcMode == LcMode::LENSFUNAUTOMATCH || lcMode == LcMode::LENSFUNMANUAL;
}

bool LensProfParams::lfAutoMatch() const
{
    return lcMode == LcMode::LENSFUNAUTOMATCH;
}

bool LensProfParams::useLcp() const
{
    return lcMode == LcMode::LCP && lcpFile.length() > 0;
}

bool LensProfParams::lfManual() const
{
    return lcMode == LcMode::LENSFUNMANUAL;
}

const std::vector<const char*>& LensProfParams::getMethodStrings() const
{
    static const std::vector<const char*> method_strings = {
        "none",
        "lfauto",
        "lfmanual",
        "lcp"
    };
    return method_strings;
}

Glib::ustring LensProfParams::getMethodString(LcMode mode) const
{
    return getMethodStrings()[toUnderlying(mode)];
}

LensProfParams::LcMode LensProfParams::getMethodNumber(const Glib::ustring& mode) const
{
    for (std::vector<const char*>::size_type i = 0; i < getMethodStrings().size(); ++i) {
        if (getMethodStrings()[i] == mode) {
            return static_cast<LcMode>(i);
        }
    }

    return LcMode::NONE;
}

PerspectiveParams::PerspectiveParams() :
    enabled(false),
    horizontal(0.0),
    vertical(0.0),
    angle(0.0),
    shear(0.0),
    flength(0),
    cropfactor(1),
    aspect(1)
{
}

bool PerspectiveParams::operator ==(const PerspectiveParams& other) const
{
    return
        enabled == other.enabled
        && horizontal == other.horizontal
        && vertical == other.vertical
        && angle == other.angle
        && shear == other.shear
        && flength == other.flength
        && cropfactor == other.cropfactor
        && aspect == other.aspect;
}

bool PerspectiveParams::operator !=(const PerspectiveParams& other) const
{
    return !(*this == other);
}

GradientParams::GradientParams() :
    enabled(false),
    degree(0.0),
    feather(25),
    strength(0.60),
    centerX(0),
    centerY(0)
{
}

bool GradientParams::operator ==(const GradientParams& other) const
{
    return
        enabled == other.enabled
        && degree == other.degree
        && feather == other.feather
        && strength == other.strength
        && centerX == other.centerX
        && centerY == other.centerY;
}

bool GradientParams::operator !=(const GradientParams& other) const
{
    return !(*this == other);
}

PCVignetteParams::PCVignetteParams() :
    enabled(false),
    strength(0.60),
    feather(50),
    roundness(50)
{
}

bool PCVignetteParams::operator ==(const PCVignetteParams& other) const
{
    return
        enabled == other.enabled
        && strength == other.strength
        && feather == other.feather
        && roundness == other.roundness;
}

bool PCVignetteParams::operator !=(const PCVignetteParams& other) const
{
    return !(*this == other);
}

VignettingParams::VignettingParams() :
    enabled(false),
    amount(0),
    radius(50),
    strength(1),
    centerX(0),
    centerY(0)
{
}

bool VignettingParams::operator ==(const VignettingParams& other) const
{
    return
        enabled == other.enabled
        && amount == other.amount
        && radius == other.radius
        && strength == other.strength
        && centerX == other.centerX
        && centerY == other.centerY;
}

bool VignettingParams::operator !=(const VignettingParams& other) const
{
    return !(*this == other);
}

ChannelMixerParams::ChannelMixerParams() :
    enabled(false),
    red{
        1000,
        0,
        0
    },
    green{
        0,
        1000,
        0
    },
    blue{
        0,
        0,
        1000
    }
{
}

bool ChannelMixerParams::operator ==(const ChannelMixerParams& other) const
{
    if (enabled != other.enabled) {
        return false;
    }

    for (unsigned int i = 0; i < 3; ++i) {
        if (
            red[i] != other.red[i]
            || green[i] != other.green[i]
            || blue[i] != other.blue[i]
        ) {
            return false;
        }
    }

    return true;
}

bool ChannelMixerParams::operator !=(const ChannelMixerParams& other) const
{
    return !(*this == other);
}

BlackWhiteParams::BlackWhiteParams() :
    enabled(false),
    method("Desaturation"),
    luminanceCurve{
        FCT_Linear
    },
    filter("None"),
    setting("RGB-Rel"),
    mixerRed(33),
    mixerGreen(33),
    mixerBlue(33),
    gammaRed(0),
    gammaGreen(0),
    gammaBlue(0)
{
}

bool BlackWhiteParams::operator ==(const BlackWhiteParams& other) const
{
    return
        enabled == other.enabled
        && method == other.method
        && luminanceCurve == other.luminanceCurve
        && filter == other.filter
        && setting == other.setting
        && mixerRed == other.mixerRed
        && mixerGreen == other.mixerGreen
        && mixerBlue == other.mixerBlue
        && gammaRed == other.gammaRed
        && gammaGreen == other.gammaGreen
        && gammaBlue == other.gammaBlue;
}

bool BlackWhiteParams::operator !=(const BlackWhiteParams& other) const
{
    return !(*this == other);
}

CACorrParams::CACorrParams() :
    enabled(false),
    red(0.0),
    blue(0.0)
{
}

bool CACorrParams::operator ==(const CACorrParams& other) const
{
    return
        enabled == other.enabled
        && red == other.red
        && blue == other.blue;
}

bool CACorrParams::operator !=(const CACorrParams& other) const
{
    return !(*this == other);
}

ResizeParams::ResizeParams() :
    enabled(false),
    scale(1.0),
    appliesTo("Cropped area"),
    method("Lanczos"),
    dataspec(3),
    width(900),
    height(900),
    allowUpscaling(false)
{
}

bool ResizeParams::operator ==(const ResizeParams& other) const
{
    return
        enabled == other.enabled
        && scale == other.scale
        && appliesTo == other.appliesTo
        && method == other.method
        && dataspec == other.dataspec
        && width == other.width
        && height == other.height
        && allowUpscaling == other.allowUpscaling;
}

bool ResizeParams::operator !=(const ResizeParams& other) const
{
    return !(*this == other);
}

const Glib::ustring ColorManagementParams::NoICMString = Glib::ustring("No ICM: sRGB output");

ColorManagementParams::ColorManagementParams() :
    inputProfile("(cameraICC)"),
    toneCurve(false),
    applyLookTable(false),
    applyBaselineExposureOffset(true),
    applyHueSatMap(true),
    dcpIlluminant(0),
    workingProfile("ProPhoto"),
    outputProfile(options.rtSettings.srgb),
    outputIntent(RI_RELATIVE),
    outputBPC(true)
{
}

bool ColorManagementParams::operator ==(const ColorManagementParams& other) const
{
    return
        inputProfile == other.inputProfile
        && toneCurve == other.toneCurve
        && applyLookTable == other.applyLookTable
        && applyBaselineExposureOffset == other.applyBaselineExposureOffset
        && applyHueSatMap == other.applyHueSatMap
        && dcpIlluminant == other.dcpIlluminant
        && workingProfile == other.workingProfile
        && outputProfile == other.outputProfile
        && outputIntent == other.outputIntent
        && outputBPC == other.outputBPC;
}

bool ColorManagementParams::operator !=(const ColorManagementParams& other) const
{
    return !(*this == other);
}


DirPyrEqualizerParams::Levels::Levels():
    mult{
        1.0,
        1.0,
        1.0,
        1.0,
        1.0,
        1.0
    },
    threshold(0.2)
{
}


bool DirPyrEqualizerParams::Levels::operator==(const Levels &other) const
{
    for (int i = 0; i < 6; ++i) {
        if (mult[i] != other.mult[i]) {
            return false;
        }
    }
    return threshold == other.threshold;
}


bool DirPyrEqualizerParams::Levels::operator!=(const Levels &other) const
{
    return !(*this == other);
}


DirPyrEqualizerParams::DirPyrEqualizerParams() :
    enabled(false),
    levels{Levels()},
    labmasks{LabCorrectionMask()},
    showMask(-1)
    // threshold(0.2),
    // skinprotect(0.0),
    // hueskin (-5, 25, 170, 120, false),
    // cbdlMethod("bef")
{
}


bool DirPyrEqualizerParams::operator ==(const DirPyrEqualizerParams& other) const
{
    return
        enabled == other.enabled
        && levels == other.levels
        && labmasks == other.labmasks
        && showMask == other.showMask;
        // && gamutlab == other.gamutlab
        // && [this, &other]() -> bool
        //     {
        //         for (unsigned int i = 0; i < 6; ++i) {
        //             if (mult[i] != other.mult[i]) {
        //                 return false;
        //             }
        //         }
        //         return true;
        //     }()
        // && threshold == other.threshold
        // && skinprotect == other.skinprotect
        // && hueskin == other.hueskin
        // && cbdlMethod == other.cbdlMethod;
}

bool DirPyrEqualizerParams::operator !=(const DirPyrEqualizerParams& other) const
{
    return !(*this == other);
}


FilmSimulationParams::FilmSimulationParams() :
    enabled(false),
    strength(100)
{
}

bool FilmSimulationParams::operator ==(const FilmSimulationParams& other) const
{
    return
        enabled == other.enabled
        && clutFilename == other.clutFilename
        && strength == other.strength;
}

bool FilmSimulationParams::operator !=(const FilmSimulationParams& other) const
{
    return !(*this == other);
}


SoftLightParams::SoftLightParams() :
    enabled(false),
    strength(30)
{
}

bool SoftLightParams::operator ==(const SoftLightParams& other) const
{
    return
        enabled == other.enabled
        && strength == other.strength;
}

bool SoftLightParams::operator !=(const SoftLightParams& other) const
{
    return !(*this == other);
}


DehazeParams::DehazeParams() :
    enabled(false),
    strength(50),
    showDepthMap(false),
    depth(25),
    luminance(false)
{
}

bool DehazeParams::operator ==(const DehazeParams& other) const
{
    return
        enabled == other.enabled
        && strength == other.strength
        && showDepthMap == other.showDepthMap
        && depth == other.depth
        && luminance == other.luminance;
}

bool DehazeParams::operator !=(const DehazeParams& other) const
{
    return !(*this == other);
}


GrainParams::GrainParams():
    enabled(false),
    iso(400),
    strength(25),
    scale(100)
{
}

bool GrainParams::operator==(const GrainParams &other) const
{
    return enabled == other.enabled
        && iso == other.iso
        && strength == other.strength
        && scale == other.scale;
}

bool GrainParams::operator!=(const GrainParams &other) const
{
    return !(*this == other);
}


GuidedSmoothingParams::Region::Region():
    channel(Channel::RGB),
    radius(0),
    epsilon(0)
{
}


bool GuidedSmoothingParams::Region::operator==(const Region &other) const
{
    return channel == other.channel
        && radius == other.radius
        && epsilon == other.epsilon;
}


bool GuidedSmoothingParams::Region::operator!=(const Region &other) const
{
    return !(*this == other);
}


GuidedSmoothingParams::GuidedSmoothingParams():
    enabled(false),
    regions{Region()},
    labmasks{LabCorrectionMask()},
    showMask(-1)
{
}


bool GuidedSmoothingParams::operator==(const GuidedSmoothingParams &other) const
{
    return enabled == other.enabled
        && regions == other.regions
        && labmasks == other.labmasks
        && showMask == other.showMask;
}


bool GuidedSmoothingParams::operator!=(const GuidedSmoothingParams &other) const
{
    return !(*this == other);
}


ColorCorrectionParams::LabCorrectionRegion::LabCorrectionRegion():
    a(0),
    b(0),
    saturation(0),
    slope(1),
    offset(0),
    power(1),
    channel(ColorCorrectionParams::LabCorrectionRegion::CHAN_ALL)
{
}


bool ColorCorrectionParams::LabCorrectionRegion::operator==(const LabCorrectionRegion &other) const
{
    return a == other.a
        && b == other.b
        && saturation == other.saturation
        && slope == other.slope
        && offset == other.offset
        && power == other.power
        && channel == other.channel;
}


bool ColorCorrectionParams::LabCorrectionRegion::operator!=(const LabCorrectionRegion &other) const
{
    return !(*this == other);
}


ColorCorrectionParams::ColorCorrectionParams():
    enabled(false),
    regions{LabCorrectionRegion()},
    labmasks{LabCorrectionMask()},
    showMask(-1)
{
}


bool ColorCorrectionParams::operator==(const ColorCorrectionParams &other) const
{
    return enabled == other.enabled
        && regions == other.regions
        && labmasks == other.labmasks
        && showMask == other.showMask;
}


bool ColorCorrectionParams::operator!=(const ColorCorrectionParams &other) const
{
    return !(*this == other);
}


RAWParams::BayerSensor::BayerSensor() :
    method(getMethodString(Method::AMAZE)),
    border(4),
    imageNum(0),
    ccSteps(0),
    black0(0.0),
    black1(0.0),
    black2(0.0),
    black3(0.0),
    twogreen(true),
    linenoise(0),
    linenoiseDirection(LineNoiseDirection::BOTH),
    greenthresh(0),
    dcb_iterations(2),
    lmmse_iterations(2),
    dualDemosaicAutoContrast(true),
    dualDemosaicContrast(20),
    pixelShiftMotionCorrectionMethod(PSMotionCorrectionMethod::AUTO),
    pixelShiftEperIso(0.0),
    pixelShiftSigma(1.0),
    pixelShiftShowMotion(false),
    pixelShiftShowMotionMaskOnly(false),
    pixelShiftHoleFill(true),
    pixelShiftMedian(false),
    pixelShiftGreen(true),
    pixelShiftBlur(true),
    pixelShiftSmoothFactor(0.7),
    pixelShiftEqualBright(false),
    pixelShiftEqualBrightChannel(false),
    pixelShiftNonGreenCross(true),
    pixelShiftDemosaicMethod(getPSDemosaicMethodString(PSDemosaicMethod::AMAZE)),
    dcb_enhance(true),
    pdafLinesFilter(false),
    enable_black(false),
    enable_preproc(false)
{
}

bool RAWParams::BayerSensor::operator ==(const BayerSensor& other) const
{
    return
        method == other.method
        && border == other.border
        && imageNum == other.imageNum
        && ccSteps == other.ccSteps
        && black0 == other.black0
        && black1 == other.black1
        && black2 == other.black2
        && black3 == other.black3
        && twogreen == other.twogreen
        && linenoise == other.linenoise
        && linenoiseDirection == other.linenoiseDirection
        && greenthresh == other.greenthresh
        && dcb_iterations == other.dcb_iterations
        && lmmse_iterations == other.lmmse_iterations
        && dualDemosaicAutoContrast == other.dualDemosaicAutoContrast
        && dualDemosaicContrast == other.dualDemosaicContrast
        && pixelShiftMotionCorrectionMethod == other.pixelShiftMotionCorrectionMethod
        && pixelShiftEperIso == other.pixelShiftEperIso
        && pixelShiftSigma == other.pixelShiftSigma
        && pixelShiftShowMotion == other.pixelShiftShowMotion
        && pixelShiftShowMotionMaskOnly == other.pixelShiftShowMotionMaskOnly
        && pixelShiftHoleFill == other.pixelShiftHoleFill
        && pixelShiftMedian == other.pixelShiftMedian
        && pixelShiftGreen == other.pixelShiftGreen
        && pixelShiftBlur == other.pixelShiftBlur
        && pixelShiftSmoothFactor == other.pixelShiftSmoothFactor
        && pixelShiftEqualBright == other.pixelShiftEqualBright
        && pixelShiftEqualBrightChannel == other.pixelShiftEqualBrightChannel
        && pixelShiftNonGreenCross == other.pixelShiftNonGreenCross
        && pixelShiftDemosaicMethod == other.pixelShiftDemosaicMethod
        && dcb_enhance == other.dcb_enhance
        && pdafLinesFilter == other.pdafLinesFilter
        && enable_black == other.enable_black
        && enable_preproc == other.enable_preproc;
}

bool RAWParams::BayerSensor::operator !=(const BayerSensor& other) const
{
    return !(*this == other);
}

void RAWParams::BayerSensor::setPixelShiftDefaults()
{
    pixelShiftMotionCorrectionMethod = RAWParams::BayerSensor::PSMotionCorrectionMethod::AUTO;
    pixelShiftEperIso = 0.0;
    pixelShiftSigma = 1.0;
    pixelShiftHoleFill = true;
    pixelShiftMedian = false;
    pixelShiftGreen = true;
    pixelShiftBlur = true;
    pixelShiftSmoothFactor = 0.7;
    pixelShiftEqualBright = false;
    pixelShiftEqualBrightChannel = false;
    pixelShiftNonGreenCross = true;
    pixelShiftDemosaicMethod = getPSDemosaicMethodString(PSDemosaicMethod::AMAZE);
}

const std::vector<const char*>& RAWParams::BayerSensor::getMethodStrings()
{
    static const std::vector<const char*> method_strings {
        "amaze",
        "amazevng4",
        "rcd",
        "rcdvng4",
        "dcb",
        "dcbvng4",
        "lmmse",
        "igv",
        "ahd",
        "eahd",
        "hphd",
        "vng4",
        "fast",
        "mono",
        "pixelshift",
        "none"
    };
    return method_strings;
}

Glib::ustring RAWParams::BayerSensor::getMethodString(Method method)
{
    return getMethodStrings()[toUnderlying(method)];
}

const std::vector<const char*>& RAWParams::BayerSensor::getPSDemosaicMethodStrings()
{
    static const std::vector<const char*> method_strings {
        "amaze",
        "amazevng4",
        "lmmse"
    };
    return method_strings;
}

Glib::ustring RAWParams::BayerSensor::getPSDemosaicMethodString(PSDemosaicMethod method)
{
    return getPSDemosaicMethodStrings()[toUnderlying(method)];
}



RAWParams::XTransSensor::XTransSensor() :
    method(getMethodString(Method::THREE_PASS)),
    dualDemosaicAutoContrast(true),
    dualDemosaicContrast(20),
    border(7),
    ccSteps(0),
    blackred(0.0),
    blackgreen(0.0),
    blackblue(0.0),
    enable_black(false)
{
}

bool RAWParams::XTransSensor::operator ==(const XTransSensor& other) const
{
    return
        method == other.method
        && dualDemosaicAutoContrast == other.dualDemosaicAutoContrast
        && dualDemosaicContrast == other.dualDemosaicContrast
        && border == other.border
        && ccSteps == other.ccSteps
        && blackred == other.blackred
        && blackgreen == other.blackgreen
        && blackblue == other.blackblue
        && enable_black == other.enable_black;
}

bool RAWParams::XTransSensor::operator !=(const XTransSensor& other) const
{
    return !(*this == other);
}

const std::vector<const char*>& RAWParams::XTransSensor::getMethodStrings()
{
    static const std::vector<const char*> method_strings {
        "4-pass",
        "3-pass (best)",
        "2-pass",
        "1-pass (medium)",
        "fast",
        "mono",
        "none"
    };
    return method_strings;
}

Glib::ustring RAWParams::XTransSensor::getMethodString(Method method)
{
    return getMethodStrings()[toUnderlying(method)];
}

RAWParams::RAWParams() :
    df_autoselect(false),
    ff_AutoSelect(false),
    ff_BlurRadius(32),
    ff_BlurType(getFlatFieldBlurTypeString(FlatFieldBlurType::AREA)),
    ff_AutoClipControl(false),
    ff_clipControl(0),
    ca_autocorrect(false),
    ca_avoidcolourshift(true),
    caautoiterations(2),
    cared(0.0),
    cablue(0.0),
    expos(1.0),
    hotPixelFilter(false),
    deadPixelFilter(false),
    hotdeadpix_thresh(100),
    enable_darkframe(false),
    enable_flatfield(false),
    enable_ca(false),
    enable_hotdeadpix(false),
    enable_whitepoint(false)
{
}

bool RAWParams::operator ==(const RAWParams& other) const
{
    return
        bayersensor == other.bayersensor
        && xtranssensor == other.xtranssensor
        && dark_frame == other.dark_frame
        && df_autoselect == other.df_autoselect
        && ff_file == other.ff_file
        && ff_AutoSelect == other.ff_AutoSelect
        && ff_BlurRadius == other.ff_BlurRadius
        && ff_BlurType == other.ff_BlurType
        && ff_AutoClipControl == other.ff_AutoClipControl
        && ff_clipControl == other.ff_clipControl
        && ca_autocorrect == other.ca_autocorrect
        && ca_avoidcolourshift == other.ca_avoidcolourshift
        && caautoiterations == other.caautoiterations
        && cared == other.cared
        && cablue == other.cablue
        && expos == other.expos
        && hotPixelFilter == other.hotPixelFilter
        && deadPixelFilter == other.deadPixelFilter
        && hotdeadpix_thresh == other.hotdeadpix_thresh
        && enable_darkframe == other.enable_darkframe
        && enable_flatfield == other.enable_flatfield
        && enable_ca == other.enable_ca
        && enable_hotdeadpix == other.enable_hotdeadpix
        && enable_whitepoint == other.enable_whitepoint;
}

bool RAWParams::operator !=(const RAWParams& other) const
{
    return !(*this == other);
}

const std::vector<const char*>& RAWParams::getFlatFieldBlurTypeStrings()
{
    static const std::vector<const char*> blur_type_strings {
        "Area Flatfield",
        "Vertical Flatfield",
        "Horizontal Flatfield",
        "V+H Flatfield"
    };
    return blur_type_strings;
}

Glib::ustring RAWParams::getFlatFieldBlurTypeString(FlatFieldBlurType type)
{
    return getFlatFieldBlurTypeStrings()[toUnderlying(type)];
}


MetaDataParams::MetaDataParams():
    mode(MetaDataParams::TUNNEL)
{
}

bool MetaDataParams::operator==(const MetaDataParams &other) const
{
    return mode == other.mode;
}

bool MetaDataParams::operator!=(const MetaDataParams &other) const
{
    return !(*this == other);
}


ProcParams::ProcParams()
{
    setDefaults();
}

void ProcParams::setDefaults()
{
    exposure = ExposureParams();

    brightContrSat = BrightnessContrastSaturationParams();
        
    toneCurve = ToneCurveParams();

    labCurve = LCurveParams();

    rgbCurves = RGBCurvesParams();

    localContrast = LocalContrastParams();

    sharpenMicro = SharpenMicroParams();

    sharpening = SharpeningParams();

    prsharpening = SharpeningParams();
    prsharpening.contrast = 15.0;
    prsharpening.method = "rld";
    prsharpening.deconvamount = 100;
    prsharpening.deconvradius = 0.45;
    prsharpening.deconviter = 100;
    prsharpening.deconvdamping = 0;

    wb = WBParams();

    defringe = DefringeParams();

    impulseDenoise = ImpulseDenoiseParams();

    denoise = DenoiseParams();

    epd = EPDParams();

    fattal = FattalToneMappingParams();

    logenc = LogEncodingParams();

    sh = SHParams();

    toneEqualizer = ToneEqualizerParams();

    crop = CropParams();

    coarse = CoarseTransformParams();

    commonTrans = CommonTransformParams();

    rotate = RotateParams();

    distortion = DistortionParams();

    lensProf = LensProfParams();

    perspective = PerspectiveParams();

    gradient = GradientParams();

    pcvignette = PCVignetteParams();

    vignetting = VignettingParams();

    chmixer = ChannelMixerParams();

    blackwhite = BlackWhiteParams();

    cacorrection = CACorrParams();

    resize = ResizeParams();

    icm = ColorManagementParams();

    dirpyrequalizer = DirPyrEqualizerParams();

    filmSimulation = FilmSimulationParams();

    softlight = SoftLightParams();

    dehaze = DehazeParams();

    grain = GrainParams();

    smoothing = GuidedSmoothingParams();

    colorcorrection = ColorCorrectionParams();

    raw = RAWParams();

    metadata = MetaDataParams();
    exif.clear();
    iptc.clear();

    rank = 0;
    colorlabel = 0;
    inTrash = false;

    ppVersion = PPVERSION;
}


int ProcParams::save(const Glib::ustring& fname, const Glib::ustring& fname2, bool fnameAbsolute, const ParamsEdited *pedited)
{
    if (fname.empty() && fname2.empty()) {
        return 0;
    }

    Glib::ustring sPParams;

    try {
        Glib::KeyFile keyFile;
        int ret = save(keyFile, pedited, fname, fnameAbsolute);
        if (ret != 0) {
            return ret;
        }

        sPParams = keyFile.to_data();
    } catch (Glib::KeyFileError&) {}

    if (sPParams.empty()) {
        return 1;
    }

    int error1, error2;
    error1 = write(fname, sPParams);

    if (!fname2.empty()) {

        error2 = write(fname2, sPParams);
        // If at least one file has been saved, it's a success
        return error1 & error2;
    } else {
        return error1;
    }
}


int ProcParams::save(Glib::KeyFile &keyFile, const ParamsEdited *pedited,
                     const Glib::ustring &fname, bool fnameAbsolute) const
{
#define RELEVANT_(n) (!pedited || pedited->n)
    try {
// Version
        keyFile.set_string("Version", "AppVersion", RTVERSION);
        keyFile.set_integer("Version", "Version", PPVERSION);

        if (RELEVANT_(general)) {
            saveToKeyfile("General", "Rank", rank, keyFile);
            saveToKeyfile("General", "ColorLabel", colorlabel, keyFile);
            saveToKeyfile("General", "InTrash", inTrash, keyFile);
        }

// Exposure
        if (RELEVANT_(exposure)) {
            saveToKeyfile("Exposure", "Enabled", exposure.enabled, keyFile);
            saveToKeyfile("Exposure", "Auto", exposure.autoexp, keyFile);
            saveToKeyfile("Exposure", "Clip", exposure.clip, keyFile);
            saveToKeyfile("Exposure", "Compensation", exposure.expcomp, keyFile);
            saveToKeyfile("Exposure", "Black", exposure.black, keyFile);
            saveToKeyfile("Exposure", "HighlightCompr", exposure.hlcompr, keyFile);
            saveToKeyfile("Exposure", "HighlightComprThreshold", exposure.hlcomprthresh, keyFile);
            saveToKeyfile("Exposure", "ShadowCompr", exposure.shcompr, keyFile);
            saveToKeyfile("Exposure", "ClampOOG", exposure.clampOOG, keyFile);

            saveToKeyfile("Exposure", "HLRecoveryEnabled", exposure.hrenabled, keyFile);
            saveToKeyfile("Exposure", "HLRecoveryMethod", exposure.method, keyFile);
        }

// Brightness, Contrast, Saturation
        if (RELEVANT_(brightContrSat)) {
            saveToKeyfile("BrightnessContrastSaturation", "Enabled", brightContrSat.enabled, keyFile);
            saveToKeyfile("BrightnessContrastSaturation", "Brightness", brightContrSat.brightness, keyFile);
            saveToKeyfile("BrightnessContrastSaturation", "Contrast", brightContrSat.contrast, keyFile);
            saveToKeyfile("BrightnessContrastSaturation", "Saturation", brightContrSat.saturation, keyFile);
        }

// Tone curve
        if (RELEVANT_(toneCurve)) {
            saveToKeyfile("ToneCurve", "Enabled", toneCurve.enabled, keyFile);
            saveToKeyfile("ToneCurve", "HistogramMatching", toneCurve.histmatching, keyFile);
            saveToKeyfile("ToneCurve", "CurveFromHistogramMatching", toneCurve.fromHistMatching, keyFile);

            const std::map<ToneCurveParams::TcMode, const char*> tc_mapping = {
                {ToneCurveParams::TcMode::STD, "Standard"},
                {ToneCurveParams::TcMode::FILMLIKE, "FilmLike"},
                {ToneCurveParams::TcMode::SATANDVALBLENDING, "SatAndValueBlending"},
                {ToneCurveParams::TcMode::WEIGHTEDSTD, "WeightedStd"},
                {ToneCurveParams::TcMode::LUMINANCE, "Luminance"},
                {ToneCurveParams::TcMode::PERCEPTUAL, "Perceptual"}
            };

            saveToKeyfile("ToneCurve", "CurveMode", tc_mapping, toneCurve.curveMode, keyFile);
            saveToKeyfile("ToneCurve", "CurveMode2", tc_mapping, toneCurve.curveMode2, keyFile);

            saveToKeyfile("ToneCurve", "Curve", toneCurve.curve, keyFile);
            saveToKeyfile("ToneCurve", "Curve2", toneCurve.curve2, keyFile);
        }

// Local contrast
        if (RELEVANT_(localContrast)) {
            saveToKeyfile("Local Contrast", "Enabled", localContrast.enabled, keyFile);
            saveToKeyfile("Local Contrast", "Mode", int(localContrast.mode), keyFile);
            saveToKeyfile("Local Contrast", "Radius", localContrast.radius, keyFile);
            saveToKeyfile("Local Contrast", "Amount", localContrast.amount, keyFile);
            saveToKeyfile("Local Contrast", "Darkness", localContrast.darkness, keyFile);
            saveToKeyfile("Local Contrast", "Lightness", localContrast.lightness, keyFile);
            saveToKeyfile("Local Contrast", "Contrast", localContrast.contrast, keyFile);
            saveToKeyfile("Local Contrast", "Curve", localContrast.curve, keyFile);
        }


// Channel mixer
        if (RELEVANT_(chmixer)) {
            saveToKeyfile("Channel Mixer", "Enabled", chmixer.enabled, keyFile);
            Glib::ArrayHandle<int> rmix(chmixer.red, 3, Glib::OWNERSHIP_NONE);
            keyFile.set_integer_list("Channel Mixer", "Red", rmix);
            Glib::ArrayHandle<int> gmix(chmixer.green, 3, Glib::OWNERSHIP_NONE);
            keyFile.set_integer_list("Channel Mixer", "Green", gmix);
            Glib::ArrayHandle<int> bmix(chmixer.blue, 3, Glib::OWNERSHIP_NONE);
            keyFile.set_integer_list("Channel Mixer", "Blue", bmix);
        }

// Black & White
        if (RELEVANT_(blackwhite)) {
            saveToKeyfile("Black & White", "Enabled", blackwhite.enabled, keyFile);
            saveToKeyfile("Black & White", "Method", blackwhite.method, keyFile);
            saveToKeyfile("Black & White", "Setting", blackwhite.setting, keyFile);
            saveToKeyfile("Black & White", "Filter", blackwhite.filter, keyFile);
            saveToKeyfile("Black & White", "MixerRed", blackwhite.mixerRed, keyFile);
            saveToKeyfile("Black & White", "MixerGreen", blackwhite.mixerGreen, keyFile);
            saveToKeyfile("Black & White", "MixerBlue", blackwhite.mixerBlue, keyFile);
            saveToKeyfile("Black & White", "GammaRed", blackwhite.gammaRed, keyFile);
            saveToKeyfile("Black & White", "GammaGreen", blackwhite.gammaGreen, keyFile);
            saveToKeyfile("Black & White", "GammaBlue", blackwhite.gammaBlue, keyFile);
            saveToKeyfile("Black & White", "LuminanceCurve", blackwhite.luminanceCurve, keyFile);
        }

// Luma curve
        if (RELEVANT_(labCurve)) {
            saveToKeyfile("Luminance Curve", "Enabled", labCurve.enabled, keyFile);
            saveToKeyfile("Luminance Curve", "Brightness", labCurve.brightness, keyFile);
            saveToKeyfile("Luminance Curve", "Contrast", labCurve.contrast, keyFile);
            saveToKeyfile("Luminance Curve", "Chromaticity", labCurve.chromaticity, keyFile);
            saveToKeyfile("Luminance Curve", "AvoidColorShift", labCurve.avoidcolorshift, keyFile);
            saveToKeyfile("Luminance Curve", "RedAndSkinTonesProtection", labCurve.rstprotection, keyFile);
            saveToKeyfile("Luminance Curve", "LCredsk", labCurve.lcredsk, keyFile);
            saveToKeyfile("Luminance Curve", "LCurve", labCurve.lcurve, keyFile);
            saveToKeyfile("Luminance Curve", "aCurve", labCurve.acurve, keyFile);
            saveToKeyfile("Luminance Curve", "bCurve", labCurve.bcurve, keyFile);
            saveToKeyfile("Luminance Curve", "ccCurve", labCurve.cccurve, keyFile);
            saveToKeyfile("Luminance Curve", "chCurve", labCurve.chcurve, keyFile);
            saveToKeyfile("Luminance Curve", "lhCurve", labCurve.lhcurve, keyFile);
            saveToKeyfile("Luminance Curve", "hhCurve", labCurve.hhcurve, keyFile);
            saveToKeyfile("Luminance Curve", "LcCurve", labCurve.lccurve, keyFile);
            saveToKeyfile("Luminance Curve", "ClCurve", labCurve.clcurve, keyFile);
        }

// Sharpening
        if (RELEVANT_(sharpening)) {
            saveToKeyfile("Sharpening", "Enabled", sharpening.enabled, keyFile);
            saveToKeyfile("Sharpening", "Contrast", sharpening.contrast, keyFile);
            saveToKeyfile("Sharpening", "Method", sharpening.method, keyFile);
            saveToKeyfile("Sharpening", "Radius", sharpening.radius, keyFile);
            saveToKeyfile("Sharpening", "BlurRadius", sharpening.blurradius, keyFile);
            saveToKeyfile("Sharpening", "Amount", sharpening.amount, keyFile);
            saveToKeyfile("Sharpening", "Threshold", sharpening.threshold.toVector(), keyFile);
            saveToKeyfile("Sharpening", "OnlyEdges", sharpening.edgesonly, keyFile);
            saveToKeyfile("Sharpening", "EdgedetectionRadius", sharpening.edges_radius, keyFile);
            saveToKeyfile("Sharpening", "EdgeTolerance", sharpening.edges_tolerance, keyFile);
            saveToKeyfile("Sharpening", "HalocontrolEnabled", sharpening.halocontrol, keyFile);
            saveToKeyfile("Sharpening", "HalocontrolAmount", sharpening.halocontrol_amount, keyFile);
            saveToKeyfile("Sharpening", "DeconvRadius", sharpening.deconvradius, keyFile);
            saveToKeyfile("Sharpening", "DeconvAmount", sharpening.deconvamount, keyFile);
            saveToKeyfile("Sharpening", "DeconvDamping", sharpening.deconvdamping, keyFile);
            saveToKeyfile("Sharpening", "DeconvIterations", sharpening.deconviter, keyFile);
        }

// Micro-contrast sharpening
        if (RELEVANT_(sharpenMicro)) {
            saveToKeyfile("SharpenMicro", "Enabled", sharpenMicro.enabled, keyFile);
            saveToKeyfile("SharpenMicro", "Matrix", sharpenMicro.matrix, keyFile);
            saveToKeyfile("SharpenMicro", "Strength", sharpenMicro.amount, keyFile);
            saveToKeyfile("SharpenMicro", "Contrast", sharpenMicro.contrast, keyFile);
            saveToKeyfile("SharpenMicro", "Uniformity", sharpenMicro.uniformity, keyFile);
        }

// WB
        if (RELEVANT_(wb)) {
            saveToKeyfile("White Balance", "Enabled", wb.enabled, keyFile);
            saveToKeyfile("White Balance", "Setting", wb.method, keyFile);
            saveToKeyfile("White Balance", "Temperature", wb.temperature, keyFile);
            saveToKeyfile("White Balance", "Green", wb.green, keyFile);
            saveToKeyfile("White Balance", "Equal", wb.equal, keyFile);
            saveToKeyfile("White Balance", "TemperatureBias", wb.tempBias, keyFile);
        }


// Impulse denoise
        if (RELEVANT_(impulseDenoise)) {
            saveToKeyfile("Impulse Denoising", "Enabled", impulseDenoise.enabled, keyFile);
            saveToKeyfile("Impulse Denoising", "Threshold", impulseDenoise.thresh, keyFile);
        }

// Defringe
        if (RELEVANT_(defringe)) {
            saveToKeyfile("Defringing", "Enabled", defringe.enabled, keyFile);
            saveToKeyfile("Defringing", "Radius", defringe.radius, keyFile);
            saveToKeyfile("Defringing", "Threshold", defringe.threshold, keyFile);
            saveToKeyfile("Defringing", "HueCurve", defringe.huecurve, keyFile);
        }

// Dehaze
        if (RELEVANT_(dehaze)) {
            saveToKeyfile("Dehaze", "Enabled", dehaze.enabled, keyFile);
            saveToKeyfile("Dehaze", "Strength", dehaze.strength, keyFile);        
            saveToKeyfile("Dehaze", "ShowDepthMap", dehaze.showDepthMap, keyFile);        
            saveToKeyfile("Dehaze", "Depth", dehaze.depth, keyFile);
            saveToKeyfile("Dehaze", "Luminance", dehaze.luminance, keyFile);
        }

// Denoising
        if (RELEVANT_(denoise)) {
            saveToKeyfile("Denoise", "Enabled", denoise.enabled, keyFile);
            saveToKeyfile("Denoise", "ColorSpace", Glib::ustring(denoise.colorSpace == DenoiseParams::ColorSpace::LAB ? "Lab" : "RGB"), keyFile);
            saveToKeyfile("Denoise", "Aggressive", denoise.aggressive, keyFile);
            saveToKeyfile("Denoise", "Gamma", denoise.gamma, keyFile);
            saveToKeyfile("Denoise", "Luminance", denoise.luminance, keyFile);
            saveToKeyfile("Denoise", "LuminanceDetail", denoise.luminanceDetail, keyFile);
            saveToKeyfile("Denoise", "ChrominanceMethod", int(denoise.chrominanceMethod), keyFile);
            saveToKeyfile("Denoise", "ChrominanceAutoFactor", denoise.chrominanceAutoFactor, keyFile);
            saveToKeyfile("Denoise", "Chrominance", denoise.chrominance, keyFile);
            saveToKeyfile("Denoise", "ChrominanceRedGreen", denoise.chrominanceRedGreen, keyFile);
            saveToKeyfile("Denoise", "ChrominanceBlueYellow", denoise.chrominanceBlueYellow, keyFile);
            saveToKeyfile("Denoise", "SmoothingEnabled", denoise.smoothingEnabled, keyFile);
            saveToKeyfile("Denoise", "SmoothingMethod", int(denoise.smoothingMethod), keyFile);
            saveToKeyfile("Denoise", "MedianType", int(denoise.medianType), keyFile);
            saveToKeyfile("Denoise", "MedianMethod", int(denoise.medianMethod), keyFile);
            saveToKeyfile("Denoise", "MedianIterations", denoise.medianIterations, keyFile);
            saveToKeyfile("Denoise", "GuidedLumaRadius", denoise.guidedLumaRadius, keyFile);
            saveToKeyfile("Denoise", "GuidedChromaRadius", denoise.guidedChromaRadius, keyFile);
            saveToKeyfile("Denoise", "GuidedLumaStrength", denoise.guidedLumaStrength, keyFile);
            saveToKeyfile("Denoise", "GuidedChromaStrength", denoise.guidedChromaStrength, keyFile);
        }

// EPD
        if (RELEVANT_(epd)) {
            saveToKeyfile("EPD", "Enabled", epd.enabled, keyFile);
            for (size_t j = 0; j < epd.regions.size(); ++j) {
                std::string n = j ? std::string("_") + std::to_string(j) : std::string("");
                auto &r = epd.regions[j];
                putToKeyfile("EPD", Glib::ustring("Strength") + n, r.strength, keyFile);
                putToKeyfile("EPD", Glib::ustring("Gamma") + n, r.gamma, keyFile);
                putToKeyfile("EPD", Glib::ustring("EdgeStopping") + n, r.edgeStopping, keyFile);
                putToKeyfile("EPD", Glib::ustring("Scale") + n, r.scale, keyFile);
                putToKeyfile("EPD", Glib::ustring("ReweightingIterates") + n, r.reweightingIterates, keyFile);
                epd.labmasks[j].save(keyFile, "EPD", "", n);
            }
            saveToKeyfile("EPD", "showMask", epd.showMask, keyFile);
        }

// Fattal
        if (RELEVANT_(fattal)) {
            saveToKeyfile("FattalToneMapping", "Enabled", fattal.enabled, keyFile);
            saveToKeyfile("FattalToneMapping", "Threshold", fattal.threshold, keyFile);
            saveToKeyfile("FattalToneMapping", "Amount", fattal.amount, keyFile);
            saveToKeyfile("FattalToneMapping", "Anchor", fattal.anchor, keyFile);
        }

// Log encoding
        if (RELEVANT_(logenc)) {
            saveToKeyfile("LogEncoding", "Enabled", logenc.enabled, keyFile);
            saveToKeyfile("LogEncoding", "Auto", logenc.autocompute, keyFile);
            saveToKeyfile("LogEncoding", "AutoGray", logenc.autogray, keyFile);
            saveToKeyfile("LogEncoding", "SourceGray", logenc.sourceGray, keyFile);
            saveToKeyfile("LogEncoding", "TargetGray", logenc.targetGray, keyFile);
            saveToKeyfile("LogEncoding", "BlackEv", logenc.blackEv, keyFile);
            saveToKeyfile("LogEncoding", "WhiteEv", logenc.whiteEv, keyFile);
            saveToKeyfile("LogEncoding", "Detail", logenc.detail, keyFile);
        }

// Shadows & highlights
        if (RELEVANT_(sh)) {
            saveToKeyfile("Shadows & Highlights", "Enabled", sh.enabled, keyFile);
            saveToKeyfile("Shadows & Highlights", "Highlights", sh.highlights, keyFile);
            saveToKeyfile("Shadows & Highlights", "HighlightTonalWidth", sh.htonalwidth, keyFile);
            saveToKeyfile("Shadows & Highlights", "Shadows", sh.shadows, keyFile);
            saveToKeyfile("Shadows & Highlights", "ShadowTonalWidth", sh.stonalwidth, keyFile);
            saveToKeyfile("Shadows & Highlights", "Radius", sh.radius, keyFile);
            saveToKeyfile("Shadows & Highlights", "Lab", sh.lab, keyFile);
        }

// ToneEqualizer
        if (RELEVANT_(toneEqualizer)) {
            saveToKeyfile("ToneEqualizer", "Enabled", toneEqualizer.enabled, keyFile);
            for (size_t i = 0; i < toneEqualizer.bands.size(); ++i) {
                saveToKeyfile("ToneEqualizer", "Band" + std::to_string(i), toneEqualizer.bands[i], keyFile);
            }
            saveToKeyfile("ToneEqualizer", "Detail", toneEqualizer.detail, keyFile);
        }
        
// Crop
        if (RELEVANT_(crop)) {
            saveToKeyfile("Crop", "Enabled", crop.enabled, keyFile);
            saveToKeyfile("Crop", "X", crop.x, keyFile);
            saveToKeyfile("Crop", "Y", crop.y, keyFile);
            saveToKeyfile("Crop", "W", crop.w, keyFile);
            saveToKeyfile("Crop", "H", crop.h, keyFile);
            saveToKeyfile("Crop", "FixedRatio", crop.fixratio, keyFile);
            saveToKeyfile("Crop", "Ratio", crop.ratio, keyFile);
            saveToKeyfile("Crop", "Orientation", crop.orientation, keyFile);
            saveToKeyfile("Crop", "Guide", crop.guide, keyFile);
        }

// Coarse transformation
        if (RELEVANT_(coarse)) {
            saveToKeyfile("Coarse Transformation", "Rotate", coarse.rotate, keyFile);
            saveToKeyfile("Coarse Transformation", "HorizontalFlip", coarse.hflip, keyFile);
            saveToKeyfile("Coarse Transformation", "VerticalFlip", coarse.vflip, keyFile);
        }

// Common properties for transformations
        if (RELEVANT_(commonTrans)) {
            saveToKeyfile("Common Properties for Transformations", "AutoFill", commonTrans.autofill, keyFile);
        }

// Rotation
        if (RELEVANT_(rotate)) {
            saveToKeyfile("Rotation", "Enabled", rotate.enabled, keyFile);
            saveToKeyfile("Rotation", "Degree", rotate.degree, keyFile);
        }

// Distortion
        if (RELEVANT_(distortion)) {
            saveToKeyfile("Distortion", "Enabled", distortion.enabled, keyFile);
            saveToKeyfile("Distortion", "Amount", distortion.amount, keyFile);
        }

// Lens profile
        if (RELEVANT_(lensProf)) {
            saveToKeyfile("LensProfile", "LcMode", lensProf.getMethodString(lensProf.lcMode), keyFile);
            saveToKeyfile("LensProfile", "LCPFile", relativePathIfInside(fname, fnameAbsolute, lensProf.lcpFile), keyFile);
            saveToKeyfile("LensProfile", "UseDistortion", lensProf.useDist, keyFile);
            saveToKeyfile("LensProfile", "UseVignette", lensProf.useVign, keyFile);
            saveToKeyfile("LensProfile", "UseCA", lensProf.useCA, keyFile);
            saveToKeyfile("LensProfile", "LFCameraMake", lensProf.lfCameraMake, keyFile);
            saveToKeyfile("LensProfile", "LFCameraModel", lensProf.lfCameraModel, keyFile);
            saveToKeyfile("LensProfile", "LFLens", lensProf.lfLens, keyFile);
        }

// Perspective correction
        if (RELEVANT_(perspective)) {
            saveToKeyfile("Perspective", "Enabled", perspective.enabled, keyFile);
            saveToKeyfile("Perspective", "Horizontal", perspective.horizontal, keyFile);
            saveToKeyfile("Perspective", "Vertical", perspective.vertical, keyFile);
            saveToKeyfile("Perspective", "Angle", perspective.angle, keyFile);
            saveToKeyfile("Perspective", "Shear", perspective.shear, keyFile);
            saveToKeyfile("Perspective", "FocalLength", perspective.flength, keyFile);
            saveToKeyfile("Perspective", "CropFactor", perspective.cropfactor, keyFile);
            saveToKeyfile("Perspective", "Aspect", perspective.aspect, keyFile);
        }

// Gradient
        if (RELEVANT_(gradient)) {
            saveToKeyfile("Gradient", "Enabled", gradient.enabled, keyFile);
            saveToKeyfile("Gradient", "Degree", gradient.degree, keyFile);
            saveToKeyfile("Gradient", "Feather", gradient.feather, keyFile);
            saveToKeyfile("Gradient", "Strength", gradient.strength, keyFile);
            saveToKeyfile("Gradient", "CenterX", gradient.centerX, keyFile);
            saveToKeyfile("Gradient", "CenterY", gradient.centerY, keyFile);
        }

// Post-crop vignette
        if (RELEVANT_(pcvignette)) {
            saveToKeyfile("PCVignette", "Enabled", pcvignette.enabled, keyFile);
            saveToKeyfile("PCVignette", "Strength", pcvignette.strength, keyFile);
            saveToKeyfile("PCVignette", "Feather", pcvignette.feather, keyFile);
            saveToKeyfile("PCVignette", "Roundness", pcvignette.roundness, keyFile);
        }

// C/A correction
        if (RELEVANT_(cacorrection)) {
            saveToKeyfile("CACorrection", "Enabled", cacorrection.enabled, keyFile);
            saveToKeyfile("CACorrection", "Red", cacorrection.red, keyFile);
            saveToKeyfile("CACorrection", "Blue", cacorrection.blue, keyFile);
        }

// Vignetting correction
        if (RELEVANT_(vignetting)) {
            saveToKeyfile("Vignetting Correction", "Enabled", vignetting.enabled, keyFile);
            saveToKeyfile("Vignetting Correction", "Amount", vignetting.amount, keyFile);
            saveToKeyfile("Vignetting Correction", "Radius", vignetting.radius, keyFile);
            saveToKeyfile("Vignetting Correction", "Strength", vignetting.strength, keyFile);
            saveToKeyfile("Vignetting Correction", "CenterX", vignetting.centerX, keyFile);
            saveToKeyfile("Vignetting Correction", "CenterY", vignetting.centerY, keyFile);
        }

// Resize
        if (RELEVANT_(resize)) {
            saveToKeyfile("Resize", "Enabled", resize.enabled, keyFile);
            saveToKeyfile("Resize", "Scale", resize.scale, keyFile);
            saveToKeyfile("Resize", "AppliesTo", resize.appliesTo, keyFile);
            saveToKeyfile("Resize", "Method", resize.method, keyFile);
            saveToKeyfile("Resize", "DataSpecified", resize.dataspec, keyFile);
            saveToKeyfile("Resize", "Width", resize.width, keyFile);
            saveToKeyfile("Resize", "Height", resize.height, keyFile);
            saveToKeyfile("Resize", "AllowUpscaling", resize.allowUpscaling, keyFile);
        }

// Post resize sharpening
        if (RELEVANT_(prsharpening)) {
            saveToKeyfile("PostResizeSharpening", "Enabled", prsharpening.enabled, keyFile);
            saveToKeyfile("PostResizeSharpening", "Contrast", prsharpening.contrast, keyFile);
            saveToKeyfile("PostResizeSharpening", "Method", prsharpening.method, keyFile);
            saveToKeyfile("PostResizeSharpening", "Radius", prsharpening.radius, keyFile);
            saveToKeyfile("PostResizeSharpening", "Amount", prsharpening.amount, keyFile);
            saveToKeyfile("PostResizeSharpening", "Threshold", prsharpening.threshold.toVector(), keyFile);
            saveToKeyfile("PostResizeSharpening", "OnlyEdges", prsharpening.edgesonly, keyFile);
            saveToKeyfile("PostResizeSharpening", "EdgedetectionRadius", prsharpening.edges_radius, keyFile);
            saveToKeyfile("PostResizeSharpening", "EdgeTolerance", prsharpening.edges_tolerance, keyFile);
            saveToKeyfile("PostResizeSharpening", "HalocontrolEnabled", prsharpening.halocontrol, keyFile);
            saveToKeyfile("PostResizeSharpening", "HalocontrolAmount", prsharpening.halocontrol_amount, keyFile);
            saveToKeyfile("PostResizeSharpening", "DeconvRadius", prsharpening.deconvradius, keyFile);
            saveToKeyfile("PostResizeSharpening", "DeconvAmount", prsharpening.deconvamount, keyFile);
            saveToKeyfile("PostResizeSharpening", "DeconvDamping", prsharpening.deconvdamping, keyFile);
            saveToKeyfile("PostResizeSharpening", "DeconvIterations", prsharpening.deconviter, keyFile);
        }

// Color management
        if (RELEVANT_(icm)) {
            saveToKeyfile("Color Management", "InputProfile", relativePathIfInside(fname, fnameAbsolute, icm.inputProfile), keyFile);
            saveToKeyfile("Color Management", "ToneCurve", icm.toneCurve, keyFile);
            saveToKeyfile("Color Management", "ApplyLookTable", icm.applyLookTable, keyFile);
            saveToKeyfile("Color Management", "ApplyBaselineExposureOffset", icm.applyBaselineExposureOffset, keyFile);
            saveToKeyfile("Color Management", "ApplyHueSatMap", icm.applyHueSatMap, keyFile);
            saveToKeyfile("Color Management", "DCPIlluminant", icm.dcpIlluminant, keyFile);
            saveToKeyfile("Color Management", "WorkingProfile", icm.workingProfile, keyFile);
            saveToKeyfile("Color Management", "OutputProfile", icm.outputProfile, keyFile);
            saveToKeyfile(
                "Color Management",
                "OutputProfileIntent",
                {
                    {RI_PERCEPTUAL, "Perceptual"},
                    {RI_RELATIVE, "Relative"},
                    {RI_SATURATION, "Saturation"},
                    {RI_ABSOLUTE, "Absolute"}

                },
                icm.outputIntent,
                keyFile
                );
            saveToKeyfile("Color Management", "OutputBPC", icm.outputBPC, keyFile);
        }


// Directional pyramid equalizer
        if (RELEVANT_(dirpyrequalizer)) {
            saveToKeyfile("Directional Pyramid Equalizer", "Enabled", dirpyrequalizer.enabled, keyFile);
            for (size_t j = 0; j < dirpyrequalizer.levels.size(); ++j) {
                std::string n = std::to_string(j+1);
                auto &l = dirpyrequalizer.levels[j];
                for (int i = 0; i < 6; ++i) {
                    putToKeyfile("Directional Pyramid Equalizer", Glib::ustring("Mult") + std::to_string(i) + "_" + n, l.mult[i], keyFile);
                }
                putToKeyfile("Directional Pyramid Equalizer", Glib::ustring("Threshold_") + n, l.threshold, keyFile);
                dirpyrequalizer.labmasks[j].save(keyFile, "Directional Pyramid Equalizer", "", Glib::ustring("_") + n);
            }
            saveToKeyfile("Directional Pyramid Equalizer", "ShowMask", dirpyrequalizer.showMask, keyFile);
        }

// Soft Light
        if (RELEVANT_(softlight)) {
            saveToKeyfile("SoftLight", "Enabled", softlight.enabled, keyFile);
            saveToKeyfile("SoftLight", "Strength", softlight.strength, keyFile);
        }

// Film simulation
        if (RELEVANT_(filmSimulation)) {
            saveToKeyfile("Film Simulation", "Enabled", filmSimulation.enabled, keyFile);
            saveToKeyfile("Film Simulation", "ClutFilename", filmSimulation.clutFilename, keyFile);
            saveToKeyfile("Film Simulation", "Strength", filmSimulation.strength, keyFile);
        }

// RGB curves        
        if (RELEVANT_(rgbCurves)) {
            saveToKeyfile("RGB Curves", "Enabled", rgbCurves.enabled, keyFile);
            saveToKeyfile("RGB Curves", "LumaMode", rgbCurves.lumamode, keyFile);
            saveToKeyfile("RGB Curves", "rCurve", rgbCurves.rcurve, keyFile);
            saveToKeyfile("RGB Curves", "gCurve", rgbCurves.gcurve, keyFile);
            saveToKeyfile("RGB Curves", "bCurve", rgbCurves.bcurve, keyFile);
        }

// Grain
        if (RELEVANT_(grain)) {
            saveToKeyfile("Grain", "Enabled", grain.enabled, keyFile);
            saveToKeyfile("Grain", "ISO", grain.iso, keyFile);
            saveToKeyfile("Grain", "Strength", grain.strength, keyFile);
            saveToKeyfile("Grain", "Scale", grain.scale, keyFile);
        }


// Smoothing
        if (RELEVANT_(smoothing)) {
            saveToKeyfile("GuidedSmoothing", "Enabled", smoothing.enabled, keyFile);
            for (size_t j = 0; j < smoothing.regions.size(); ++j) {
                std::string n = std::to_string(j+1);
                auto &r = smoothing.regions[j];
                putToKeyfile("GuidedSmoothing", Glib::ustring("Channel_") + n, int(r.channel), keyFile);
                putToKeyfile("GuidedSmoothing", Glib::ustring("Radius_") + n, r.radius, keyFile);
                putToKeyfile("GuidedSmoothing", Glib::ustring("Epsilon_") + n, r.epsilon, keyFile);
                smoothing.labmasks[j].save(keyFile, "GuidedSmoothing", "", Glib::ustring("_") + n);
            }
            saveToKeyfile("GuidedSmoothing", "ShowMask", smoothing.showMask, keyFile);
        }

// ColorCorrection
        if (RELEVANT_(colorcorrection)) {
            saveToKeyfile("ColorCorrection", "Enabled", colorcorrection.enabled, keyFile);
            for (size_t j = 0; j < colorcorrection.regions.size(); ++j) {
                std::string n = std::to_string(j+1);
                auto &l = colorcorrection.regions[j];
                putToKeyfile("ColorCorrection", Glib::ustring("A_") + n, l.a, keyFile);
                putToKeyfile("ColorCorrection", Glib::ustring("B_") + n, l.b, keyFile);
                putToKeyfile("ColorCorrection", Glib::ustring("Saturation_") + n, l.saturation, keyFile);
                putToKeyfile("ColorCorrection", Glib::ustring("Slope_") + n, l.slope, keyFile);
                putToKeyfile("ColorCorrection", Glib::ustring("Offset_") + n, l.offset, keyFile);
                putToKeyfile("ColorCorrection", Glib::ustring("Power_") + n, l.power, keyFile);
                putToKeyfile("ColorCorrection", Glib::ustring("Channel_") + n, l.channel, keyFile);
                colorcorrection.labmasks[j].save(keyFile, "ColorCorrection", "", Glib::ustring("_") + n);
            }
            saveToKeyfile("ColorCorrection", "showMask", colorcorrection.showMask, keyFile);
        }
        
// Raw
        if (RELEVANT_(darkframe)) {
            saveToKeyfile("RAW", "DarkFrameEnabled", raw.enable_darkframe, keyFile);
            saveToKeyfile("RAW", "DarkFrame", relativePathIfInside(fname, fnameAbsolute, raw.dark_frame), keyFile);
            saveToKeyfile("RAW", "DarkFrameAuto", raw.df_autoselect, keyFile);
        }
        if (RELEVANT_(flatfield)) {
            saveToKeyfile("RAW", "FlatFieldEnabled", raw.enable_flatfield, keyFile);
            saveToKeyfile("RAW", "FlatFieldFile", relativePathIfInside(fname, fnameAbsolute, raw.ff_file), keyFile);
            saveToKeyfile("RAW", "FlatFieldAutoSelect", raw.ff_AutoSelect, keyFile);
            saveToKeyfile("RAW", "FlatFieldBlurRadius", raw.ff_BlurRadius, keyFile);
            saveToKeyfile("RAW", "FlatFieldBlurType", raw.ff_BlurType, keyFile);
            saveToKeyfile("RAW", "FlatFieldAutoClipControl", raw.ff_AutoClipControl, keyFile);
            saveToKeyfile("RAW", "FlatFieldClipControl", raw.ff_clipControl, keyFile);
        }
        if (RELEVANT_(rawCA)) {
            saveToKeyfile("RAW", "CAEnabled", raw.enable_ca, keyFile);
            saveToKeyfile("RAW", "CA", raw.ca_autocorrect, keyFile);
            saveToKeyfile("RAW", "CAAvoidColourshift", raw.ca_avoidcolourshift, keyFile);
            saveToKeyfile("RAW", "CAAutoIterations", raw.caautoiterations, keyFile);
            saveToKeyfile("RAW", "CARed", raw.cared, keyFile);
            saveToKeyfile("RAW", "CABlue", raw.cablue, keyFile);
        }
        if (RELEVANT_(hotDeadPixelFilter)) {
            saveToKeyfile("RAW", "HotDeadPixelEnabled", raw.enable_hotdeadpix, keyFile);
            saveToKeyfile("RAW", "HotPixelFilter", raw.hotPixelFilter, keyFile);
            saveToKeyfile("RAW", "DeadPixelFilter", raw.deadPixelFilter, keyFile);
            saveToKeyfile("RAW", "HotDeadPixelThresh", raw.hotdeadpix_thresh, keyFile);
        }
        if (RELEVANT_(demosaic)) {
            saveToKeyfile("RAW Bayer", "Method", raw.bayersensor.method, keyFile);
            saveToKeyfile("RAW Bayer", "Border", raw.bayersensor.border, keyFile);
            saveToKeyfile("RAW Bayer", "ImageNum", raw.bayersensor.imageNum + 1, keyFile);
            saveToKeyfile("RAW Bayer", "CcSteps", raw.bayersensor.ccSteps, keyFile);
        }
        if (RELEVANT_(rawBlack)) {
            saveToKeyfile("RAW Bayer", "PreBlackEnabled", raw.bayersensor.enable_black, keyFile);
            saveToKeyfile("RAW Bayer", "PreBlack0", raw.bayersensor.black0, keyFile);
            saveToKeyfile("RAW Bayer", "PreBlack1", raw.bayersensor.black1, keyFile);
            saveToKeyfile("RAW Bayer", "PreBlack2", raw.bayersensor.black2, keyFile);
            saveToKeyfile("RAW Bayer", "PreBlack3", raw.bayersensor.black3, keyFile);
            saveToKeyfile("RAW Bayer", "PreTwoGreen", raw.bayersensor.twogreen, keyFile);
        }
        if (RELEVANT_(rawPreprocessing)) {
            saveToKeyfile("RAW Bayer", "PreprocessingEnabled", raw.bayersensor.enable_preproc, keyFile);
            saveToKeyfile("RAW Bayer", "LineDenoise", raw.bayersensor.linenoise, keyFile);
            saveToKeyfile("RAW Bayer", "LineDenoiseDirection", toUnderlying(raw.bayersensor.linenoiseDirection), keyFile);
            saveToKeyfile("RAW Bayer", "GreenEqThreshold", raw.bayersensor.greenthresh, keyFile);
        }
        if (RELEVANT_(demosaic)) {
            saveToKeyfile("RAW Bayer", "DCBIterations", raw.bayersensor.dcb_iterations, keyFile);
            saveToKeyfile("RAW Bayer", "DCBEnhance", raw.bayersensor.dcb_enhance, keyFile);
            saveToKeyfile("RAW Bayer", "LMMSEIterations", raw.bayersensor.lmmse_iterations, keyFile);
            saveToKeyfile("RAW Bayer", "DualDemosaicAutoContrast", raw.bayersensor.dualDemosaicAutoContrast, keyFile);
            saveToKeyfile("RAW Bayer", "DualDemosaicContrast", raw.bayersensor.dualDemosaicContrast, keyFile);
            saveToKeyfile("RAW Bayer", "PixelShiftMotionCorrectionMethod", toUnderlying(raw.bayersensor.pixelShiftMotionCorrectionMethod), keyFile);
            saveToKeyfile("RAW Bayer", "PixelShiftEperIso", raw.bayersensor.pixelShiftEperIso, keyFile);
            saveToKeyfile("RAW Bayer", "PixelShiftSigma", raw.bayersensor.pixelShiftSigma, keyFile);
            saveToKeyfile("RAW Bayer", "PixelShiftShowMotion", raw.bayersensor.pixelShiftShowMotion, keyFile);
            saveToKeyfile("RAW Bayer", "PixelShiftShowMotionMaskOnly", raw.bayersensor.pixelShiftShowMotionMaskOnly, keyFile);
            saveToKeyfile("RAW Bayer", "pixelShiftHoleFill", raw.bayersensor.pixelShiftHoleFill, keyFile);
            saveToKeyfile("RAW Bayer", "pixelShiftMedian", raw.bayersensor.pixelShiftMedian, keyFile);
            saveToKeyfile("RAW Bayer", "pixelShiftGreen", raw.bayersensor.pixelShiftGreen, keyFile);
            saveToKeyfile("RAW Bayer", "pixelShiftBlur", raw.bayersensor.pixelShiftBlur, keyFile);
            saveToKeyfile("RAW Bayer", "pixelShiftSmoothFactor", raw.bayersensor.pixelShiftSmoothFactor, keyFile);
            saveToKeyfile("RAW Bayer", "pixelShiftEqualBright", raw.bayersensor.pixelShiftEqualBright, keyFile);
            saveToKeyfile("RAW Bayer", "pixelShiftEqualBrightChannel", raw.bayersensor.pixelShiftEqualBrightChannel, keyFile);
            saveToKeyfile("RAW Bayer", "pixelShiftNonGreenCross", raw.bayersensor.pixelShiftNonGreenCross, keyFile);
            saveToKeyfile("RAW Bayer", "pixelShiftDemosaicMethod", raw.bayersensor.pixelShiftDemosaicMethod, keyFile);
        }
        if (RELEVANT_(rawPreprocessing)) {
            saveToKeyfile("RAW Bayer", "PDAFLinesFilter", raw.bayersensor.pdafLinesFilter, keyFile);
        }
        if (RELEVANT_(demosaic)) {
            saveToKeyfile("RAW X-Trans", "Method", raw.xtranssensor.method, keyFile);
            saveToKeyfile("RAW X-Trans", "DualDemosaicAutoContrast", raw.xtranssensor.dualDemosaicAutoContrast, keyFile);
            saveToKeyfile("RAW X-Trans", "DualDemosaicContrast", raw.xtranssensor.dualDemosaicContrast, keyFile);
            saveToKeyfile("RAW X-Trans", "Border", raw.xtranssensor.border, keyFile);
            saveToKeyfile("RAW X-Trans", "CcSteps", raw.xtranssensor.ccSteps, keyFile);
        }
        if (RELEVANT_(rawBlack)) {
            saveToKeyfile("RAW X-Trans", "PreBlackEnabled", raw.xtranssensor.enable_black, keyFile);
            saveToKeyfile("RAW X-Trans", "PreBlackRed", raw.xtranssensor.blackred, keyFile);
            saveToKeyfile("RAW X-Trans", "PreBlackGreen", raw.xtranssensor.blackgreen, keyFile);
            saveToKeyfile("RAW X-Trans", "PreBlackBlue", raw.xtranssensor.blackblue, keyFile);
        }

// Raw exposition
        if (RELEVANT_(rawWhite)) {
            saveToKeyfile("RAW", "PreExposureEnabled", raw.enable_whitepoint, keyFile);
            saveToKeyfile("RAW", "PreExposure", raw.expos, keyFile);
        }

// MetaData
        if (RELEVANT_(metadata)) {
            saveToKeyfile("MetaData", "Mode", metadata.mode, keyFile);
        }

// EXIF change list
        if (RELEVANT_(exif)) {
            std::map<Glib::ustring, Glib::ustring> m;
            for (auto &p : exif_keys) {
                m[p.second] = p.first;
            }
            for (ExifPairs::const_iterator i = exif.begin(); i != exif.end(); ++i) {
                auto it = m.find(i->first);
                if (it != m.end()) {
                    keyFile.set_string("Exif", it->second, i->second);
                }
            }
        }

// IPTC change list
        if (RELEVANT_(iptc)) {
            std::map<Glib::ustring, Glib::ustring> m;
            for (auto &p : iptc_keys) {
                m[p.second] = p.first;
            }
            for (IPTCPairs::const_iterator i = iptc.begin(); i != iptc.end(); ++i) {
                auto it = m.find(i->first);
                if (it != m.end()) {
                    Glib::ArrayHandle<Glib::ustring> values = i->second;
                    keyFile.set_string_list("IPTC", it->second, values);
                }
            }
        }
    } catch (Glib::KeyFileError&) {
        return 1;
    }

    return 0;
#undef RELEVANT_
}


int ProcParams::load(const Glib::ustring& fname, const ParamsEdited *pedited)
{
    setlocale(LC_NUMERIC, "C");  // to set decimal point to "."

    if (fname.empty()) {
        return 1;
    }

    Glib::KeyFile keyFile;

    try {
        if (!Glib::file_test(fname, Glib::FILE_TEST_EXISTS) ||
                !keyFile.load_from_file(fname)) {
            return 1;
        }

        return load(keyFile, pedited, true, fname);
    } catch (const Glib::Error& e) {
        printf("-->%s\n", e.what().c_str());
        setDefaults();
        return 1;
    } catch (...) {
        printf("-->unknown exception!\n");
        setDefaults();
        return 1;
    }
}


int ProcParams::load(const Glib::KeyFile &keyFile, const ParamsEdited *pedited,
                     bool resetOnError, const Glib::ustring &fname)
{
#define RELEVANT_(n) (!pedited || pedited->n)
    
    try {
        ppVersion = PPVERSION;
        appVersion = RTVERSION;

        if (keyFile.has_group("Version")) {
            if (keyFile.has_key("Version", "AppVersion")) {
                appVersion = keyFile.get_string("Version", "AppVersion");
            }

            if (keyFile.has_key("Version", "Version")) {
                ppVersion = keyFile.get_integer("Version", "Version");
            }
        }

        if (keyFile.has_group("General") && RELEVANT_(general)) {
            assignFromKeyfile(keyFile, "General", "Rank", rank);
            assignFromKeyfile(keyFile, "General", "ColorLabel", colorlabel);
            assignFromKeyfile(keyFile, "General", "InTrash", inTrash);
        }

        const std::map<std::string, ToneCurveParams::TcMode> tc_mapping = {
            {"Standard", ToneCurveParams::TcMode::STD},
            {"FilmLike", ToneCurveParams::TcMode::FILMLIKE},
            {"SatAndValueBlending", ToneCurveParams::TcMode::SATANDVALBLENDING},
            {"WeightedStd", ToneCurveParams::TcMode::WEIGHTEDSTD},
            {"Luminance", ToneCurveParams::TcMode::LUMINANCE},
            {"Perceptual", ToneCurveParams::TcMode::PERCEPTUAL}
        };

        if (ppVersion < 350) {
            if (keyFile.has_group("Exposure")) {
                if (RELEVANT_(exposure)) {
                    exposure.enabled = true;
                    
                    if (ppVersion < PPVERSION_AEXP) {
                        exposure.autoexp = false; // prevent execution of autoexp when opening file created with earlier versions of autoexp algorithm
                    } else {
                        assignFromKeyfile(keyFile, "Exposure", "Auto", exposure.autoexp);
                    }

                    assignFromKeyfile(keyFile, "Exposure", "Clip", exposure.clip);
                    assignFromKeyfile(keyFile, "Exposure", "Compensation", exposure.expcomp);
                    assignFromKeyfile(keyFile, "Exposure", "Black", exposure.black);
                    assignFromKeyfile(keyFile, "Exposure", "HighlightCompr", exposure.hlcompr);
                    assignFromKeyfile(keyFile, "Exposure", "HighlightComprThreshold", exposure.hlcomprthresh);
                    assignFromKeyfile(keyFile, "Exposure", "ShadowCompr", exposure.shcompr);
                    assignFromKeyfile(keyFile, "Exposure", "ClampOOG", exposure.clampOOG);

                    if (exposure.shcompr > 100) {
                        exposure.shcompr = 100; // older pp3 files can have values above 100.
                    }
                }

                if (RELEVANT_(toneCurve)) {
                    toneCurve.enabled = true;
                    
                    assignFromKeyfile(keyFile, "Exposure", "CurveMode", tc_mapping, toneCurve.curveMode);
                    assignFromKeyfile(keyFile, "Exposure", "CurveMode2", tc_mapping, toneCurve.curveMode2);

                    if (ppVersion > 200) {
                        assignFromKeyfile(keyFile, "Exposure", "Curve", toneCurve.curve);
                        assignFromKeyfile(keyFile, "Exposure", "Curve2", toneCurve.curve2);
                    }

                    assignFromKeyfile(keyFile, "Exposure", "HistogramMatching", toneCurve.histmatching);
                    if (ppVersion < 340) {
                        toneCurve.fromHistMatching = false;
                    } else {
                        assignFromKeyfile(keyFile, "Exposure", "CurveFromHistogramMatching", toneCurve.fromHistMatching);
                    }
                }
                if (RELEVANT_(brightContrSat)) {
                    brightContrSat.enabled = true;
                    assignFromKeyfile(keyFile, "Exposure", "Brightness", brightContrSat.brightness);
                    assignFromKeyfile(keyFile, "Exposure", "Contrast", brightContrSat.contrast);
                    assignFromKeyfile(keyFile, "Exposure", "Saturation", brightContrSat.saturation);
                }
            }
            if (keyFile.has_group("HLRecovery") && RELEVANT_(exposure)) {
                assignFromKeyfile(keyFile, "HLRecovery", "Enabled", exposure.hrenabled);
                assignFromKeyfile(keyFile, "HLRecovery", "Method", exposure.method);
            }
        } else {
            if (keyFile.has_group("Exposure") && RELEVANT_(exposure)) {
                assignFromKeyfile(keyFile, "Exposure", "Enabled", exposure.enabled);
                assignFromKeyfile(keyFile, "Exposure", "Auto", exposure.autoexp);
                assignFromKeyfile(keyFile, "Exposure", "Clip", exposure.clip);
                assignFromKeyfile(keyFile, "Exposure", "Compensation", exposure.expcomp);
                assignFromKeyfile(keyFile, "Exposure", "Black", exposure.black);
                assignFromKeyfile(keyFile, "Exposure", "HighlightCompr", exposure.hlcompr);
                assignFromKeyfile(keyFile, "Exposure", "HighlightComprThreshold", exposure.hlcomprthresh);
                assignFromKeyfile(keyFile, "Exposure", "ShadowCompr", exposure.shcompr);
                assignFromKeyfile(keyFile, "Exposure", "ClampOOG", exposure.clampOOG);
                assignFromKeyfile(keyFile, "Exposure", "HLRecoveryEnabled", exposure.hrenabled);
                assignFromKeyfile(keyFile, "Exposure", "HLRecoveryMethod", exposure.method);
            }
            if (keyFile.has_group("BrightnessContrastSaturation") && RELEVANT_(brightContrSat)) {
                assignFromKeyfile(keyFile, "BrightnessContrastSaturation", "Enabled", brightContrSat.enabled);
                assignFromKeyfile(keyFile, "BrightnessContrastSaturation", "Brightness", brightContrSat.brightness);
                assignFromKeyfile(keyFile, "BrightnessContrastSaturation", "Contrast", brightContrSat.contrast);
                assignFromKeyfile(keyFile, "BrightnessContrastSaturation", "Saturation", brightContrSat.saturation);
            }
            if (keyFile.has_group("ToneCurve") && RELEVANT_(toneCurve)) {
                assignFromKeyfile(keyFile, "ToneCurve", "Enabled", toneCurve.enabled);
                assignFromKeyfile(keyFile, "ToneCurve", "CurveMode", tc_mapping, toneCurve.curveMode);
                assignFromKeyfile(keyFile, "ToneCurve", "CurveMode2", tc_mapping, toneCurve.curveMode2);

                assignFromKeyfile(keyFile, "ToneCurve", "Curve", toneCurve.curve);
                assignFromKeyfile(keyFile, "ToneCurve", "Curve2", toneCurve.curve2);
                assignFromKeyfile(keyFile, "ToneCurve", "HistogramMatching", toneCurve.histmatching);
                assignFromKeyfile(keyFile, "ToneCurve", "CurveFromHistogramMatching", toneCurve.fromHistMatching);
            }
        }

        if (keyFile.has_group("Channel Mixer") && RELEVANT_(chmixer)) {
            if (ppVersion >= 329) {
                assignFromKeyfile(keyFile, "Channel Mixer", "Enabled", chmixer.enabled);
            } else {
                chmixer.enabled = true;
            }

            if (keyFile.has_key("Channel Mixer", "Red") && keyFile.has_key("Channel Mixer", "Green") && keyFile.has_key("Channel Mixer", "Blue")) {
                const std::vector<int> rmix = keyFile.get_integer_list("Channel Mixer", "Red");
                const std::vector<int> gmix = keyFile.get_integer_list("Channel Mixer", "Green");
                const std::vector<int> bmix = keyFile.get_integer_list("Channel Mixer", "Blue");

                if (rmix.size() == 3 && gmix.size() == 3 && bmix.size() == 3) {
                    memcpy(chmixer.red,   rmix.data(), 3 * sizeof(int));
                    memcpy(chmixer.green, gmix.data(), 3 * sizeof(int));
                    memcpy(chmixer.blue,  bmix.data(), 3 * sizeof(int));
                }
                if (ppVersion < 338) {
                    for (int i = 0; i < 3; ++i) {
                        chmixer.red[i] *= 10;
                        chmixer.green[i] *= 10;
                        chmixer.blue[i] *= 10;
                    }
                }
            }
        }

        if (keyFile.has_group("Black & White") && RELEVANT_(blackwhite)) {
            assignFromKeyfile(keyFile, "Black & White", "Enabled", blackwhite.enabled);
            assignFromKeyfile(keyFile, "Black & White", "Method", blackwhite.method);
            assignFromKeyfile(keyFile, "Black & White", "MixerRed", blackwhite.mixerRed);
            assignFromKeyfile(keyFile, "Black & White", "MixerGreen", blackwhite.mixerGreen);
            assignFromKeyfile(keyFile, "Black & White", "MixerBlue", blackwhite.mixerBlue);
            assignFromKeyfile(keyFile, "Black & White", "GammaRed", blackwhite.gammaRed);
            assignFromKeyfile(keyFile, "Black & White", "GammaGreen", blackwhite.gammaGreen);
            assignFromKeyfile(keyFile, "Black & White", "GammaBlue", blackwhite.gammaBlue);
            assignFromKeyfile(keyFile, "Black & White", "Filter", blackwhite.filter);
            assignFromKeyfile(keyFile, "Black & White", "Setting", blackwhite.setting);
            assignFromKeyfile(keyFile, "Black & White", "LuminanceCurve", blackwhite.luminanceCurve);
        }

        if (keyFile.has_group("Local Contrast") && RELEVANT_(localContrast)) {
            assignFromKeyfile(keyFile, "Local Contrast", "Enabled", localContrast.enabled);
            int m = static_cast<int>(LocalContrastParams::USM);
            assignFromKeyfile(keyFile, "Local Contrast", "Mode", m);
            localContrast.mode = LocalContrastParams::Mode(min(max(m, 0), 1));
            assignFromKeyfile(keyFile, "Local Contrast", "Radius", localContrast.radius);
            assignFromKeyfile(keyFile, "Local Contrast", "Amount", localContrast.amount);
            assignFromKeyfile(keyFile, "Local Contrast", "Darkness", localContrast.darkness);
            assignFromKeyfile(keyFile, "Local Contrast", "Lightness", localContrast.lightness);
            assignFromKeyfile(keyFile, "Local Contrast", "Contrast", localContrast.contrast);
            assignFromKeyfile(keyFile, "Local Contrast", "Curve", localContrast.curve);
        }

        if (keyFile.has_group("Luminance Curve") && RELEVANT_(labCurve)) {
            if (ppVersion >= 329) {
                assignFromKeyfile(keyFile, "Luminance Curve", "Enabled", labCurve.enabled);
            } else {
                labCurve.enabled = true;
            }

            assignFromKeyfile(keyFile, "Luminance Curve", "Brightness", labCurve.brightness);
            assignFromKeyfile(keyFile, "Luminance Curve", "Contrast", labCurve.contrast);

            if (ppVersion < 303) {
                // transform Saturation into Chromaticity
                // if Saturation == 0, should we set BWToning on?
                assignFromKeyfile(keyFile, "Luminance Curve", "Saturation", labCurve.chromaticity);
                // transform AvoidColorClipping into AvoidColorShift
                assignFromKeyfile(keyFile, "Luminance Curve", "AvoidColorClipping", labCurve.avoidcolorshift);
            } else {
                if (keyFile.has_key("Luminance Curve", "Chromaticity")) {
                    labCurve.chromaticity = keyFile.get_integer("Luminance Curve", "Chromaticity");

                    if (ppVersion >= 303 && ppVersion < 314 && labCurve.chromaticity == -100) {
                        blackwhite.enabled = true;
                    }
                }

                assignFromKeyfile(keyFile, "Luminance Curve", "AvoidColorShift", labCurve.avoidcolorshift);
                assignFromKeyfile(keyFile, "Luminance Curve", "RedAndSkinTonesProtection", labCurve.rstprotection);
            }

            assignFromKeyfile(keyFile, "Luminance Curve", "LCredsk", labCurve.lcredsk);

            if (ppVersion < 314) {
                // Backward compatibility: If BWtoning is true, Chromaticity has to be set to -100, which will produce the same effect
                // and will enable the b&w toning mode ('a' & 'b' curves)
                if (keyFile.has_key("Luminance Curve", "BWtoning")) {
                    if (keyFile.get_boolean("Luminance Curve", "BWtoning")) {
                        labCurve.chromaticity = -100;
                    }
                }
            }

            assignFromKeyfile(keyFile, "Luminance Curve", "LCurve", labCurve.lcurve);
            assignFromKeyfile(keyFile, "Luminance Curve", "aCurve", labCurve.acurve);
            assignFromKeyfile(keyFile, "Luminance Curve", "bCurve", labCurve.bcurve);
            assignFromKeyfile(keyFile, "Luminance Curve", "ccCurve", labCurve.cccurve);
            assignFromKeyfile(keyFile, "Luminance Curve", "chCurve", labCurve.chcurve);
            assignFromKeyfile(keyFile, "Luminance Curve", "lhCurve", labCurve.lhcurve);
            assignFromKeyfile(keyFile, "Luminance Curve", "hhCurve", labCurve.hhcurve);
            assignFromKeyfile(keyFile, "Luminance Curve", "LcCurve", labCurve.lccurve);
            assignFromKeyfile(keyFile, "Luminance Curve", "ClCurve", labCurve.clcurve);
        }

        if (keyFile.has_group("Sharpening") && RELEVANT_(sharpening)) {
            assignFromKeyfile(keyFile, "Sharpening", "Enabled", sharpening.enabled);
            if (ppVersion >= 334) {
                assignFromKeyfile(keyFile, "Sharpening", "Contrast", sharpening.contrast);
            } else {
                sharpening.contrast = 0;
            }
            assignFromKeyfile(keyFile, "Sharpening", "Radius", sharpening.radius);
            assignFromKeyfile(keyFile, "Sharpening", "BlurRadius", sharpening.blurradius);
            assignFromKeyfile(keyFile, "Sharpening", "Amount", sharpening.amount);

            if (keyFile.has_key("Sharpening", "Threshold")) {
                if (ppVersion < 302) {
                    int thresh = min(keyFile.get_integer("Sharpening", "Threshold"), 2000);
                    sharpening.threshold.setValues(thresh, thresh, 2000, 2000);  // TODO: 2000 is the maximum value and is taken of rtgui/sharpening.cc ; should be changed by the tool modularization
                } else {
                    const std::vector<int> thresh = keyFile.get_integer_list("Sharpening", "Threshold");

                    if (thresh.size() >= 4) {
                        sharpening.threshold.setValues(thresh[0], thresh[1], min(thresh[2], 2000), min(thresh[3], 2000));
                    }
                }
            }

            assignFromKeyfile(keyFile, "Sharpening", "OnlyEdges", sharpening.edgesonly);
            assignFromKeyfile(keyFile, "Sharpening", "EdgedetectionRadius", sharpening.edges_radius);
            assignFromKeyfile(keyFile, "Sharpening", "EdgeTolerance", sharpening.edges_tolerance);
            assignFromKeyfile(keyFile, "Sharpening", "HalocontrolEnabled", sharpening.halocontrol);
            assignFromKeyfile(keyFile, "Sharpening", "HalocontrolAmount", sharpening.halocontrol_amount);
            assignFromKeyfile(keyFile, "Sharpening", "Method", sharpening.method);
            assignFromKeyfile(keyFile, "Sharpening", "DeconvRadius", sharpening.deconvradius);
            assignFromKeyfile(keyFile, "Sharpening", "DeconvAmount", sharpening.deconvamount);
            assignFromKeyfile(keyFile, "Sharpening", "DeconvDamping", sharpening.deconvdamping);
            assignFromKeyfile(keyFile, "Sharpening", "DeconvIterations", sharpening.deconviter);
        }

        if (keyFile.has_group("SharpenMicro") && RELEVANT_(sharpenMicro)) {
            assignFromKeyfile(keyFile, "SharpenMicro", "Enabled", sharpenMicro.enabled);
            assignFromKeyfile(keyFile, "SharpenMicro", "Matrix", sharpenMicro.matrix);
            assignFromKeyfile(keyFile, "SharpenMicro", "Strength", sharpenMicro.amount);
            if (ppVersion >= 334) {
                assignFromKeyfile(keyFile, "SharpenMicro", "Contrast", sharpenMicro.contrast);
            } else {
                sharpenMicro.contrast = 0;
            }
            if (ppVersion >= 346) {
                assignFromKeyfile(keyFile, "SharpenMicro", "Uniformity", sharpenMicro.uniformity);
            } else {
                double temp;
                assignFromKeyfile(keyFile, "SharpenMicro", "Uniformity", temp);
                sharpenMicro.uniformity = temp / 10;
            }
        }

        if (keyFile.has_group("White Balance") && RELEVANT_(wb)) {
            assignFromKeyfile(keyFile, "White Balance", "Enabled", wb.enabled);
            assignFromKeyfile(keyFile, "White Balance", "Setting", wb.method);
            assignFromKeyfile(keyFile, "White Balance", "Temperature", wb.temperature);
            assignFromKeyfile(keyFile, "White Balance", "Green", wb.green);
            assignFromKeyfile(keyFile, "White Balance", "Equal", wb.equal);
            assignFromKeyfile(keyFile, "White Balance", "TemperatureBias", wb.tempBias);
        }

        if (keyFile.has_group("Defringing") && RELEVANT_(defringe)) {
            assignFromKeyfile(keyFile, "Defringing", "Enabled", defringe.enabled);
            assignFromKeyfile(keyFile, "Defringing", "Radius", defringe.radius);

            if (keyFile.has_key("Defringing", "Threshold")) {
                defringe.threshold = (float)keyFile.get_integer("Defringing", "Threshold");
            }

            if (ppVersion < 310) {
                defringe.threshold = sqrt(defringe.threshold * 33.f / 5.f);
            }

            assignFromKeyfile(keyFile, "Defringing", "HueCurve", defringe.huecurve);
        }

        if (keyFile.has_group("Impulse Denoising") && RELEVANT_(impulseDenoise)) {
            assignFromKeyfile(keyFile, "Impulse Denoising", "Enabled", impulseDenoise.enabled);
            assignFromKeyfile(keyFile, "Impulse Denoising", "Threshold", impulseDenoise.thresh);
        }

        if (ppVersion < 346) {
            if (keyFile.has_group("Directional Pyramid Denoising") && RELEVANT_(denoise)) { //TODO: No longer an accurate description for FT denoise
                assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "Enabled", denoise.enabled);
                if (assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "Median", denoise.smoothingEnabled)) {
                    denoise.smoothingMethod = DenoiseParams::SmoothingMethod::MEDIAN;
                }
                assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "Luma", denoise.luminance);
                assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "Ldetail", denoise.luminanceDetail);
                assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "Chroma", denoise.chrominance);
                Glib::ustring val;
                if (assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "Method", val)) {
                    if (val == "RGB") {
                        denoise.colorSpace = DenoiseParams::ColorSpace::RGB;
                    } else {
                        denoise.colorSpace = DenoiseParams::ColorSpace::LAB;
                    }
                }
                if (assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "C2Method", val)) {
                    if (val == "MANU") {
                        denoise.chrominanceMethod = DenoiseParams::ChrominanceMethod::MANUAL;
                    } else {
                        denoise.chrominanceMethod = DenoiseParams::ChrominanceMethod::AUTOMATIC;
                    }
                }
                if (assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "SMethod", val)) {
                    denoise.aggressive = (val == "shalbi");
                }
                if (assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "MedMethod", val)) {
                    const std::vector<Glib::ustring> medtps = { "soft", "33", "55soft", "55", "77", "99" };
                    auto it = std::find(medtps.begin(), medtps.end(), val);
                    if (it != medtps.end()) {
                        denoise.medianType = static_cast<DenoiseParams::MedianType>(it - medtps.begin());
                    }
                }
                if (assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "MethodMed", val)) {
                    const std::vector<Glib::ustring> med = { "Lonly", "ab", "Lpab", "Lab", "RGB" };
                    auto it = std::find(med.begin(), med.end(), val);
                    if (it != med.end()) {
                        denoise.medianMethod = static_cast<DenoiseParams::MedianMethod>(it - med.begin());
                    }
                }
                assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "Redchro", denoise.chrominanceRedGreen);
                assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "Bluechro", denoise.chrominanceBlueYellow);
                assignFromKeyfile(keyFile, "Directional Pyramid Denoising", "Passes", denoise.medianIterations);
            }
        } else {
            if (keyFile.has_group("Denoise") && RELEVANT_(denoise)) {
                assignFromKeyfile(keyFile, "Denoise", "Enabled", denoise.enabled);
                int val;
                Glib::ustring sval;
                if (assignFromKeyfile(keyFile, "Denoise", "ColorSpace", sval)) {
                    denoise.colorSpace = (sval == "RGB" ? DenoiseParams::ColorSpace::RGB : DenoiseParams::ColorSpace::LAB);
                }
                assignFromKeyfile(keyFile, "Denoise", "Aggressive", denoise.aggressive);
                assignFromKeyfile(keyFile, "Denoise", "Gamma", denoise.gamma);
                assignFromKeyfile(keyFile, "Denoise", "Luminance", denoise.luminance);
                assignFromKeyfile(keyFile, "Denoise", "LuminanceDetail", denoise.luminanceDetail);
                if (assignFromKeyfile(keyFile, "Denoise", "ChrominanceMethod", val)) {
                    denoise.chrominanceMethod = static_cast<DenoiseParams::ChrominanceMethod>(val);
                }
                assignFromKeyfile(keyFile, "Denoise", "ChrominanceAutoFactor", denoise.chrominanceAutoFactor);
                assignFromKeyfile(keyFile, "Denoise", "Chrominance", denoise.chrominance);
                assignFromKeyfile(keyFile, "Denoise", "ChrominanceRedGreen", denoise.chrominanceRedGreen);
                assignFromKeyfile(keyFile, "Denoise", "ChrominanceBlueYellow", denoise.chrominanceBlueYellow);
                assignFromKeyfile(keyFile, "Denoise", "SmoothingEnabled", denoise.smoothingEnabled);
                if (assignFromKeyfile(keyFile, "Denoise", "SmoothingMethod", val)) {
                    denoise.smoothingMethod = static_cast<DenoiseParams::SmoothingMethod>(val);
                }
                if (assignFromKeyfile(keyFile, "Denoise", "MedianType", val)) {
                    denoise.medianType = static_cast<DenoiseParams::MedianType>(val);
                }
                if (assignFromKeyfile(keyFile, "Denoise", "MedianMethod", val)) {
                    denoise.medianMethod = static_cast<DenoiseParams::MedianMethod>(val);
                }
                assignFromKeyfile(keyFile, "Denoise", "MedianIterations", denoise.medianIterations);
                assignFromKeyfile(keyFile, "Denoise", "GuidedLumaRadius", denoise.guidedLumaRadius);
                assignFromKeyfile(keyFile, "Denoise", "GuidedChromaRadius", denoise.guidedChromaRadius);
                assignFromKeyfile(keyFile, "Denoise", "GuidedLumaStrength", denoise.guidedLumaStrength);
                assignFromKeyfile(keyFile, "Denoise", "GuidedChromaStrength", denoise.guidedChromaStrength);
            }
        }            

        if (keyFile.has_group("EPD") && RELEVANT_(epd)) {
            assignFromKeyfile(keyFile, "EPD", "Enabled", epd.enabled);
                
            std::vector<EPDParams::Region> ll;
            std::vector<LabCorrectionMask> lm;
            bool found = false;
            bool done = false;
            for (int i = 0; !done; ++i) {
                EPDParams::Region cur;
                LabCorrectionMask curmask;
                done = true;
                std::string n = i ? std::string("_") + std::to_string(i) : std::string("");
                if (assignFromKeyfile(keyFile, "EPD", Glib::ustring("Strength") + n, cur.strength)) {
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, "EPD", Glib::ustring("Gamma") + n, cur.gamma)) {
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, "EPD", Glib::ustring("EdgeStopping") + n, cur.edgeStopping)) {
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, "EPD", Glib::ustring("Scale") + n, cur.scale)) {
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, "EPD", Glib::ustring("ReweightingIterates") + n, cur.reweightingIterates)) {
                    found = true;
                    done = false;
                }
                if (curmask.load(keyFile, "EPD", "", n)) {
                    found = true;
                    done = false;
                }
                if (!done) {
                    ll.emplace_back(cur);
                    lm.emplace_back(curmask);
                }
            }
            if (found) {
                epd.regions = std::move(ll);
                epd.labmasks = std::move(lm);
            }
            assert(epd.regions.size() == epd.labmasks.size());
            assignFromKeyfile(keyFile, "EPD", "ShowMask", epd.showMask);
        }

        if (keyFile.has_group("FattalToneMapping") && RELEVANT_(fattal)) {
            assignFromKeyfile(keyFile, "FattalToneMapping", "Enabled", fattal.enabled);
            assignFromKeyfile(keyFile, "FattalToneMapping", "Threshold", fattal.threshold);
            assignFromKeyfile(keyFile, "FattalToneMapping", "Amount", fattal.amount);
            assignFromKeyfile(keyFile, "FattalToneMapping", "Anchor", fattal.anchor);
        }

        if (keyFile.has_group("LogEncoding") && RELEVANT_(logenc)) {
            assignFromKeyfile(keyFile, "LogEncoding", "Enabled", logenc.enabled);
            assignFromKeyfile(keyFile, "LogEncoding", "Auto", logenc.autocompute);
            assignFromKeyfile(keyFile, "LogEncoding", "AutoGray", logenc.autogray);
            if (ppVersion < 349) {
                assignFromKeyfile(keyFile, "LogEncoding", "GrayPoint", logenc.sourceGray);
            } else {
                assignFromKeyfile(keyFile, "LogEncoding", "SourceGray", logenc.sourceGray);
                assignFromKeyfile(keyFile, "LogEncoding", "TargetGray", logenc.targetGray);
            }
            assignFromKeyfile(keyFile, "LogEncoding", "BlackEv", logenc.blackEv);
            assignFromKeyfile(keyFile, "LogEncoding", "WhiteEv", logenc.whiteEv);
            assignFromKeyfile(keyFile, "LogEncoding", "Detail", logenc.detail);
        }

        if (keyFile.has_group ("Shadows & Highlights") && ppVersion >= 333 && RELEVANT_(sh)) {
            assignFromKeyfile(keyFile, "Shadows & Highlights", "Enabled", sh.enabled);
            assignFromKeyfile(keyFile, "Shadows & Highlights", "Highlights", sh.highlights);
            assignFromKeyfile(keyFile, "Shadows & Highlights", "HighlightTonalWidth", sh.htonalwidth);
            assignFromKeyfile(keyFile, "Shadows & Highlights", "Shadows", sh.shadows);
            assignFromKeyfile(keyFile, "Shadows & Highlights", "ShadowTonalWidth", sh.stonalwidth);
            assignFromKeyfile(keyFile, "Shadows & Highlights", "Radius", sh.radius);
            if (ppVersion >= 344) {
                assignFromKeyfile(keyFile, "Shadows & Highlights", "Lab", sh.lab);
            } else {
                sh.lab = true;
            }
        }

        if (keyFile.has_group("ToneEqualizer") && RELEVANT_(toneEqualizer)) {
            assignFromKeyfile(keyFile, "ToneEqualizer", "Enabled", toneEqualizer.enabled);
            for (size_t i = 0; i < toneEqualizer.bands.size(); ++i) {
                assignFromKeyfile(keyFile, "ToneEqualizer", "Band" + std::to_string(i), toneEqualizer.bands[i]);
            }
            assignFromKeyfile(keyFile, "ToneEqualizer", "Detail", toneEqualizer.detail);
        }
        
        if (keyFile.has_group("Crop") && RELEVANT_(crop)) {
            assignFromKeyfile(keyFile, "Crop", "Enabled", crop.enabled);
            assignFromKeyfile(keyFile, "Crop", "X", crop.x);
            assignFromKeyfile(keyFile, "Crop", "Y", crop.y);

            if (keyFile.has_key("Crop", "W")) {
                crop.w = std::max(keyFile.get_integer("Crop", "W"), 1);
            }

            if (keyFile.has_key("Crop", "H")) {
                crop.h = std::max(keyFile.get_integer("Crop", "H"), 1);
            }

            assignFromKeyfile(keyFile, "Crop", "FixedRatio", crop.fixratio);

            if (assignFromKeyfile(keyFile, "Crop", "Ratio", crop.ratio)) {
                //backwards compatibility for crop.ratio
                if (crop.ratio == "DIN") {
                    crop.ratio = "1.414 - DIN EN ISO 216";
                }

                if (crop.ratio == "8.5:11") {
                    crop.ratio = "8.5:11 - US Letter";
                }

                if (crop.ratio == "11:17") {
                    crop.ratio = "11:17 - Tabloid";
                }
            }

            assignFromKeyfile(keyFile, "Crop", "Orientation", crop.orientation);
            assignFromKeyfile(keyFile, "Crop", "Guide", crop.guide);
        }

        if (keyFile.has_group("Coarse Transformation") && RELEVANT_(coarse)) {
            assignFromKeyfile(keyFile, "Coarse Transformation", "Rotate", coarse.rotate);
            assignFromKeyfile(keyFile, "Coarse Transformation", "HorizontalFlip", coarse.hflip);
            assignFromKeyfile(keyFile, "Coarse Transformation", "VerticalFlip", coarse.vflip);
        }

        if (keyFile.has_group("Rotation") && RELEVANT_(rotate)) {
            assignFromKeyfile(keyFile, "Rotation", "Enabled", rotate.enabled);
            assignFromKeyfile(keyFile, "Rotation", "Degree", rotate.degree);
        }

        if (keyFile.has_group("Common Properties for Transformations") && RELEVANT_(commonTrans)) {
            assignFromKeyfile(keyFile, "Common Properties for Transformations", "AutoFill", commonTrans.autofill);
        }

        if (keyFile.has_group("Distortion") && RELEVANT_(distortion)) {
            assignFromKeyfile(keyFile, "Distortion", "Enabled", distortion.enabled);
            assignFromKeyfile(keyFile, "Distortion", "Amount", distortion.amount);
        }

        if (keyFile.has_group("LensProfile") && RELEVANT_(lensProf)) {
            if (keyFile.has_key("LensProfile", "LcMode")) {
                lensProf.lcMode = lensProf.getMethodNumber(keyFile.get_string("LensProfile", "LcMode"));
            }

            if (keyFile.has_key("LensProfile", "LCPFile")) {
                lensProf.lcpFile = expandRelativePath(fname, "", keyFile.get_string("LensProfile", "LCPFile"));

                if (ppVersion < 327 && !lensProf.lcpFile.empty()) {
                    lensProf.lcMode = LensProfParams::LcMode::LCP;
                }
            }

            assignFromKeyfile(keyFile, "LensProfile", "UseDistortion", lensProf.useDist);
            assignFromKeyfile(keyFile, "LensProfile", "UseVignette", lensProf.useVign);
            assignFromKeyfile(keyFile, "LensProfile", "UseCA", lensProf.useCA);

            if (keyFile.has_key("LensProfile", "LFCameraMake")) {
                lensProf.lfCameraMake = keyFile.get_string("LensProfile", "LFCameraMake");
            }

            if (keyFile.has_key("LensProfile", "LFCameraModel")) {
                lensProf.lfCameraModel = keyFile.get_string("LensProfile", "LFCameraModel");
            }

            if (keyFile.has_key("LensProfile", "LFLens")) {
                lensProf.lfLens = keyFile.get_string("LensProfile", "LFLens");
            }
        }

        if (keyFile.has_group("Perspective") && RELEVANT_(perspective)) {
            assignFromKeyfile(keyFile, "Perspective", "Enabled", perspective.enabled);
            assignFromKeyfile(keyFile, "Perspective", "Horizontal", perspective.horizontal);
            assignFromKeyfile(keyFile, "Perspective", "Vertical", perspective.vertical);
            assignFromKeyfile(keyFile, "Perspective", "Angle", perspective.angle);
            assignFromKeyfile(keyFile, "Perspective", "Shear", perspective.shear);
            assignFromKeyfile(keyFile, "Perspective", "FocalLength", perspective.flength);
            assignFromKeyfile(keyFile, "Perspective", "CropFactor", perspective.cropfactor);
            assignFromKeyfile(keyFile, "Perspective", "Aspect", perspective.aspect);
        }

        if (keyFile.has_group("Gradient") && RELEVANT_(gradient)) {
            assignFromKeyfile(keyFile, "Gradient", "Enabled", gradient.enabled);
            assignFromKeyfile(keyFile, "Gradient", "Degree", gradient.degree);
            assignFromKeyfile(keyFile, "Gradient", "Feather", gradient.feather);
            assignFromKeyfile(keyFile, "Gradient", "Strength", gradient.strength);
            assignFromKeyfile(keyFile, "Gradient", "CenterX", gradient.centerX);
            assignFromKeyfile(keyFile, "Gradient", "CenterY", gradient.centerY);
        }

        if (keyFile.has_group("PCVignette") && RELEVANT_(pcvignette)) {
            assignFromKeyfile(keyFile, "PCVignette", "Enabled", pcvignette.enabled);
            assignFromKeyfile(keyFile, "PCVignette", "Strength", pcvignette.strength);
            assignFromKeyfile(keyFile, "PCVignette", "Feather", pcvignette.feather);
            assignFromKeyfile(keyFile, "PCVignette", "Roundness", pcvignette.roundness);
        }

        if (keyFile.has_group("CACorrection") && RELEVANT_(cacorrection)) {
            assignFromKeyfile(keyFile, "CACorrection", "Enabled", cacorrection.enabled);
            assignFromKeyfile(keyFile, "CACorrection", "Red", cacorrection.red);
            assignFromKeyfile(keyFile, "CACorrection", "Blue", cacorrection.blue);
        }

        if (keyFile.has_group("Vignetting Correction") && RELEVANT_(vignetting)) {
            assignFromKeyfile(keyFile, "Vignetting Correction", "Enabled", vignetting.enabled);
            assignFromKeyfile(keyFile, "Vignetting Correction", "Amount", vignetting.amount);
            assignFromKeyfile(keyFile, "Vignetting Correction", "Radius", vignetting.radius);
            assignFromKeyfile(keyFile, "Vignetting Correction", "Strength", vignetting.strength);
            assignFromKeyfile(keyFile, "Vignetting Correction", "CenterX", vignetting.centerX);
            assignFromKeyfile(keyFile, "Vignetting Correction", "CenterY", vignetting.centerY);
        }

        if (keyFile.has_group("Resize") && RELEVANT_(resize)) {
            assignFromKeyfile(keyFile, "Resize", "Enabled", resize.enabled);
            assignFromKeyfile(keyFile, "Resize", "Scale", resize.scale);
            assignFromKeyfile(keyFile, "Resize", "AppliesTo", resize.appliesTo);
            assignFromKeyfile(keyFile, "Resize", "Method", resize.method);
            assignFromKeyfile(keyFile, "Resize", "DataSpecified", resize.dataspec);
            assignFromKeyfile(keyFile, "Resize", "Width", resize.width);
            assignFromKeyfile(keyFile, "Resize", "Height", resize.height);
            if (ppVersion >= 339) {
                assignFromKeyfile(keyFile, "Resize", "AllowUpscaling", resize.allowUpscaling);
            } else {
                resize.allowUpscaling = false;
            }
        }

        if (keyFile.has_group("PostResizeSharpening") && RELEVANT_(prsharpening)) {
            assignFromKeyfile(keyFile, "PostResizeSharpening", "Enabled", prsharpening.enabled);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "Contrast", prsharpening.contrast);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "Radius", prsharpening.radius);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "Amount", prsharpening.amount);

            if (keyFile.has_key("PostResizeSharpening", "Threshold")) {
                if (ppVersion < 302) {
                    int thresh = min(keyFile.get_integer("PostResizeSharpening", "Threshold"), 2000);
                    prsharpening.threshold.setValues(thresh, thresh, 2000, 2000);  // TODO: 2000 is the maximum value and is taken of rtgui/sharpening.cc ; should be changed by the tool modularization
                } else {
                    const std::vector<int> thresh = keyFile.get_integer_list("PostResizeSharpening", "Threshold");

                    if (thresh.size() >= 4) {
                        prsharpening.threshold.setValues(thresh[0], thresh[1], min(thresh[2], 2000), min(thresh[3], 2000));
                    }
                }
            }

            assignFromKeyfile(keyFile, "PostResizeSharpening", "OnlyEdges", prsharpening.edgesonly);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "EdgedetectionRadius", prsharpening.edges_radius);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "EdgeTolerance", prsharpening.edges_tolerance);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "HalocontrolEnabled", prsharpening.halocontrol);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "HalocontrolAmount", prsharpening.halocontrol_amount);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "Method", prsharpening.method);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "DeconvRadius", prsharpening.deconvradius);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "DeconvAmount", prsharpening.deconvamount);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "DeconvDamping", prsharpening.deconvdamping);
            assignFromKeyfile(keyFile, "PostResizeSharpening", "DeconvIterations", prsharpening.deconviter);
        }

        if (keyFile.has_group("Color Management") && RELEVANT_(icm)) {
            if (keyFile.has_key("Color Management", "InputProfile")) {
                icm.inputProfile = expandRelativePath(fname, "file:", keyFile.get_string("Color Management", "InputProfile"));
            }

            assignFromKeyfile(keyFile, "Color Management", "ToneCurve", icm.toneCurve);
            assignFromKeyfile(keyFile, "Color Management", "ApplyLookTable", icm.applyLookTable);
            assignFromKeyfile(keyFile, "Color Management", "ApplyBaselineExposureOffset", icm.applyBaselineExposureOffset);
            assignFromKeyfile(keyFile, "Color Management", "ApplyHueSatMap", icm.applyHueSatMap);
            assignFromKeyfile(keyFile, "Color Management", "DCPIlluminant", icm.dcpIlluminant);
            assignFromKeyfile(keyFile, "Color Management", "WorkingProfile", icm.workingProfile);

            assignFromKeyfile(keyFile, "Color Management", "OutputProfile", icm.outputProfile);
            if (ppVersion < 341) {
                if (icm.outputProfile == "RT_Medium_gsRGB") {
                    icm.outputProfile = "RTv4_Medium";
                } else if (icm.outputProfile == "RT_Large_gBT709" || icm.outputProfile == "RT_Large_g10" || icm.outputProfile == "RT_Large_gsRGB") {
                    icm.outputProfile = "RTv4_Large";
                } else if (icm.outputProfile == "WideGamutRGB") {
                    icm.outputProfile = "RTv4_Wide";
                } else if (icm.outputProfile == "RT_sRGB_gBT709" || icm.outputProfile == "RT_sRGB_g10" || icm.outputProfile == "RT_sRGB") {
                    icm.outputProfile = "RTv4_sRGB";
                } else if (icm.outputProfile == "BetaRGB") { // Have we ever provided this profile ? Should we convert this filename ?
                    icm.outputProfile = "RTv4_Beta";
                } else if (icm.outputProfile == "BestRGB") { // Have we ever provided this profile ? Should we convert this filename ?
                    icm.outputProfile = "RTv4_Best";
                } else if (icm.outputProfile == "Rec2020") {
                    icm.outputProfile = "RTv4_Rec2020";
                } else if (icm.outputProfile == "Bruce") { // Have we ever provided this profile ? Should we convert this filename ?
                    icm.outputProfile = "RTv4_Bruce";
                } else if (icm.outputProfile == "ACES") {
                    icm.outputProfile = "RTv4_ACES-AP0";
                }
            }
            if (keyFile.has_key("Color Management", "OutputProfileIntent")) {
                Glib::ustring intent = keyFile.get_string("Color Management", "OutputProfileIntent");

                if (intent == "Perceptual") {
                    icm.outputIntent = RI_PERCEPTUAL;
                } else if (intent == "Relative") {
                    icm.outputIntent = RI_RELATIVE;
                } else if (intent == "Saturation") {
                    icm.outputIntent = RI_SATURATION;
                } else if (intent == "Absolute") {
                    icm.outputIntent = RI_ABSOLUTE;
                }
            }
            assignFromKeyfile(keyFile, "Color Management", "OutputBPC", icm.outputBPC);
        }

        if (keyFile.has_group("Directional Pyramid Equalizer") && RELEVANT_(dirpyrequalizer)) {
            if (ppVersion < 347) {
                assignFromKeyfile(keyFile, "Directional Pyramid Equalizer", "Enabled", dirpyrequalizer.enabled);

                dirpyrequalizer.levels = {DirPyrEqualizerParams::Levels()};
                dirpyrequalizer.labmasks = {LabCorrectionMask()};

                auto &l = dirpyrequalizer.levels[0];

                if (ppVersion < 316) {
                    for (int i = 0; i < 5; i ++) {
                        std::stringstream ss;
                        ss << "Mult" << i;

                        if (keyFile.has_key("Directional Pyramid Equalizer", ss.str())) {
                            if (i == 4) {
                                l.threshold = keyFile.get_double("Directional Pyramid Equalizer", ss.str());
                            } else {
                                l.mult[i] = keyFile.get_double("Directional Pyramid Equalizer", ss.str());
                            }
                        }
                    }

                    l.mult[4] = 1.0;
                } else {
                    // 5 level wavelet + dedicated threshold parameter
                    for (int i = 0; i < 6; i ++) {
                        std::stringstream ss;
                        ss << "Mult" << i;

                        if (keyFile.has_key("Directional Pyramid Equalizer", ss.str())) {
                            l.mult[i] = keyFile.get_double("Directional Pyramid Equalizer", ss.str());
                        }
                    }

                    assignFromKeyfile(keyFile, "Directional Pyramid Equalizer", "Threshold", l.threshold);
//                    assignFromKeyfile(keyFile, "Directional Pyramid Equalizer", "Skinprotect", dirpyrequalizer.skinprotect);
                    // TODO - port skin protection from old pp3's
                }
            } else {
                assignFromKeyfile(keyFile, "Directional Pyramid Equalizer", "Enabled", dirpyrequalizer.enabled);
                
                std::vector<DirPyrEqualizerParams::Levels> ll;
                std::vector<LabCorrectionMask> lm;
                bool found = false;
                bool done = false;
                for (int i = 1; !done; ++i) {
                    DirPyrEqualizerParams::Levels cur;
                    LabCorrectionMask curmask;
                    done = true;
                    std::string n = std::to_string(i);
                    for (int j = 0; j < 6; ++j) {
                        if (assignFromKeyfile(keyFile, "Directional Pyramid Equalizer", Glib::ustring("Mult") + std::to_string(j) + "_" + n, cur.mult[j])) {
                            found = true;
                            done = false;
                        }
                    }
                    if (assignFromKeyfile(keyFile, "Directional Pyramid Equalizer", Glib::ustring("Threshold_") + n, cur.threshold)) {
                        found = true;
                        done = false;
                    }
                    if (curmask.load(keyFile, "Directional Pyramid Equalizer", "", Glib::ustring("_") + n)) {
                        found = true;
                        done = false;
                    }
                    if (!done) {
                        ll.emplace_back(cur);
                        lm.emplace_back(curmask);
                    }
                }
                if (found) {
                    dirpyrequalizer.levels = std::move(ll);
                    dirpyrequalizer.labmasks = std::move(lm);
                }
                assert(dirpyrequalizer.levels.size() == dirpyrequalizer.labmasks.size());
                assignFromKeyfile(keyFile, "Directional Pyramid Equalizer", "ShowMask", dirpyrequalizer.showMask);
            }
        }

        if (keyFile.has_group("SoftLight") && RELEVANT_(softlight)) {
            assignFromKeyfile(keyFile, "SoftLight", "Enabled", softlight.enabled);
            assignFromKeyfile(keyFile, "SoftLight", "Strength", softlight.strength);
        }

        if (keyFile.has_group("Dehaze") && RELEVANT_(dehaze)) {
            assignFromKeyfile(keyFile, "Dehaze", "Enabled", dehaze.enabled);
            assignFromKeyfile(keyFile, "Dehaze", "Strength", dehaze.strength);
            assignFromKeyfile(keyFile, "Dehaze", "ShowDepthMap", dehaze.showDepthMap);
            assignFromKeyfile(keyFile, "Dehaze", "Depth", dehaze.depth);
            assignFromKeyfile(keyFile, "Dehaze", "Luminance", dehaze.luminance);
        }
        
        if (keyFile.has_group("Film Simulation") && RELEVANT_(filmSimulation)) {
            assignFromKeyfile(keyFile, "Film Simulation", "Enabled", filmSimulation.enabled);
            assignFromKeyfile(keyFile, "Film Simulation", "ClutFilename", filmSimulation.clutFilename);

            if (keyFile.has_key("Film Simulation", "Strength")) {
                if (ppVersion < 321) {
                    filmSimulation.strength = keyFile.get_double("Film Simulation", "Strength") * 100 + 0.1;
                } else {
                    filmSimulation.strength = keyFile.get_integer("Film Simulation", "Strength");
                }
            }
        }

        if (keyFile.has_group("RGB Curves") && RELEVANT_(rgbCurves)) {
            if (ppVersion >= 329) {
                assignFromKeyfile(keyFile, "RGB Curves", "Enabled", rgbCurves.enabled);
            } else {
                rgbCurves.enabled = true;
            }

            assignFromKeyfile(keyFile, "RGB Curves", "LumaMode", rgbCurves.lumamode);
            assignFromKeyfile(keyFile, "RGB Curves", "rCurve", rgbCurves.rcurve);
            assignFromKeyfile(keyFile, "RGB Curves", "gCurve", rgbCurves.gcurve);
            assignFromKeyfile(keyFile, "RGB Curves", "bCurve", rgbCurves.bcurve);
        }

        if (keyFile.has_group("Grain") && RELEVANT_(grain)) {
            assignFromKeyfile(keyFile, "Grain", "Enabled", grain.enabled);
            assignFromKeyfile(keyFile, "Grain", "ISO", grain.iso);
            assignFromKeyfile(keyFile, "Grain", "Strength", grain.strength);
            assignFromKeyfile(keyFile, "Grain", "Scale", grain.scale);
        }

        if (keyFile.has_group("GuidedSmoothing") && RELEVANT_(smoothing)) {
            assignFromKeyfile(keyFile, "GuidedSmoothing", "Enabled", smoothing.enabled);
                
            std::vector<GuidedSmoothingParams::Region> ll;
            std::vector<LabCorrectionMask> lm;
            bool found = false;
            bool done = false;
            for (int i = 1; !done; ++i) {
                GuidedSmoothingParams::Region cur;
                LabCorrectionMask curmask;
                done = true;
                std::string n = std::to_string(i);
                int c;
                if (assignFromKeyfile(keyFile, "GuidedSmoothing", Glib::ustring("Channel_") + n, c)) {
                    cur.channel = GuidedSmoothingParams::Region::Channel(c);
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, "GuidedSmoothing", Glib::ustring("Radius_") + n, cur.radius)) {
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, "GuidedSmoothing", Glib::ustring("Epsilon_") + n, cur.epsilon)) {
                    found = true;
                    done = false;
                }
                if (curmask.load(keyFile, "GuidedSmoothing", "", Glib::ustring("_") + n)) {
                    found = true;
                    done = false;
                }
                if (!done) {
                    ll.emplace_back(cur);
                    lm.emplace_back(curmask);
                }
            }
            if (found) {
                smoothing.regions = std::move(ll);
                smoothing.labmasks = std::move(lm);
            }
            assert(smoothing.regions.size() == smoothing.labmasks.size());
            assignFromKeyfile(keyFile, "GuidedSmoothing", "ShowMask", smoothing.showMask);
        }

        const char *ccgroup = "ColorCorrection";
        if (keyFile.has_group(ccgroup) && RELEVANT_(colorcorrection)) {
            const Glib::ustring prefix = "";
            assignFromKeyfile(keyFile, ccgroup, "Enabled", colorcorrection.enabled);
            std::vector<ColorCorrectionParams::LabCorrectionRegion> lg;
            std::vector<LabCorrectionMask> lm;
            bool found = false;
            bool done = false;
            for (int i = 1; !done; ++i) {
                ColorCorrectionParams::LabCorrectionRegion cur;
                LabCorrectionMask curmask;
                done = true;
                std::string n = std::to_string(i);
                if (assignFromKeyfile(keyFile, ccgroup, prefix + Glib::ustring("A_") + n, cur.a)) {
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, ccgroup, prefix + Glib::ustring("B_") + n, cur.b)) {
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, ccgroup, prefix + Glib::ustring("Saturation_") + n, cur.saturation)) {
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, ccgroup, prefix + Glib::ustring("Slope_") + n, cur.slope)) {
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, ccgroup, prefix + Glib::ustring("Offset_") + n, cur.offset)) {
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, ccgroup, prefix + Glib::ustring("Power_") + n, cur.power)) {
                    found = true;
                    done = false;
                }
                if (assignFromKeyfile(keyFile, ccgroup, prefix + Glib::ustring("Channel_") + n, cur.channel)) {
                    found = true;
                    done = false;
                }
                if (curmask.load(keyFile, ccgroup, prefix, Glib::ustring("_") + n)) {
                    found = true;
                    done = false;
                }
                if (!done) {
                    lg.emplace_back(cur);
                    lm.emplace_back(curmask);
                }
            }
            if (found) {
                colorcorrection.regions = std::move(lg);
                colorcorrection.labmasks = std::move(lm);
            }
            assert(colorcorrection.regions.size() == colorcorrection.labmasks.size());
            assignFromKeyfile(keyFile, ccgroup, ppVersion < 348 ? "showMask" : "LabRegionsShowMask", colorcorrection.showMask);
        }

        if (keyFile.has_group("RAW")) {
            if (RELEVANT_(darkframe)) {
                assignFromKeyfile(keyFile, "RAW", "DarkFrameEnabled", raw.enable_darkframe);
                if (keyFile.has_key("RAW", "DarkFrame")) {
                    raw.dark_frame = expandRelativePath(fname, "", keyFile.get_string("RAW", "DarkFrame"));
                }

                assignFromKeyfile(keyFile, "RAW", "DarkFrameAuto", raw.df_autoselect);
            }

            if (RELEVANT_(flatfield)) {
                assignFromKeyfile(keyFile, "RAW", "FlatFieldEnabled", raw.enable_flatfield);
                if (keyFile.has_key("RAW", "FlatFieldFile")) {
                    raw.ff_file = expandRelativePath(fname, "", keyFile.get_string("RAW", "FlatFieldFile"));
                }

                assignFromKeyfile(keyFile, "RAW", "FlatFieldAutoSelect", raw.ff_AutoSelect);
                assignFromKeyfile(keyFile, "RAW", "FlatFieldBlurRadius", raw.ff_BlurRadius);
                assignFromKeyfile(keyFile, "RAW", "FlatFieldBlurType", raw.ff_BlurType);
                assignFromKeyfile(keyFile, "RAW", "FlatFieldAutoClipControl", raw.ff_AutoClipControl);

                if (ppVersion < 328) {
                    // With ppversion < 328 this value was stored as a boolean, which is nonsense.
                    // To avoid annoying warnings we skip reading and assume 0.
                    raw.ff_clipControl = 0;
                } else {
                    assignFromKeyfile(keyFile, "RAW", "FlatFieldClipControl", raw.ff_clipControl);
                }
            }

            if (RELEVANT_(rawCA)) {
                assignFromKeyfile(keyFile, "RAW", "CAEnabled", raw.enable_ca);
                assignFromKeyfile(keyFile, "RAW", "CA", raw.ca_autocorrect);
                if (ppVersion >= 342) {
                    assignFromKeyfile(keyFile, "RAW", "CAAutoIterations", raw.caautoiterations);
                } else {
                    raw.caautoiterations = 1;
                }

                if (ppVersion >= 343) {
                    assignFromKeyfile(keyFile, "RAW", "CAAvoidColourshift", raw.ca_avoidcolourshift);
                } else {
                    raw.ca_avoidcolourshift = false;
                }
                assignFromKeyfile(keyFile, "RAW", "CARed", raw.cared);
                assignFromKeyfile(keyFile, "RAW", "CABlue", raw.cablue);
            }

            if (RELEVANT_(hotDeadPixelFilter)) {
                assignFromKeyfile(keyFile, "RAW", "HotDeadPixelEnabled", raw.enable_hotdeadpix);
                // For compatibility to elder pp3 versions
                assignFromKeyfile(keyFile, "RAW", "HotDeadPixels", raw.hotPixelFilter);
                raw.deadPixelFilter = raw.hotPixelFilter;

                assignFromKeyfile(keyFile, "RAW", "HotPixelFilter", raw.hotPixelFilter);
                assignFromKeyfile(keyFile, "RAW", "DeadPixelFilter", raw.deadPixelFilter);
                assignFromKeyfile(keyFile, "RAW", "HotDeadPixelThresh", raw.hotdeadpix_thresh);
            }
            if (RELEVANT_(rawWhite)) {
                assignFromKeyfile(keyFile, "RAW", "PreExposureEnabled", raw.enable_whitepoint);
                assignFromKeyfile(keyFile, "RAW", "PreExposure", raw.expos);
            }

            if (ppVersion < 320) {
                if (RELEVANT_(demosaic)) {
                    assignFromKeyfile(keyFile, "RAW", "Method", raw.bayersensor.method);
                    assignFromKeyfile(keyFile, "RAW", "CcSteps", raw.bayersensor.ccSteps);
                    assignFromKeyfile(keyFile, "RAW", "DCBIterations", raw.bayersensor.dcb_iterations);
                    assignFromKeyfile(keyFile, "RAW", "DCBEnhance", raw.bayersensor.dcb_enhance);
                    assignFromKeyfile(keyFile, "RAW", "LMMSEIterations", raw.bayersensor.lmmse_iterations);
                }
                if (RELEVANT_(rawPreprocessing)) {
                    assignFromKeyfile(keyFile, "RAW", "LineDenoise", raw.bayersensor.linenoise);
                    assignFromKeyfile(keyFile, "RAW", "GreenEqThreshold", raw.bayersensor.greenthresh);
                }
                if (RELEVANT_(rawBlack)) {
                    assignFromKeyfile(keyFile, "RAW", "PreBlackzero", raw.bayersensor.black0);
                    assignFromKeyfile(keyFile, "RAW", "PreBlackone", raw.bayersensor.black1);
                    assignFromKeyfile(keyFile, "RAW", "PreBlacktwo", raw.bayersensor.black2);
                    assignFromKeyfile(keyFile, "RAW", "PreBlackthree", raw.bayersensor.black3);
                    assignFromKeyfile(keyFile, "RAW", "PreTwoGreen", raw.bayersensor.twogreen);
                }
            }
        }

        if (keyFile.has_group("RAW Bayer")) {
            if (RELEVANT_(demosaic)) {
                assignFromKeyfile(keyFile, "RAW Bayer", "Method", raw.bayersensor.method);
                assignFromKeyfile(keyFile, "RAW Bayer", "Border", raw.bayersensor.border);

                if (keyFile.has_key("RAW Bayer", "ImageNum")) {
                    raw.bayersensor.imageNum = keyFile.get_integer("RAW Bayer", "ImageNum") - 1;
                }

                assignFromKeyfile(keyFile, "RAW Bayer", "CcSteps", raw.bayersensor.ccSteps);
            }

            if (RELEVANT_(rawBlack)) {
                assignFromKeyfile(keyFile, "RAW Bayer", "PreBlackEnabled", raw.bayersensor.enable_black);
                assignFromKeyfile(keyFile, "RAW Bayer", "PreBlack0", raw.bayersensor.black0);
                assignFromKeyfile(keyFile, "RAW Bayer", "PreBlack1", raw.bayersensor.black1);
                assignFromKeyfile(keyFile, "RAW Bayer", "PreBlack2", raw.bayersensor.black2);
                assignFromKeyfile(keyFile, "RAW Bayer", "PreBlack3", raw.bayersensor.black3);
                assignFromKeyfile(keyFile, "RAW Bayer", "PreTwoGreen", raw.bayersensor.twogreen);
            }

            if (RELEVANT_(rawPreprocessing)) {
                assignFromKeyfile(keyFile, "RAW Bayer", "PreprocessingEnabled", raw.bayersensor.enable_preproc);
                assignFromKeyfile(keyFile, "RAW Bayer", "LineDenoise", raw.bayersensor.linenoise);

                if (keyFile.has_key("RAW Bayer", "LineDenoiseDirection")) {
                    raw.bayersensor.linenoiseDirection = RAWParams::BayerSensor::LineNoiseDirection(keyFile.get_integer("RAW Bayer", "LineDenoiseDirection"));
                }

                assignFromKeyfile(keyFile, "RAW Bayer", "GreenEqThreshold", raw.bayersensor.greenthresh);
            }

            if (RELEVANT_(demosaic)) {
                assignFromKeyfile(keyFile, "RAW Bayer", "DCBIterations", raw.bayersensor.dcb_iterations);
                assignFromKeyfile(keyFile, "RAW Bayer", "DCBEnhance", raw.bayersensor.dcb_enhance);
                assignFromKeyfile(keyFile, "RAW Bayer", "LMMSEIterations", raw.bayersensor.lmmse_iterations);
                assignFromKeyfile(keyFile, "RAW Bayer", "DualDemosaicAutoContrast", raw.bayersensor.dualDemosaicAutoContrast);
                if (ppVersion < 345) {
                    raw.bayersensor.dualDemosaicAutoContrast = false;
                }
                assignFromKeyfile(keyFile, "RAW Bayer", "DualDemosaicContrast", raw.bayersensor.dualDemosaicContrast);

                if (keyFile.has_key("RAW Bayer", "PixelShiftMotionCorrectionMethod")) {
                    raw.bayersensor.pixelShiftMotionCorrectionMethod = (RAWParams::BayerSensor::PSMotionCorrectionMethod)keyFile.get_integer("RAW Bayer", "PixelShiftMotionCorrectionMethod");
                }

                assignFromKeyfile(keyFile, "RAW Bayer", "PixelShiftEperIso", raw.bayersensor.pixelShiftEperIso);
                if (ppVersion < 332) {
                    raw.bayersensor.pixelShiftEperIso += 1.0;
                }
                assignFromKeyfile(keyFile, "RAW Bayer", "PixelShiftSigma", raw.bayersensor.pixelShiftSigma);
                assignFromKeyfile(keyFile, "RAW Bayer", "PixelShiftShowMotion", raw.bayersensor.pixelShiftShowMotion);
                assignFromKeyfile(keyFile, "RAW Bayer", "PixelShiftShowMotionMaskOnly", raw.bayersensor.pixelShiftShowMotionMaskOnly);
                assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftHoleFill", raw.bayersensor.pixelShiftHoleFill);
                assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftMedian", raw.bayersensor.pixelShiftMedian);
                assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftGreen", raw.bayersensor.pixelShiftGreen);
                assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftBlur", raw.bayersensor.pixelShiftBlur);
                assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftSmoothFactor", raw.bayersensor.pixelShiftSmoothFactor);
                assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftEqualBright", raw.bayersensor.pixelShiftEqualBright);
                assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftEqualBrightChannel", raw.bayersensor.pixelShiftEqualBrightChannel);
                assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftNonGreenCross", raw.bayersensor.pixelShiftNonGreenCross);

                if (ppVersion < 336) {
                    if (keyFile.has_key("RAW Bayer", "pixelShiftLmmse")) {
                        bool useLmmse = keyFile.get_boolean ("RAW Bayer", "pixelShiftLmmse");
                        if (useLmmse) {
                            raw.bayersensor.pixelShiftDemosaicMethod = raw.bayersensor.getPSDemosaicMethodString(RAWParams::BayerSensor::PSDemosaicMethod::LMMSE);
                        } else {
                            raw.bayersensor.pixelShiftDemosaicMethod = raw.bayersensor.getPSDemosaicMethodString(RAWParams::BayerSensor::PSDemosaicMethod::AMAZE);
                        }
                    }
                } else {
                    assignFromKeyfile(keyFile, "RAW Bayer", "pixelShiftDemosaicMethod", raw.bayersensor.pixelShiftDemosaicMethod);
                }
            }

            if (RELEVANT_(rawPreprocessing)) {
                assignFromKeyfile(keyFile, "RAW Bayer", "PDAFLinesFilter", raw.bayersensor.pdafLinesFilter);
            }
        }

        if (keyFile.has_group("RAW X-Trans")) {
            if (RELEVANT_(demosaic)) {
                assignFromKeyfile(keyFile, "RAW X-Trans", "Method", raw.xtranssensor.method);
                assignFromKeyfile(keyFile, "RAW X-Trans", "DualDemosaicAutoContrast", raw.xtranssensor.dualDemosaicAutoContrast);
                if (ppVersion < 345) {
                    raw.xtranssensor.dualDemosaicAutoContrast = false;
                }
                assignFromKeyfile(keyFile, "RAW X-Trans", "DualDemosaicContrast", raw.xtranssensor.dualDemosaicContrast);
                assignFromKeyfile(keyFile, "RAW X-Trans", "Border", raw.xtranssensor.border);
                assignFromKeyfile(keyFile, "RAW X-Trans", "CcSteps", raw.xtranssensor.ccSteps);
            }
            if (RELEVANT_(rawBlack)) {
                assignFromKeyfile(keyFile, "RAW X-Trans", "PreBlackEnabled", raw.xtranssensor.enable_black);
                assignFromKeyfile(keyFile, "RAW X-Trans", "PreBlackRed", raw.xtranssensor.blackred);
                assignFromKeyfile(keyFile, "RAW X-Trans", "PreBlackGreen", raw.xtranssensor.blackgreen);
                assignFromKeyfile(keyFile, "RAW X-Trans", "PreBlackBlue", raw.xtranssensor.blackblue);
            }
        }

        if (keyFile.has_group("MetaData") && RELEVANT_(metadata)) {
            int mode = int(MetaDataParams::TUNNEL);
            assignFromKeyfile(keyFile, "MetaData", "Mode", mode);

            if (mode >= int(MetaDataParams::TUNNEL) && mode <= int(MetaDataParams::STRIP)) {
                metadata.mode = static_cast<MetaDataParams::Mode>(mode);
            }
        }

        if (keyFile.has_group("Exif") && RELEVANT_(exif)) {
            for (const auto& key : keyFile.get_keys("Exif")) {
                auto it = exif_keys.find(key);
                if (it != exif_keys.end()) {
                    exif[it->second] = keyFile.get_string("Exif", key);
                }
            }
        }

        /*
         * Load iptc change settings
         *
         * Existing values are preserved, and the stored values
         * are added to the list. To reset a field, the user has to
         * save the profile with the field leaved empty, but still
         * terminated by a semi-column ";"
         *
         * Please note that the old Keywords and SupplementalCategories
         * tag content is fully replaced by the new one,
         * i.e. they don't merge
         */
        if (keyFile.has_group("IPTC") && RELEVANT_(iptc)) {
            for (const auto& key : keyFile.get_keys("IPTC")) {
                // does this key already exist?
                auto it = iptc_keys.find(key);
                if (it == iptc_keys.end()) {
                    continue;
                }
                auto kk = it->second;
                const IPTCPairs::iterator element = iptc.find(kk);

                if (element != iptc.end()) {
                    // it already exist so we cleanup the values
                    element->second.clear();
                }

                // TODO: look out if merging Keywords and SupplementalCategories from the procparams chain would be interesting
                for (const auto& currLoadedTagValue : keyFile.get_string_list("IPTC", key)) {
                    iptc[kk].push_back(currLoadedTagValue);
                }
            }
        }

        return 0;
    } catch (const Glib::Error& e) {
        printf("-->%s\n", e.what().c_str());
        if (resetOnError) {
            setDefaults();
        }
        return 1;
    } catch (...) {
        printf("-->unknown exception!\n");
        if (resetOnError) {
            setDefaults();
        }
        return 1;
    }

    return 0;

#undef RELEVANT_
}


ProcParams* ProcParams::create()
{
    return new ProcParams();
}

void ProcParams::destroy(ProcParams* pp)
{
    delete pp;
}

bool ProcParams::operator ==(const ProcParams& other) const
{
    return
        exposure == other.exposure
        && brightContrSat == other.brightContrSat
        && toneCurve == other.toneCurve
        && localContrast == other.localContrast
        && labCurve == other.labCurve
        && sharpenMicro == other.sharpenMicro
        && sharpening == other.sharpening
        && prsharpening == other.prsharpening
        && wb == other.wb
        && impulseDenoise == other.impulseDenoise
        && denoise == other.denoise
        && epd == other.epd
        && fattal == other.fattal
        && logenc == other.logenc
        && defringe == other.defringe
        && sh == other.sh
        && toneEqualizer == other.toneEqualizer
        && crop == other.crop
        && coarse == other.coarse
        && rotate == other.rotate
        && commonTrans == other.commonTrans
        && distortion == other.distortion
        && lensProf == other.lensProf
        && perspective == other.perspective
        && gradient == other.gradient
        && pcvignette == other.pcvignette
        && cacorrection == other.cacorrection
        && vignetting == other.vignetting
        && chmixer == other.chmixer
        && blackwhite == other.blackwhite
        && resize == other.resize
        && raw == other.raw
        && icm == other.icm
        && dirpyrequalizer == other.dirpyrequalizer
        && filmSimulation == other.filmSimulation
        && softlight == other.softlight
        && rgbCurves == other.rgbCurves
        && metadata == other.metadata
        && exif == other.exif
        && iptc == other.iptc
        && dehaze == other.dehaze
        && grain == other.grain
        && smoothing == other.smoothing
        && colorcorrection == other.colorcorrection;
}

bool ProcParams::operator !=(const ProcParams& other) const
{
    return !(*this == other);
}

void ProcParams::init()
{
}

void ProcParams::cleanup()
{
}

int ProcParams::write(const Glib::ustring& fname, const Glib::ustring& content) const
{
    int error = 0;

    if (fname.length()) {
        FILE *f;
        f = g_fopen(fname.c_str(), "wt");

        if (f == nullptr) {
            error = 1;
        } else {
            fprintf(f, "%s", content.c_str());
            fclose(f);
        }
    }

    return error;
}


FullPartialProfile::FullPartialProfile(const ProcParams &pp):
    pp_(pp)
{
}


bool FullPartialProfile::applyTo(ProcParams &pp) const
{
    pp = pp_;
    return true;
}


FilePartialProfile::FilePartialProfile(const Glib::ustring &fname, bool full):
    fname_(fname),
    full_(full)
{
}


bool FilePartialProfile::applyTo(ProcParams &pp) const
{
    if (full_) {
        pp.setDefaults();
    }
    return fname_.empty() || (pp.load(fname_) == 0);
}


PEditedPartialProfile::PEditedPartialProfile(const Glib::ustring &fname, const ParamsEdited &pe):
    fname_(fname),
    pp_(),
    pe_(pe)
{
}


PEditedPartialProfile::PEditedPartialProfile(const ProcParams &pp, const ParamsEdited &pe):
    fname_(""),
    pp_(pp),
    pe_(pe)
{
}


bool PEditedPartialProfile::applyTo(ProcParams &pp) const
{
    if (!fname_.empty()) {
        Glib::KeyFile keyfile;
        try {
            if (!Glib::file_test(fname_, Glib::FILE_TEST_EXISTS) ||
                !keyfile.load_from_file(fname_)) {
                return false;
            }
        } catch (const Glib::Error& e) {
            printf("-->%s\n", e.what().c_str());
            return false;
        }
        return pp.load(keyfile, &pe_, false) == 0;
    } else {
        Glib::KeyFile keyfile;
        if (pp_.save(keyfile, &pe_) == 0) {
            return pp.load(keyfile, &pe_, false) == 0;
        }
    }
    return false;
}

}} // namespace rtengine::procparams

