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
#include <functional>
#include <strings.h>
#include <glib/gstdio.h>
#include <tiff.h>
#include <exiv2/exiv2.hpp>

#include "imagedata.h"
//#include "iptcpairs.h"
#include "imagesource.h"
#include "rt_math.h"
#pragma GCC diagnostic warning "-Wextra"
#define PRINT_HDR_PS_DETECTION 0

using namespace rtengine;

// extern "C" IptcData *iptc_data_new_from_jpeg_file (FILE* infile);

namespace {

Glib::ustring to_utf8 (const std::string& str)
{
    try {
        return Glib::locale_to_utf8 (str);
    } catch (Glib::Error&) {
        return Glib::convert_with_fallback (str, "UTF-8", "ISO-8859-1", "?");
    }
}

template<typename T>
T getFromFrame(
    const std::vector<std::unique_ptr<FrameData>>& frames,
    std::size_t frame,
    const std::function<T (const FrameData&)>& function
)
{
    if (frame < frames.size()) {
        return function(*frames[frame]);
    }
    if (!frames.empty()) {
        return function(*frames[0]);
    }
    return {};
}

} // namespace


namespace rtengine {

extern const Settings *settings;

Exiv2::Image::AutoPtr open_exiv2(const Glib::ustring &fname)
{
#ifdef EXV_UNICODE_PATH
    auto *ws = g_utf8_to_utf16(fname.c_str(), -1, NULL, NULL, NULL);
    std::wstring wfname(ws);
    g_free(ws);
    auto image = Exiv2::ImageFactory::open(wfname);
#else
    auto image = Exiv2::ImageFactory::open(fname);
#endif
    return image;
}

} // namespace rtengine



FramesMetaData* FramesMetaData::fromFile (const Glib::ustring& fname, std::unique_ptr<RawMetaDataLocation> rml, bool firstFrameOnly)
{
    return new FramesData (fname, std::move(rml), firstFrameOnly);
}

// FrameData::FrameData(rtexif::TagDirectory* frameRootDir_, rtexif::TagDirectory* rootDir, rtexif::TagDirectory* firstRootDir)
FrameData::FrameData(const Glib::ustring &fname):
    ok_(false),
    time(),
    timeStamp(),
    iso_speed(0),
    aperture(0.),
    focal_len(0.),
    focal_len35mm(0.),
    focus_dist(0.f),
    shutter(0.),
    expcomp(0.),
    make("Unknown"),
    model("Unknown"),
    orientation("Unknown"),
    lens("Unknown"),
    sampleFormat(IIOSF_UNKNOWN),
    isPixelShift(false),
    isHDR(false)
{
    memset(&time, 0, sizeof(time));

    // if (!frameRootDir) {
    //     return;
    // }

    // rtexif::Tag* tag;
    // rtexif::TagDirectory* newFrameRootDir = frameRootDir;

    memset(&time, 0, sizeof(time));
    timeStamp = 0;
    iso_speed = 0;
    aperture = 0.0;
    focal_len = 0.0;
    focal_len35mm = 0.0;
    focus_dist = 0.0f;
    shutter = 0.0;
    expcomp = 0.0;
    make.clear();
    model.clear();
    serial.clear();
    orientation.clear();
    lens.clear();

    try {
        auto image = open_exiv2(fname);
        image->readMetadata();
        auto &exif = image->exifData();
        ok_ = true;

        // taken and adapted from darktable (src/common/exif.cc)
/*
   This file is part of darktable,
   copyright (c) 2009--2013 johannes hanika.
   copyright (c) 2011 henrik andersson.
   copyright (c) 2012-2017 tobias ellinghaus.

   darktable is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   darktable is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with darktable.  If not, see <http://www.gnu.org/licenses/>.
 */
        
        Exiv2::ExifData::const_iterator pos;
    
        const auto find_exif_tag =
            [&](const std::string &name) -> bool
            {
                pos = exif.findKey(Exiv2::ExifKey(name));
                return (pos != exif.end() && pos->size());
            };

        const auto find_tag =
            [&](decltype(Exiv2::make) func) -> bool
            {
                pos = func(exif);
                return pos != exif.end() && pos->size();
            };

        /* List of tag names taken from exiv2's printSummary() in actions.cpp */

        if (find_tag(Exiv2::make)) {
            make = pos->print(&exif);
        }
        
        if (find_tag(Exiv2::model)) {
            model = pos->print(&exif);
        }

        if (make.size() > 0) {
            for (const auto& corp : {
                    "Canon",
                    "NIKON",
                    "EPSON",
                    "KODAK",
                    "Kodak",
                    "OLYMPUS",
                    "PENTAX",
                    "RICOH",
                    "MINOLTA",
                    "Minolta",
                    "Konica",
                    "CASIO",
                    "Sinar",
                    "Phase One",
                    "SAMSUNG",
                    "Mamiya",
                    "MOTOROLA",
                    "Leaf",
                    "Panasonic"
                  }) {
                if (make.find(corp) != std::string::npos) { // Simplify company names
                    make = corp;
                    break;
                }
            }
        }            
        make.erase(make.find_last_not_of(' ') + 1);
        model.erase(model.find_last_not_of(' ') + 1);

        if (make.length() > 0 && model.find(make + " ") == 0) {
            model = model.substr(make.length() + 1);
        }

        if (find_tag(Exiv2::exposureTime)) {
            shutter = pos->toFloat();
        }

        if (find_tag(Exiv2::fNumber)) {
            aperture = pos->toFloat();
        }

        /* Read ISO speed - Nikon happens to return a pair for Lo and Hi modes */
        if (find_tag(Exiv2::isoSpeed)) {
            // if standard exif iso tag, use the old way of interpreting the return value to be more regression-save
            if (strcmp(pos->key().c_str(), "Exif.Photo.ISOSpeedRatings") == 0) {
                int isofield = pos->count() > 1 ? 1 : 0;
                iso_speed = pos->toFloat(isofield);
            } else {
                std::string str = pos->print();
                iso_speed = std::atof(str.c_str());
            }
        }
        // some newer cameras support iso settings that exceed the 16 bit of exif's ISOSpeedRatings
        if (iso_speed == 65535 || iso_speed == 0) {
            if (find_exif_tag("Exif.PentaxDng.ISO") || find_exif_tag("Exif.Pentax.ISO")) {
                std::string str = pos->print();
                iso_speed = std::atof(str.c_str());
            } else if((!g_strcmp0(make.c_str(), "SONY") || !g_strcmp0(make.c_str(), "Canon"))
                    && find_exif_tag("Exif.Photo.RecommendedExposureIndex")) {
                iso_speed = pos->toFloat();
            }
        }

        if (find_tag(Exiv2::focalLength)) {
            // This works around a bug in exiv2 the developers refuse to fix
            // For details see http://dev.exiv2.org/issues/1083
            if (pos->key() == "Exif.Canon.FocalLength" && pos->count() == 4) {
                focal_len = pos->toFloat(1);
            } else {
                focal_len = pos->toFloat();
            }
        }

        if (find_exif_tag("Exif.Photo.FocalLengthIn35mmFilm")) {
            focal_len35mm = pos->toFloat();
        }

        if (find_tag(Exiv2::subjectDistance)) {
            focus_dist = (0.01 * pow(10, pos->toFloat() / 40));
        }
        
        if (find_tag(Exiv2::orientation)) {
            orientation = pos->print(&exif);
        }

        if (find_tag(Exiv2::lensName)) {
            lens = pos->print(&exif);
        }

        // /* Read lens name */
        // if ((find_exif_tag("Exif.CanonCs.LensType") && pos->print(&exif) != "(0)"
        //     && pos->print(&exif) != "(65535)")
        //     || find_exif_tag("Exif.Canon.0x0095")) {
        //     lens = pos->print(&exif);
        // } else if (EXIV2_MAKE_VERSION(0,25,0) <= Exiv2::versionNumber() && find_exif_tag("Exif.PentaxDng.LensType")) {
        //     lens = pos->print(&exif);
        // } else if (find_exif_tag("Exif.Panasonic.LensType")) {
        //     lens = pos->print(&exif);
        // } else if(find_exif_tag("Exif.OlympusEq.LensType")) {
        //     /* For every Olympus camera Exif.OlympusEq.LensType is present. */
        //     lens = pos->print(&exif);

        //     /* We have to check if Exif.OlympusEq.LensType has been translated by
        //      * exiv2. If it hasn't, fall back to Exif.OlympusEq.LensModel. */
        //     if(std::string::npos == lens.find_first_not_of(" 1234567890")) {
        //         /* Exif.OlympusEq.LensType contains only digits and spaces.
        //          * This means that exiv2 couldn't convert it to human readable
        //          * form. */
        //         if (find_exif_tag("Exif.OlympusEq.LensModel")) {
        //             lens = pos->print(&exif);
        //         } else if(find_exif_tag("Exif.Photo.LensModel")) {
        //         /* Just in case Exif.OlympusEq.LensModel hasn't been found */
        //             lens = pos->print(&exif);
        //         }
        //         //fprintf(stderr, "[exif] Warning: lens \"%s\" unknown as \"%s\"\n", img->exif_lens, lens.c_str());
        //     }
        // } else if ((pos = Exiv2::lensName(exif)) != exif.end() && pos->size()) {
        //     lens = pos->print(&exif);
        // } else if (find_exif_tag("Exif.Photo.LensModel")) {
        //     lens = pos->print(&exif);
        // }

        std::string datetime_taken;
        if (find_exif_tag("Exif.Image.DateTimeOriginal")) {
            datetime_taken = pos->print(&exif);
        } else if(find_exif_tag("Exif.Photo.DateTimeOriginal")) {
            datetime_taken = pos->print(&exif);
        }
        if (sscanf(datetime_taken.c_str(), "%d:%d:%d %d:%d:%d", &time.tm_year, &time.tm_mon, &time.tm_mday, &time.tm_hour, &time.tm_min, &time.tm_sec) == 6) {
            time.tm_year -= 1900;
            time.tm_mon -= 1;
            time.tm_isdst = -1;
            timeStamp = mktime(&time);
        }

        // Improve lens detection for Sony lenses.
        if (find_exif_tag("Exif.Sony2.LensID") && pos->toLong() != 65535 && pos->print().find('|') == std::string::npos) {
            lens = pos->print(&exif);
        } else if((!strncmp(model.c_str(), "NEX", 3)) || (!strncmp(model.c_str(), "ILCE", 4))) {
            // Workaround for an issue on newer Sony NEX cams.
            // The default EXIF field is not used by Sony to store lens data
            // http://dev.exiv2.org/issues/833
            // http://darktable.org/redmine/issues/8813
            // FIXME: This is still a workaround
            lens = "Unknown";
            if (find_exif_tag("Exif.Photo.LensModel")) {
                lens = pos->print(&exif);
            }
        }

        if (find_exif_tag("Exif.Image.ExposureBiasValue")) {
            expcomp = pos->toFloat();
        }

        // -----------------------
        // Special file type detection (HDR, PixelShift)
        // ------------------------
        uint16 bitspersample = 0, samplesperpixel = 0, sampleformat = 0, photometric = 0, compression = 0;
        auto bps = exif.findKey(Exiv2::ExifKey("Exif.Image.BitsPerSample"));
        auto spp = exif.findKey(Exiv2::ExifKey("Exif.Image.SamplesPerPixel"));
        auto sf = exif.findKey(Exiv2::ExifKey("Exif.Image.SampleFormat"));
        auto pi = exif.findKey(Exiv2::ExifKey("Exif.Image.PhotometricInterpretation"));
        auto c = exif.findKey(Exiv2::ExifKey("Exif.Image.Compression"));

        if ((!make.compare (0, 6, "PENTAX") || (!make.compare (0, 5, "RICOH") && !model.compare (0, 6, "PENTAX")))) {
//             if (find_exif_tag("Exif.Pentax.HDR") && pos->toLong() > 0) {
//                 isHDR = true;
// #if PRINT_HDR_PS_DETECTION
//                 printf("HDR detected ! -> \"HDR\" tag found\n");
// #endif
//             } else
            if (find_exif_tag("Exif.Pentax.DriveMode")) {
                std::string buf = pos->toString(3);
                buf[3] = 0;
                if (!strcmp(buf.c_str(), "HDR")) {
                    isHDR = true;
#if PRINT_HDR_PS_DETECTION
                    printf("HDR detected ! -> DriveMode = \"HDR\"\n");
#endif
                }
            }

            if (!isHDR && find_exif_tag("Exif.Pentax.Quality") &&
                (pos->toLong() == 7 || pos->toLong() == 8)) {
                isPixelShift = true;
#if PRINT_HDR_PS_DETECTION
                printf("PixelShift detected ! -> \"Quality\" = 7\n");
#endif
            }
        }

        sampleFormat = IIOSF_UNKNOWN;

        if (sf == exif.end())
            /*
             * WARNING: This is a dirty hack!
             * We assume that files which doesn't contain the TIFFTAG_SAMPLEFORMAT tag
             * (which is the case with uncompressed TIFFs produced by RT!) are RGB files,
             * but that may be not true.   --- Hombre
             */
        {
            sampleformat = SAMPLEFORMAT_UINT;
        } else {
            sampleformat = sf->toLong();
        }

        if (bps == exif.end() || spp == exif.end() || pi == exif.end()) {
            return;
        }

        bitspersample = bps->toLong();
        samplesperpixel = spp->toLong();

        photometric = pi->toLong();
        if (photometric == PHOTOMETRIC_LOGLUV) {
            if (c == exif.end()) {
                compression = COMPRESSION_NONE;
            } else {
                compression = c->toLong();
            }
        }

        if (photometric == PHOTOMETRIC_RGB || photometric == PHOTOMETRIC_MINISBLACK) {
            if (sampleformat == SAMPLEFORMAT_INT || sampleformat == SAMPLEFORMAT_UINT) {
                if (bitspersample == 8) {
                    sampleFormat = IIOSF_UNSIGNED_CHAR;
                } else if (bitspersample <= 16) {
                    sampleFormat = IIOSF_UNSIGNED_SHORT;
                }
            } else if (sampleformat == SAMPLEFORMAT_IEEEFP) {
                if (bitspersample==16) {
                    sampleFormat = IIOSF_FLOAT16;
                    isHDR = true;
#if PRINT_HDR_PS_DETECTION
                    printf("HDR detected ! -> sampleFormat = %d   (16-bit)\n", sampleFormat);
#endif
                }
                else if (bitspersample == 24) {
                    sampleFormat = IIOSF_FLOAT24;
                    isHDR = true;
#if PRINT_HDR_PS_DETECTION
                    printf("HDR detected ! -> sampleFormat = %d   (24-bit)\n", sampleFormat);
#endif
                }
                else if (bitspersample == 32) {
                    sampleFormat = IIOSF_FLOAT32;
                    isHDR = true;
#if PRINT_HDR_PS_DETECTION
                    printf("HDR detected ! -> sampleFormat = %d   (32-bit)\n", sampleFormat);
#endif
                }
            }
        } else if (photometric == PHOTOMETRIC_CFA) {
            if (sampleformat == SAMPLEFORMAT_IEEEFP) {
                if (bitspersample == 16) {
                    sampleFormat = IIOSF_FLOAT16;
                    isHDR = true;
#if PRINT_HDR_PS_DETECTION
                    printf("HDR detected ! -> sampleFormat = %d   (16-bit)\n", sampleFormat);
#endif
                }
                else if (bitspersample == 24) {
                    sampleFormat = IIOSF_FLOAT24;
                    isHDR = true;
#if PRINT_HDR_PS_DETECTION
                    printf("HDR detected ! -> sampleFormat = %d   (24-bit)\n", sampleFormat);
#endif
                }
                else if (bitspersample == 32) {
                    sampleFormat = IIOSF_FLOAT32;
                    isHDR = true;
#if PRINT_HDR_PS_DETECTION
                    printf("HDR detected ! -> sampleFormat = %d   (32-bit)\n", sampleFormat);
#endif
                }
            } else if (sampleformat == SAMPLEFORMAT_INT || sampleformat == SAMPLEFORMAT_UINT) {
                if (bitspersample == 8) {   // shouldn't occur...
                    sampleFormat = IIOSF_UNSIGNED_CHAR;
                } else if (bitspersample <= 16) {
                    sampleFormat = IIOSF_UNSIGNED_SHORT;
                }
            }
        } else if (photometric == 34892 || photometric == 32892  /* Linear RAW (see DNG spec ; 32892 seem to be a flaw from Sony's ARQ files) */) {
            if (sampleformat == SAMPLEFORMAT_IEEEFP) {
                sampleFormat = IIOSF_FLOAT32;
                isHDR = true;
#if PRINT_HDR_PS_DETECTION
                printf("HDR detected ! -> sampleFormat = %d\n", sampleFormat);
#endif
            } else if (sampleformat == SAMPLEFORMAT_INT || sampleformat == SAMPLEFORMAT_UINT) {
                if (bitspersample == 8) {   // shouldn't occur...
                    sampleFormat = IIOSF_UNSIGNED_CHAR;
                } else if (bitspersample <= 16) {
                    sampleFormat = IIOSF_UNSIGNED_SHORT;
                    if (find_exif_tag("Exif.Photo.MakerNote") && (!make.compare (0, 4, "SONY")) && bitspersample >= 12 && samplesperpixel == 4) {
                        isPixelShift = true;
#if PRINT_HDR_PS_DETECTION
                        printf("PixelShift detected ! -> \"Make\" = SONY, bitsPerPixel > 8, samplesPerPixel == 4\n");
#endif
                    }
                }
            }
        } else if (photometric == PHOTOMETRIC_LOGLUV) {
            if (compression == COMPRESSION_SGILOG24) {
                sampleFormat = IIOSF_LOGLUV24;
                isHDR = true;
#if PRINT_HDR_PS_DETECTION
                printf("HDR detected ! -> sampleFormat = %d\n", sampleFormat);
#endif
            } else if (compression == COMPRESSION_SGILOG) {
                sampleFormat = IIOSF_LOGLUV32;
                isHDR = true;
#if PRINT_HDR_PS_DETECTION
                printf("HDR detected ! -> sampleFormat = %d\n", sampleFormat);
#endif
            }
        }
    } catch(Exiv2::AnyError &e) {
        if (settings->verbose) {
            std::cerr << "EXIV2 ERROR: " << e.what() << std::endl;
        }
        ok_ = false;
    }
}

FrameData::~FrameData ()
{

    // if (iptc) {
    //     iptc_data_free (iptc);
    // }
}

// procparams::IPTCPairs FrameData::getIPTCData () const
// {
//     return getIPTCData(iptc);
// }

// procparams::IPTCPairs FrameData::getIPTCData (IptcData* iptc_)
// {

//     procparams::IPTCPairs iptcc;

//     if (!iptc_) {
//         return iptcc;
//     }

//     unsigned char buffer[2100];

//     for (int i = 0; i < 16; i++) {
//         IptcDataSet* ds = iptc_data_get_next_dataset (iptc_, nullptr, IPTC_RECORD_APP_2, strTags[i].tag);

//         if (ds) {
//             iptc_dataset_get_data (ds, buffer, 2100);
//             std::vector<Glib::ustring> icValues;
//             icValues.push_back (to_utf8((char*)buffer));

//             iptcc[strTags[i].field] = icValues;
//             iptc_dataset_unref (ds);
//         }
//     }

//     IptcDataSet* ds = nullptr;
//     std::vector<Glib::ustring> keywords;

//     while ((ds = iptc_data_get_next_dataset (iptc_, ds, IPTC_RECORD_APP_2, IPTC_TAG_KEYWORDS))) {
//         iptc_dataset_get_data (ds, buffer, 2100);
//         keywords.push_back (to_utf8((char*)buffer));
//     }

//     iptcc["Keywords"] = keywords;
//     ds = nullptr;
//     std::vector<Glib::ustring> suppCategories;

//     while ((ds = iptc_data_get_next_dataset (iptc_, ds, IPTC_RECORD_APP_2, IPTC_TAG_SUPPL_CATEGORY))) {
//         iptc_dataset_get_data (ds, buffer, 2100);
//         suppCategories.push_back (to_utf8((char*)buffer));
//         iptc_dataset_unref (ds);
//     }

//     iptcc["SupplementalCategories"] = suppCategories;
//     return iptcc;
// }


bool FrameData::getPixelShift () const
{
    return isPixelShift;
}
bool FrameData::getHDR () const
{
    return isHDR;
}
std::string FrameData::getImageType () const
{
    return isPixelShift ? "PS" : isHDR ? "HDR" : "STD";
}
IIOSampleFormat FrameData::getSampleFormat () const
{
    return sampleFormat;
}
// rtexif::TagDirectory* FrameData::getExifData () const
// {
//     return frameRootDir;
// }
bool FrameData::hasExif () const
{
    return ok_; //frameRootDir && frameRootDir->getCount();
}
// bool FrameData::hasIPTC () const
// {
//     return iptc;
// }
tm FrameData::getDateTime () const
{
    return time;
}
time_t FrameData::getDateTimeAsTS () const
{
    return timeStamp;
}
int FrameData::getISOSpeed () const
{
    return iso_speed;
}
double FrameData::getFNumber () const
{
    return aperture;
}
double FrameData::getFocalLen () const
{
    return focal_len;
}
double FrameData::getFocalLen35mm () const
{
    return focal_len35mm;
}
float FrameData::getFocusDist () const
{
    return focus_dist;
}
double FrameData::getShutterSpeed () const
{
    return shutter;
}
double FrameData::getExpComp () const
{
    return expcomp;
}
std::string FrameData::getMake () const
{
    return make;
}
std::string FrameData::getModel () const
{
    return model;
}
std::string FrameData::getLens () const
{
    return lens;
}
std::string FrameData::getSerialNumber () const
{
    return serial;
}
std::string FrameData::getOrientation () const
{
    return orientation;
}



void FramesData::setDCRawFrameCount (unsigned int frameCount)
{
    dcrawFrameCount = frameCount;
}

// unsigned int FramesData::getRootCount () const
// {
//     return roots.size();
// }

unsigned int FramesData::getFrameCount () const
{
    return dcrawFrameCount ? dcrawFrameCount : frames.size();
}

bool FramesData::getPixelShift () const
{
    // So far only Pentax and Sony provide multi-frame Pixel Shift files.
    // Only the first frame contains the Pixel Shift tag
    // If more brand have to be supported, this rule may need
    // to evolve

    return frames.empty() ? false : frames.at(0)->getPixelShift ();
}
bool FramesData::getHDR (unsigned int frame) const
{
    // So far only Pentax provides multi-frame HDR file.
    // Only the first frame contains the HDR tag
    // If more brand have to be supported, this rule may need
    // to evolve

    return frames.empty() || frame >= frames.size()  ? false : frames.at(0)->getHDR ();
}

std::string FramesData::getImageType (unsigned int frame) const
{
    return frames.empty() || frame >= frames.size() ? "STD" : frames.at(0)->getImageType();
}

IIOSampleFormat FramesData::getSampleFormat (unsigned int frame) const
{
    return frames.empty() || frame >= frames.size()  ? IIOSF_UNKNOWN : frames.at(frame)->getSampleFormat ();
}

// rtexif::TagDirectory* FramesData::getFrameExifData (unsigned int frame) const
// {
//     return frames.empty() || frame >= frames.size()  ? nullptr : frames.at(frame)->getExifData ();
// }

// rtexif::TagDirectory* FramesData::getBestExifData (ImageSource *imgSource, procparams::RAWParams *rawParams) const
// {
//     rtexif::TagDirectory *td = nullptr;
//     if (frames.empty()) {
//         return nullptr;
//     }
//     if (imgSource && rawParams) {
//         eSensorType sensorType = imgSource->getSensorType();
//         unsigned int imgNum = 0;
//         if (sensorType == ST_BAYER) {
//             imgNum = rtengine::LIM<unsigned int>(rawParams->bayersensor.imageNum, 0, frames.size() - 1);
//         /*
//         // might exist someday ?
//         } else if (sensorType == ST_FUJI_XTRANS) {
//             imgNum = rtengine::LIM<unsigned int>(rawParams->xtranssensor.imageNum, 0, frames.size() - 1);
//         } else if (sensorType == ST_NONE && !imgSource->isRAW()) {
//             // standard image multiframe support should come here (when implemented in GUI)
//         */
//         }

//         td = getFrameExifData (imgNum);
//         rtexif::Tag* makeTag;
//         if (td && (makeTag = td->findTag("Make", true))) {
//             td = makeTag->getParent();
//         } else {
//             td = getRootExifData(0);
//         }
//     }
//     return td;
// }

// rtexif::TagDirectory* FramesData::getRootExifData (unsigned int root) const
// {
//     return roots.empty() || root >= roots.size()  ? nullptr : roots.at(root);
// }

// procparams::IPTCPairs FramesData::getIPTCData (unsigned int frame) const
// {
//     if (frame < frames.size() && frames.at(frame)->hasIPTC()) {
//         return frames.at(frame)->getIPTCData();
//     } else {
//         if (iptc) {
//             return FrameData::getIPTCData(iptc);
//         } else {
//             procparams::IPTCPairs emptyPairs;
//             return emptyPairs;
//         }
//     }
// }

bool FramesData::hasExif(unsigned int frame) const
{
    return getFromFrame<bool>(
        frames,
        frame,
        [](const FrameData& frame_data)
        {
            return frame_data.hasExif();
        }
    );
}

// bool FramesData::hasIPTC(unsigned int frame) const
// {
//     return getFromFrame<bool>(
//         frames,
//         frame,
//         [](const FrameData& frame_data)
//         {
//             return frame_data.hasIPTC();
//         }
//     );
// }

tm FramesData::getDateTime(unsigned int frame) const
{
    return getFromFrame<tm>(
        frames,
        frame,
        [](const FrameData& frame_data)
        {
            return frame_data.getDateTime();
        }
    );
}

time_t FramesData::getDateTimeAsTS(unsigned int frame) const
{
    return getFromFrame<time_t>(
        frames,
        frame,
        [](const FrameData& frame_data)
        {
            return frame_data.getDateTimeAsTS();
        }
    );
}

int FramesData::getISOSpeed(unsigned int frame) const
{
    return getFromFrame<int>(
        frames,
        frame,
        [](const FrameData& frame_data)
        {
            return frame_data.getISOSpeed();
        }
    );
}

double FramesData::getFNumber(unsigned int frame) const
{
    return getFromFrame<double>(
        frames,
        frame,
        [](const FrameData& frame_data)
        {
            return frame_data.getFNumber();
        }
    );
}

double FramesData::getFocalLen(unsigned int frame) const
{
    return getFromFrame<double>(
        frames,
        frame,
        [](const FrameData& frame_data)
        {
            return frame_data.getFocalLen();
        }
    );
}

double FramesData::getFocalLen35mm(unsigned int frame) const
{
    return getFromFrame<double>(
        frames,
        frame,
        [](const FrameData& frame_data)
        {
            return frame_data.getFocalLen35mm();
        }
    );
}

float FramesData::getFocusDist(unsigned int frame) const
{
    return getFromFrame<float>(
        frames,
        frame,
        [](const FrameData& frame_data)
        {
            return frame_data.getFocusDist();
        }
    );
}

double FramesData::getShutterSpeed(unsigned int frame) const
{
    return getFromFrame<double>(
        frames,
        frame,
        [](const FrameData& frame_data)
        {
            return frame_data.getShutterSpeed();
        }
    );
}

double FramesData::getExpComp(unsigned int frame) const
{
    return getFromFrame<double>(
        frames,
        frame,
        [](const FrameData& frame_data)
        {
            return frame_data.getExpComp();
        }
    );
}

std::string FramesData::getMake(unsigned int frame) const
{
    return getFromFrame<std::string>(
        frames,
        frame,
        [](const FrameData& frame_data)
        {
            return frame_data.getMake();
        }
    );
}

std::string FramesData::getModel(unsigned int frame) const
{
    return getFromFrame<std::string>(
        frames,
        frame,
        [](const FrameData& frame_data)
        {
            return frame_data.getModel();
        }
    );
}

std::string FramesData::getLens(unsigned int frame) const
{
    return getFromFrame<std::string>(
        frames,
        frame,
        [](const FrameData& frame_data)
        {
            return frame_data.getLens();
        }
    );
}

std::string FramesData::getSerialNumber(unsigned int frame) const
{
    return getFromFrame<std::string>(
        frames,
        frame,
        [](const FrameData& frame_data)
        {
            return frame_data.getSerialNumber();
        }
    );
}

std::string FramesData::getOrientation(unsigned int frame) const
{
    return getFromFrame<std::string>(
        frames,
        frame,
        [](const FrameData& frame_data)
        {
            return frame_data.getOrientation();
        }
    );
}


//------inherited functions--------------//


std::string FramesMetaData::apertureToString (double aperture)
{

    char buffer[256];
    sprintf (buffer, "%0.1f", aperture);
    return buffer;
}

std::string FramesMetaData::shutterToString (double shutter)
{

    char buffer[256];

    if (shutter > 0.0 && shutter <= 0.5) {
        sprintf (buffer, "1/%0.0f", 1.0 / shutter);
    } else {
        sprintf (buffer, "%0.1f", shutter);
    }

    return buffer;
}

std::string FramesMetaData::expcompToString (double expcomp, bool maskZeroexpcomp)
{

    char buffer[256];

    if (maskZeroexpcomp) {
        if (expcomp != 0.0) {
            sprintf (buffer, "%0.2f", expcomp);
            return buffer;
        } else {
            return "";
        }
    } else {
        sprintf (buffer, "%0.2f", expcomp);
        return buffer;
    }
}

double FramesMetaData::shutterFromString (std::string s)
{

    size_t i = s.find_first_of ('/');

    if (i == std::string::npos) {
        return atof (s.c_str());
    } else {
        return atof (s.substr(0, i).c_str()) / atof (s.substr(i + 1).c_str());
    }
}

double FramesMetaData::apertureFromString (std::string s)
{

    return atof (s.c_str());
}

// extern "C" {

// #include <libiptcdata/iptc-data.h>
// #include <libiptcdata/iptc-jpeg.h>

//     struct _IptcDataPrivate {
//         unsigned int ref_count;

//         IptcLog *log;
//         IptcMem *mem;
//     };

//     IptcData *
//     iptc_data_new_from_jpeg_file (FILE *infile)
//     {
//         IptcData *d;
//         unsigned char * buf;
//         int buf_len = 256 * 256;
//         int len, offset;
//         unsigned int iptc_len;

//         if (!infile) {
//             return nullptr;
//         }

//         d = iptc_data_new ();

//         if (!d) {
//             return nullptr;
//         }

//         buf = (unsigned char*)iptc_mem_alloc (d->priv->mem, buf_len);

//         if (!buf) {
//             iptc_data_unref (d);
//             return nullptr;
//         }

//         len = iptc_jpeg_read_ps3 (infile, buf, buf_len);

//         if (len <= 0) {
//             goto failure;
//         }

//         offset = iptc_jpeg_ps3_find_iptc (buf, len, &iptc_len);

//         if (offset <= 0) {
//             goto failure;
//         }

//         iptc_data_load (d, buf + offset, iptc_len);

//         iptc_mem_free (d->priv->mem, buf);
//         return d;

// failure:
//         iptc_mem_free (d->priv->mem, buf);
//         iptc_data_unref (d);
//         return nullptr;
//     }

// }

FramesData::FramesData (const Glib::ustring& fname, std::unique_ptr<RawMetaDataLocation> rml, bool firstFrameOnly) :
    //iptc(nullptr),
    fname_(fname),
    dcrawFrameCount (0)
{
    frames.push_back(std::unique_ptr<FrameData>(new FrameData(fname)));
    // if (rml && (rml->exifBase >= 0 || rml->ciffBase >= 0)) {
    //     FILE* f = g_fopen (fname.c_str (), "rb");

    //     if (f) {
    //         const bool has_rml_exif_base = rml->exifBase >= 0;
    //         rtexif::ExifManager exifManager (f, std::move(rml), firstFrameOnly);

    //         if (has_rml_exif_base) {
    //             if (exifManager.f && exifManager.rml) {
    //                 if (exifManager.rml->exifBase >= 0) {
    //                     exifManager.parseRaw ();

    //                 } else if (exifManager.rml->ciffBase >= 0) {
    //                     exifManager.parseCIFF ();
    //                 }
    //             }

    //             // copying roots
    //             roots = exifManager.roots;

    //             // creating FrameData
    //             for (auto currFrame : exifManager.frames) {
    //                 frames.push_back(std::unique_ptr<FrameData>(new FrameData(currFrame, currFrame->getRoot(), roots.at(0))));
    //             }
    //             for (auto currRoot : roots) {
    //                 rtexif::Tag* t = currRoot->getTag(0x83BB);

    //                 if (t && !iptc) {
    //                     iptc = iptc_data_new_from_data ((unsigned char*)t->getValue (), (unsigned)t->getValueSize ());
    //                     break;
    //                 }
    //             }
    //         }
    //         fclose (f);
    //     }
    // } else if (hasJpegExtension(fname)) {
    //     FILE* f = g_fopen (fname.c_str (), "rb");

    //     if (f) {
    //         rtexif::ExifManager exifManager (f, std::move(rml), true);
    //         if (exifManager.f) {
    //             exifManager.parseJPEG ();
    //             roots = exifManager.roots;
    //             for (auto currFrame : exifManager.frames) {
    //                 frames.push_back(std::unique_ptr<FrameData>(new FrameData(currFrame, currFrame->getRoot(), roots.at(0))));
    //             }
    //             rewind (exifManager.f); // Not sure this is necessary
    //             iptc = iptc_data_new_from_jpeg_file (exifManager.f);
    //         }
    //         fclose (f);
    //     }
    // } else if (hasTiffExtension(fname)) {
    //     FILE* f = g_fopen (fname.c_str (), "rb");

    //     if (f) {
    //         rtexif::ExifManager exifManager (f, std::move(rml), firstFrameOnly);

    //         exifManager.parseTIFF();
    //         roots = exifManager.roots;

    //         // creating FrameData
    //         for (auto currFrame : exifManager.frames) {
    //             frames.push_back(std::unique_ptr<FrameData>(new FrameData(currFrame, currFrame->getRoot(), roots.at(0))));
    //         }
    //         for (auto currRoot : roots) {
    //             rtexif::Tag* t = currRoot->getTag(0x83BB);

    //             if (t && !iptc) {
    //                 iptc = iptc_data_new_from_data ((unsigned char*)t->getValue (), (unsigned)t->getValueSize ());
    //                 break;
    //             }
    //         }
    //         fclose (f);
    //     }
    // }
}

FramesData::~FramesData ()
{
    // for (auto currRoot : roots) {
    //     delete currRoot;
    // }

    // if (iptc) {
    //     iptc_data_free (iptc);
    // }
}


Glib::ustring FramesData::getFileName() const
{
    return fname_;
}
