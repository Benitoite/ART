/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *  Copyright (c) 2010 Oliver Duis <www.oliverduis.de>
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
#include <png.h>
#include <glib/gstdio.h>
#include <tiff.h>
#include <tiffio.h>
#include <cstdio>
#include <cstring>
#include <fcntl.h>
#include "rt_math.h"
#include "../rtgui/options.h"
#include "../rtgui/version.h"
#include "../rtgui/multilangmgr.h"
#include "myfile.h"

#ifdef WIN32
#include <winsock2.h>
#else
#include <netinet/in.h>
#endif

#include "imageio.h"
//#include "iptcpairs.h"
#include "iccjpeg.h"
#include "color.h"
#include "imagedata.h"
#include "settings.h"

#include "rtjpeg.h"

using namespace std;
using namespace rtengine;
using namespace rtengine::procparams;

namespace rtengine { extern const Settings *settings; }

namespace
{

// Opens a file for binary writing and request exclusive lock (cases were you need "wb" mode plus locking)
FILE* g_fopen_withBinaryAndLock(const Glib::ustring& fname)
{

#ifdef WIN32

    // Use native function to disallow sharing, i.e. lock the file for exclusive access.
    // This is important to e.g. prevent Windows Explorer from crashing RT due to concurrently scanning an image file.
    std::unique_ptr<wchar_t, GFreeFunc> wfname (reinterpret_cast<wchar_t*>(g_utf8_to_utf16 (fname.c_str (), -1, NULL, NULL, NULL)), g_free);

    HANDLE hFile = CreateFileW ( wfname.get (), GENERIC_READ | GENERIC_WRITE, 0 /* no sharing allowed */, NULL, CREATE_ALWAYS, FILE_ATTRIBUTE_NORMAL, NULL);
    FILE* f = nullptr;

    if (hFile != INVALID_HANDLE_VALUE) {
        f = _fdopen (_open_osfhandle ((intptr_t)hFile, 0), "wb");
    }

#else

    FILE* f = ::g_fopen (fname.c_str (), "wb");

#endif

    return f;
}

}

Glib::ustring ImageIO::errorMsg[6] = {"Success", "Cannot read file.", "Invalid header.", "Error while reading header.", "File reading error", "Image format not supported."};

void ImageIO::setOutputProfile(const char* pdata, int plen)
{

    delete [] profileData;

    if (pdata) {
        profileData = new char [plen];
        memcpy (profileData, pdata, plen);
    } else {
        profileData = nullptr;
    }

    profileLength = plen;
}

ImageIO::~ImageIO ()
{

    if (embProfile) {
        cmsCloseProfile(embProfile);
    }

    deleteLoadedProfileData();
    // delete exifRoot;
    delete [] profileData;
}

void png_read_data(png_struct_def  *png_ptr, unsigned char *data, size_t length);
void png_write_data(png_struct_def *png_ptr, unsigned char *data, size_t length);
void png_flush(png_struct_def *png_ptr);

int ImageIO::getPNGSampleFormat (const Glib::ustring &fname, IIOSampleFormat &sFormat, IIOSampleArrangement &sArrangement)
{
    FILE *file = g_fopen (fname.c_str (), "rb");

    if (!file) {
        return IMIO_CANNOTREADFILE;
    }

    //reading PNG header
    unsigned char header[8];

    if (fread (header, 1, 8, file) != 8 || png_sig_cmp (header, 0, 8)) {
        fclose(file);
        return IMIO_HEADERERROR;
    }

    //initializing main structures
    png_structp png = png_create_read_struct (PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);

    if (!png) {
        fclose (file);
        return IMIO_HEADERERROR;
    }

    png_infop info = png_create_info_struct (png);
    png_infop end_info = png_create_info_struct (png);

    if (!end_info || !info) {
        png_destroy_read_struct (&png, &info, &end_info);
        fclose (file);
        return IMIO_HEADERERROR;
    }

    if (setjmp (png_jmpbuf(png))) {
        png_destroy_read_struct (&png, &info, &end_info);
        fclose (file);
        return IMIO_READERROR;
    }

    //set up png read
    png_set_read_fn (png, file, png_read_data);
    png_set_sig_bytes (png, 8);

    png_read_info(png, info);

    //retrieving image information
    png_uint_32 width, height;
    int bit_depth, color_type, interlace_type, compression_type, filter_method;
    png_get_IHDR(png, info, &width, &height, &bit_depth, &color_type, &interlace_type, &compression_type, &filter_method);

    png_destroy_read_struct (&png, &info, &end_info);
    fclose (file);

    if (interlace_type != PNG_INTERLACE_NONE) {
        return IMIO_VARIANTNOTSUPPORTED;
    }

    if (bit_depth == 8) {
        sArrangement = IIOSA_CHUNKY;
        sFormat = IIOSF_UNSIGNED_CHAR;
        return IMIO_SUCCESS;
    } else if (bit_depth == 16) {
        sArrangement = IIOSA_CHUNKY;
        sFormat = IIOSF_UNSIGNED_SHORT;
        return IMIO_SUCCESS;
    } else {
        sArrangement = IIOSA_UNKNOWN;
        sFormat = IIOSF_UNKNOWN;
        return IMIO_VARIANTNOTSUPPORTED;
    }
}


int ImageIO::loadPNG(const Glib::ustring &fname)
{

    FILE *file = g_fopen (fname.c_str (), "rb");

    if (!file) {
        return IMIO_CANNOTREADFILE;
    }

    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_LOADPNG");
        pl->setProgress (0.0);
    }

    //reading PNG header
    unsigned char header[8];

    if (fread (header, 1, 8, file) != 8 || png_sig_cmp (header, 0, 8)) {
        fclose(file);
        return IMIO_HEADERERROR;
    }

    //initializing main structures
    png_structp png = png_create_read_struct (PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);

    if (!png) {
        fclose (file);
        return IMIO_HEADERERROR;
    }

    // silence the warning about "invalid" sRGB profiles -- see #4260
#if defined(PNG_SKIP_sRGB_CHECK_PROFILE) && defined(PNG_SET_OPTION_SUPPORTED)
    png_set_option(png, PNG_SKIP_sRGB_CHECK_PROFILE, PNG_OPTION_ON);
#endif

    png_infop info = png_create_info_struct (png);
    png_infop end_info = png_create_info_struct (png);

    if (!end_info || !info) {
        png_destroy_read_struct (&png, &info, &end_info);
        fclose (file);
        return IMIO_HEADERERROR;
    }

    if (setjmp (png_jmpbuf(png))) {
        png_destroy_read_struct (&png, &info, &end_info);
        fclose (file);
        return IMIO_READERROR;
    }

    //set up png read
    png_set_read_fn (png, file, png_read_data);
    png_set_sig_bytes (png, 8);

    png_read_info(png, info);

    embProfile = nullptr;

    //retrieving image information
    png_uint_32 width, height;
    int bit_depth, color_type, interlace_type, compression_type, filter_method;
    png_get_IHDR(png, info, &width, &height, &bit_depth, &color_type, &interlace_type, &compression_type, &filter_method);

    if (color_type == PNG_COLOR_TYPE_PALETTE || interlace_type != PNG_INTERLACE_NONE )  {
        // we don't support interlaced png or png with palette
        png_destroy_read_struct (&png, &info, &end_info);
        fclose (file);
        printf("%s uses an unsupported feature: <palette-indexed colors|interlacing>. Skipping.\n", fname.data());
        return IMIO_VARIANTNOTSUPPORTED;
    }

    if (color_type == PNG_COLOR_TYPE_GRAY || color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {
        png_set_gray_to_rgb(png);
    }

    if (png_get_valid(png, info, PNG_INFO_tRNS)) {
        png_set_tRNS_to_alpha(png);
    }

    if (color_type & PNG_COLOR_MASK_ALPHA) {
        png_set_strip_alpha(png);
    }

    // reading the embedded ICC profile if any
    if (png_get_valid(png, info, PNG_INFO_iCCP)) {
        png_charp name;
        int compression_type;
#if PNG_LIBPNG_VER < 10500
        png_charp profdata;
#else
        png_bytep profdata;
#endif
        png_uint_32 proflen;
        png_get_iCCP(png, info, &name, &compression_type, &profdata, &proflen);
        embProfile = cmsOpenProfileFromMem(profdata, proflen);
        loadedProfileData = new char[proflen];
        memcpy(loadedProfileData, profdata, proflen);
    }

    //setting gamma
    double gamma;

    if (png_get_gAMA(png, info, &gamma)) {
        png_set_gamma(png, 1.0 / gamma, gamma);    // use gamma from metadata
    } else {
        png_set_gamma(png, 2.2, 1.0 / 2.2);    // no gamma in metadata, suppose gamma 2.2
    }

//  if (bps==8 && bit_depth==16) png_set_strip_16(png);

    //updating png info struct
    png_read_update_info(png, info);
    png_get_IHDR(png, info, &width, &height, &bit_depth, &color_type, &interlace_type, &compression_type, &filter_method);

    allocate (width, height);

    int rowlen = width * 3 * bit_depth / 8;
    unsigned char *row = new unsigned char [rowlen];

    // set a new jump point to avoid memory leak
    if (setjmp (png_jmpbuf(png))) {
        png_destroy_read_struct (&png, &info, &end_info);
        fclose (file);
        delete [] row;
        return IMIO_READERROR;
    }

    for (unsigned int i = 0; i < height; i++) {

        png_read_row (png, (png_byte*)row, nullptr);

        if (bit_depth == 16) { // convert scanline to host byte order
            unsigned short* srow = (unsigned short*)row;

            for (unsigned int j = 0; j < width * 3; j++) {
                srow[j] = ntohs (srow[j]);
            }
        }

        setScanline (i, row, bit_depth);

        if (pl && !(i % 100)) {
            pl->setProgress ((double)(i + 1) / height);
        }
    }

    png_read_end (png, nullptr);
    png_destroy_read_struct (&png, &info, &end_info);

    delete [] row;
    fclose(file);

    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_READY");
        pl->setProgress (1.0);
    }

    return IMIO_SUCCESS;
}

// typedef struct  {
//     struct jpeg_error_mgr pub;  /* "public" fields */
//     jmp_buf setjmp_buffer;  /* for return to caller */
// } my_error_mgr;

// void my_error_exit (j_common_ptr cinfo)
// {
//     /* cinfo->err really points to a my_error_mgr struct, so coerce pointer */
//     my_error_mgr *myerr = (my_error_mgr*) cinfo->err;
//     /* Always display the message. */
//     /* We could postpone this until after returning, if we chose. */
//     (*cinfo->err->output_message) (cinfo);

//     /* Return control to the setjmp point */
// #if defined( WIN32 ) && defined( __x86_64__ ) && !defined(__clang__)
//     __builtin_longjmp(myerr->setjmp_buffer, 1);
// #else
//     longjmp(myerr->setjmp_buffer, 1);
// #endif
// }


int ImageIO::loadJPEGFromMemory(const char* buffer, int bufsize)
{
    jpeg_decompress_struct cinfo;
    jpeg_create_decompress(&cinfo);
    //rt_jpeg_memory_src (&cinfo, (const JOCTET*)buffer, bufsize);
    jpeg_mem_src(&cinfo, const_cast<unsigned char *>(reinterpret_cast<const unsigned char *>(buffer)), bufsize);

    /* We use our private extension JPEG error handler.
       Note that this struct must live as long as the main JPEG parameter
       struct, to avoid dangling-pointer problems.
    */
    //my_error_mgr jerr;
    rt_jpeg_error_mgr jerr;
    /* We set up the normal JPEG error routines, then override error_exit. */
    //cinfo.err = jpeg_std_error(&jerr.pub);
    //jerr.pub.error_exit = my_error_exit;
    cinfo.err = rt_jpeg_std_error(&jerr, "<MEMORY>", pl);

//     /* Establish the setjmp return context for my_error_exit to use. */
// #if defined( WIN32 ) && defined( __x86_64__ ) && !defined(__clang__)

//     if (__builtin_setjmp(jerr.setjmp_buffer)) {
// #else

//     if (setjmp(jerr.setjmp_buffer)) {
// #endif
//         /* If we get here, the JPEG code has signaled an error.
//            We need to clean up the JPEG object and return.
//         */
//         jpeg_destroy_decompress(&cinfo);
//         return IMIO_READERROR;
//     }

    try {
        if (pl) {
            pl->setProgressStr ("PROGRESSBAR_LOADJPEG");
            pl->setProgress (0.0);

        }

        setup_read_icc_profile (&cinfo);

        jpeg_read_header(&cinfo, TRUE);

        deleteLoadedProfileData();
        loadedProfileDataJpg = true;
        bool hasprofile = read_icc_profile (&cinfo, (JOCTET**)&loadedProfileData, (unsigned int*)&loadedProfileLength);

        if (hasprofile) {
            embProfile = cmsOpenProfileFromMem (loadedProfileData, loadedProfileLength);
        } else {
            embProfile = nullptr;
        }

        jpeg_start_decompress(&cinfo);

        unsigned int width = cinfo.output_width;
        unsigned int height = cinfo.output_height;

        allocate (width, height);

        std::vector<unsigned char> vrow(width * 3);
        unsigned char *row = &(vrow[0]); //new unsigned char[width * 3];

        while (cinfo.output_scanline < height) {
            if (jpeg_read_scanlines(&cinfo, &row, 1) < 1) {
                jpeg_finish_decompress(&cinfo);
                jpeg_destroy_decompress(&cinfo);
                //delete [] row;
                return IMIO_READERROR;
            }

            setScanline (cinfo.output_scanline - 1, row, 8, cinfo.num_components);

            if (pl && !(cinfo.output_scanline % 100)) {
                pl->setProgress ((double)(cinfo.output_scanline) / cinfo.output_height);
            }
        }

        //delete [] row;

        jpeg_finish_decompress(&cinfo);
        jpeg_destroy_decompress(&cinfo);

        if (pl) {
            pl->setProgressStr ("PROGRESSBAR_READY");
            pl->setProgress (1.0);
        }

        return IMIO_SUCCESS;
    } catch (rt_jpeg_error &) {
        jpeg_destroy_decompress(&cinfo);
        return IMIO_READERROR;
    }
}


int ImageIO::loadJPEG(const Glib::ustring &fname, int maxw_hint, int maxh_hint)
{
    FILE *file = g_fopen(fname.c_str (), "rb");

    if (!file) {
        return IMIO_CANNOTREADFILE;
    }

    jpeg_decompress_struct cinfo;
    rt_jpeg_error_mgr jerr;
    cinfo.err = rt_jpeg_std_error(&jerr, fname.c_str(), pl);
    jpeg_create_decompress(&cinfo);

    jpeg_stdio_src(&cinfo, file);

//    if ( setjmp((reinterpret_cast<rt_jpeg_error_mgr*>(cinfo.src))->error_jmp_buf) == 0 ) {
    try { 
        if (pl) {
            pl->setProgressStr ("PROGRESSBAR_LOADJPEG");
            pl->setProgress (0.0);
        }

        setup_read_icc_profile (&cinfo);

        //jpeg_stdio_src(&cinfo,file);
        jpeg_read_header(&cinfo, TRUE);

        //if JPEG is CMYK, then abort reading
        if (cinfo.jpeg_color_space == JCS_CMYK || cinfo.jpeg_color_space == JCS_YCCK) {
            jpeg_destroy_decompress(&cinfo);
            if (pl) {
                pl->error(M("JPEG_UNSUPPORTED_COLORSPACE_ERROR"));
            }
            return IMIO_READERROR;
        }

        cinfo.out_color_space = JCS_RGB;
        if (maxw_hint > 0 && maxh_hint > 0) {
            int w = cinfo.image_width;
            int h = cinfo.image_height;
            int d1 = w / maxw_hint;
            int d2 = h / maxh_hint;
            int d = std::min(d1, d2);
            if (d > 1) {
                cinfo.scale_num = 1;
                int l = std::min(d, 8);
                for (d = 1; (d << 1) <= l; d = d << 1) {}
                cinfo.scale_denom = d;
            }
        }

        deleteLoadedProfileData();
        loadedProfileDataJpg = true;
        bool hasprofile = read_icc_profile (&cinfo, (JOCTET**)&loadedProfileData, (unsigned int*)&loadedProfileLength);

        if (hasprofile) {
            embProfile = cmsOpenProfileFromMem (loadedProfileData, loadedProfileLength);
        } else {
            embProfile = nullptr;
        }

        jpeg_start_decompress(&cinfo);

        unsigned int width = cinfo.output_width;
        unsigned int height = cinfo.output_height;

        allocate(width, height);

        std::vector<unsigned char> vrow(width *3);
        unsigned char *row = &(vrow[0]);//new unsigned char[width * 3];

        while (cinfo.output_scanline < height) {
            if (jpeg_read_scanlines(&cinfo, &row, 1) < 1) {
                jpeg_finish_decompress(&cinfo);
                jpeg_destroy_decompress(&cinfo);
                //delete [] row;
                return IMIO_READERROR;
            }

            setScanline (cinfo.output_scanline - 1, row, 8);

            if (pl && !(cinfo.output_scanline % 100)) {
                pl->setProgress ((double)(cinfo.output_scanline) / cinfo.output_height);
            }
        }

        //delete [] row;

        jpeg_finish_decompress(&cinfo);
        jpeg_destroy_decompress(&cinfo);
        fclose(file);

        if (pl) {
            pl->setProgressStr ("PROGRESSBAR_READY");
            pl->setProgress (1.0);
        }

        return IMIO_SUCCESS;
    } catch (rt_jpeg_error &) {//else {
        jpeg_destroy_decompress(&cinfo);
        return IMIO_READERROR;
    }
}

int ImageIO::getTIFFSampleFormat (const Glib::ustring &fname, IIOSampleFormat &sFormat, IIOSampleArrangement &sArrangement)
{
#ifdef WIN32
    wchar_t *wfilename = (wchar_t*)g_utf8_to_utf16 (fname.c_str(), -1, NULL, NULL, NULL);
    TIFF* in = TIFFOpenW (wfilename, "r");
    g_free (wfilename);
#else
    TIFF* in = TIFFOpen(fname.c_str(), "r");
#endif

    if (in == nullptr) {
        return IMIO_CANNOTREADFILE;
    }

    uint16 bitspersample = 0, samplesperpixel = 0, sampleformat = 0;
    int hasTag = TIFFGetField(in, TIFFTAG_BITSPERSAMPLE, &bitspersample);
    hasTag &= TIFFGetField(in, TIFFTAG_SAMPLESPERPIXEL, &samplesperpixel);

    if (!hasTag) {
        // These are needed
        TIFFClose(in);
        sFormat = IIOSF_UNKNOWN;
        return IMIO_VARIANTNOTSUPPORTED;
    }

    if (!TIFFGetField(in, TIFFTAG_SAMPLEFORMAT, &sampleformat)) {
        /*
         * WARNING: This is a dirty hack!
         * We assume that files which doesn't contain the TIFFTAG_SAMPLEFORMAT tag
         * (which is the case with uncompressed TIFFs produced by RT!) are RGB files,
         * but that may be not true.   --- Hombre
         */
        sampleformat = SAMPLEFORMAT_UINT;
    } else if (sampleformat == SAMPLEFORMAT_VOID) {
        // according to https://www.awaresystems.be/imaging/tiff/tifftags/sampleformat.html
        // we assume SAMPLEFORMAT_UINT if SAMPLEFORMAT_VOID is set
        sampleformat = SAMPLEFORMAT_UINT;
    }

    uint16 config;
    TIFFGetField(in, TIFFTAG_PLANARCONFIG, &config);

    if (config == PLANARCONFIG_CONTIG) {
        sArrangement = IIOSA_CHUNKY;
    } else {
        sFormat = IIOSF_UNKNOWN;
        sArrangement = IIOSA_UNKNOWN;
        TIFFClose(in);
        return IMIO_VARIANTNOTSUPPORTED;
    }

    uint16 photometric;

    if (!TIFFGetField(in, TIFFTAG_PHOTOMETRIC, &photometric)) {
        TIFFClose(in);
        return IMIO_VARIANTNOTSUPPORTED;
    }

    uint16 compression;

    if (photometric == PHOTOMETRIC_LOGLUV)
        if (!TIFFGetField(in, TIFFTAG_COMPRESSION, &compression)) {
            compression = COMPRESSION_NONE;
        }

    TIFFClose(in);

    if (photometric == PHOTOMETRIC_RGB || photometric == PHOTOMETRIC_MINISBLACK) {
        if ((samplesperpixel == 1 || samplesperpixel == 3 || samplesperpixel == 4) && sampleformat == SAMPLEFORMAT_UINT) {
            if (bitspersample == 8) {
                sFormat = IIOSF_UNSIGNED_CHAR;
                return IMIO_SUCCESS;
            }

            if (bitspersample == 16) {
                sFormat = IIOSF_UNSIGNED_SHORT;
                return IMIO_SUCCESS;
            }
        } else if ((samplesperpixel == 3 || samplesperpixel == 4) && sampleformat == SAMPLEFORMAT_IEEEFP) {
            if (bitspersample==16) {
                sFormat = IIOSF_FLOAT16;
                return IMIO_SUCCESS;
            }
            if (bitspersample == 24) {
                sFormat = IIOSF_FLOAT24;
                return IMIO_SUCCESS;
            }
            if (bitspersample == 32) {
                sFormat = IIOSF_FLOAT32;
                return IMIO_SUCCESS;
            }
        }
    } else if ((samplesperpixel == 3 || samplesperpixel == 4) && photometric == PHOTOMETRIC_LOGLUV) {
        if (compression == COMPRESSION_SGILOG24) {
            sFormat = IIOSF_LOGLUV24;
            return IMIO_SUCCESS;
        } else if (compression == COMPRESSION_SGILOG) {
            sFormat = IIOSF_LOGLUV32;
            return IMIO_SUCCESS;
        }
    }

    return IMIO_VARIANTNOTSUPPORTED;
}


namespace {

tsize_t tiff_Read(thandle_t st, tdata_t buffer, tsize_t size)
{
    IMFILE *f = static_cast<IMFILE *>(st);
    return fread(buffer, 1, size, f);
}


tsize_t tiff_Write(thandle_t st, tdata_t buffer, tsize_t size)
{
    return 0;
}


int tiff_Close(thandle_t st)
{
    IMFILE *f = static_cast<IMFILE *>(st);
    fclose(f);
    return 0;
}


toff_t tiff_Seek(thandle_t st, toff_t pos, int whence)
{
    IMFILE *f = static_cast<IMFILE *>(st);
    fseek(f, pos, whence);
    return ftell(f);
}


toff_t tiff_Size(thandle_t st)
{
    IMFILE *f = static_cast<IMFILE *>(st);
    return f->size;
}


int tiff_Map(thandle_t, tdata_t*, toff_t*)
{
    return 0;
}

void tiff_Unmap(thandle_t, tdata_t, toff_t)
{
    return;
}

} // namespace


int ImageIO::loadTIFF(const Glib::ustring &fname)
{

    static MyMutex thumbMutex;
    MyMutex::MyLock lock(thumbMutex);

    if(!options.serializeTiffRead) {
        lock.release();
    }

    IMFILE *src = fopen(fname.c_str());
    TIFF *in = TIFFClientOpen(
        fname.c_str(), "r", static_cast<thandle_t>(src),
        tiff_Read, tiff_Write, tiff_Seek, tiff_Close, tiff_Size,
        tiff_Map, tiff_Unmap);

    if (in == nullptr) {
        return IMIO_CANNOTREADFILE;
    }

    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_LOADTIFF");
        pl->setProgress (0.0);
    }

    int width, height;
    TIFFGetField(in, TIFFTAG_IMAGEWIDTH, &width);
    TIFFGetField(in, TIFFTAG_IMAGELENGTH, &height);

    uint16 bitspersample, samplesperpixel;
    int hasTag = TIFFGetField(in, TIFFTAG_BITSPERSAMPLE, &bitspersample);
    hasTag &= TIFFGetField(in, TIFFTAG_SAMPLESPERPIXEL, &samplesperpixel);

    if (!hasTag) {
        // These are needed
        TIFFClose(in);
        return IMIO_VARIANTNOTSUPPORTED;
    }

    uint16 config;
    TIFFGetField(in, TIFFTAG_PLANARCONFIG, &config);

    if (config != PLANARCONFIG_CONTIG) {
        TIFFClose(in);
        return IMIO_VARIANTNOTSUPPORTED;
    }

    if (sampleFormat & (IIOSF_LOGLUV24 | IIOSF_LOGLUV32)) {
        TIFFSetField(in, TIFFTAG_SGILOGDATAFMT, SGILOGDATAFMT_FLOAT);
    }

    /*
     * We could use the min/max values set in TIFFTAG_SMINSAMPLEVALUE and
     * TIFFTAG_SMAXSAMPLEVALUE, but for now, we normalize the image to the
     * effective minimum and maximum values
     */
    if (options.rtSettings.verbose) {
        printf("Information of \"%s\":\n", fname.c_str());
        uint16 tiffDefaultScale, tiffBaselineExposure, tiffLinearResponseLimit;
        if (TIFFGetField(in, TIFFTAG_DEFAULTSCALE, &tiffDefaultScale)) {
            printf("   DefaultScale: %d\n", tiffDefaultScale);
        }
        else
            printf("   No DefaultScale value!\n");
        if (TIFFGetField(in, TIFFTAG_BASELINEEXPOSURE, &tiffBaselineExposure)) {
            printf("   BaselineExposure: %d\n", tiffBaselineExposure);
        }
        else
            printf("   No BaselineExposure value!\n");
        if (TIFFGetField(in, TIFFTAG_LINEARRESPONSELIMIT, &tiffLinearResponseLimit)) {
            printf("   LinearResponseLimit: %d\n", tiffLinearResponseLimit);
        }
        else
            printf("   No LinearResponseLimit value!\n");

        // agriggio 2020-11-06: this causes a segfault when compiled with -O3 and PROC_TARGET_NUMBER=2
        // on gcc 7.5.0, Ubuntu 18.04. Need to understand why...
        // 
        // uint16 tiffMinValue, tiffMaxValue;
        // if (TIFFGetField(in, TIFFTAG_SMINSAMPLEVALUE, &tiffMinValue)) {
        //     printf("   MinValue: %d\n", tiffMinValue);
        // }
        // else
        //     printf("   No minimum value!\n");
        // if (TIFFGetField(in, TIFFTAG_SMAXSAMPLEVALUE, &tiffMaxValue)) {
        //     printf("   MaxValue: %d\n\n", tiffMaxValue);
        // }
        // else
        //     printf("   No maximum value!\n\n");
        // printf("   Those values are not taken into account, the image data are normalized to a [0;1] range\n\n");
    }

    char* profdata;
    deleteLoadedProfileData();
    loadedProfileDataJpg = false;

    if (TIFFGetField(in, TIFFTAG_ICCPROFILE, &loadedProfileLength, &profdata)) {
        embProfile = cmsOpenProfileFromMem (profdata, loadedProfileLength);
        loadedProfileData = new char [loadedProfileLength];
        memcpy (loadedProfileData, profdata, loadedProfileLength);
    } else {
        embProfile = nullptr;
    }

    allocate (width, height);

    unsigned char* linebuffer = new unsigned char[TIFFScanlineSize(in) * (samplesperpixel == 1 ? 3 : 1)];

    for (int row = 0; row < height; row++) {
        if (TIFFReadScanline(in, linebuffer, row, 0) < 0) {
            TIFFClose(in);
            delete [] linebuffer;
            return IMIO_READERROR;
        }

        if (samplesperpixel > 3) {
            for (int i = 0; i < width; i++) {
                memmove(linebuffer + i * 3 * bitspersample / 8, linebuffer + i * samplesperpixel * bitspersample / 8, 3 * bitspersample / 8);
            }
        }
        else if (samplesperpixel == 1) {
            const size_t bytes = bitspersample / 8;
            for (int i = width - 1; i >= 0; --i) {
                const unsigned char* const src = linebuffer + i * bytes;
                unsigned char* const dest = linebuffer + i * 3 * bytes;
                memcpy(dest + 2 * bytes, src, bytes);
                memcpy(dest + 1 * bytes, src, bytes);
                memcpy(dest + 0 * bytes, src, bytes);
            }
        }

        setScanline (row, linebuffer, bitspersample);

        if (pl && !(row % 100)) {
            pl->setProgress ((double)(row + 1) / height);
        }
    }

    TIFFClose(in);
    delete [] linebuffer;

    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_READY");
        pl->setProgress (1.0);
    }

    return IMIO_SUCCESS;
}


int ImageIO::loadPPMFromMemory(const char* buffer, int width, int height, bool swap, int bps)
{
    allocate (width, height);

    int line_length(width * 3 * (bps / 8));

    if ( swap && bps > 8 ) {
        char swapped[line_length];

        for ( int row = 0; row < height; ++row ) {
            ::rtengine::swab(((char*)buffer) + (row * line_length), swapped, line_length);
            setScanline(row, (unsigned char*)&swapped[0], bps);
        }
    } else {
        for ( int row = 0; row < height; ++row ) {
            setScanline(row, ((unsigned char*)buffer) + (row * line_length), bps);
        }
    }

    return IMIO_SUCCESS;
}


int ImageIO::savePNG(const Glib::ustring &fname, int bps, bool uncompressed) const
{
    if (getWidth() < 1 || getHeight() < 1) {
        return IMIO_HEADERERROR;
    }

    FILE* const file = g_fopen_withBinaryAndLock (fname);

    if (!file) {
        return IMIO_CANNOTWRITEFILE;
    }

    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_SAVEPNG");
        pl->setProgress (0.0);
    }

    png_structp png = png_create_write_struct (PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);

    if (!png) {
        fclose (file);
        return IMIO_HEADERERROR;
    }

    // silence the warning about "invalid" sRGB profiles -- see #4260
#if defined(PNG_SKIP_sRGB_CHECK_PROFILE) && defined(PNG_SET_OPTION_SUPPORTED)
    png_set_option(png, PNG_SKIP_sRGB_CHECK_PROFILE, PNG_OPTION_ON);
#endif
    
    png_infop info = png_create_info_struct(png);

    if (!info) {
        png_destroy_write_struct (&png, nullptr);
        fclose (file);
        return IMIO_HEADERERROR;
    }

    if (setjmp(png_jmpbuf(png))) {
        png_destroy_write_struct (&png, &info);
        fclose(file);
        return IMIO_CANNOTWRITEFILE;
    }

    png_set_write_fn (png, file, png_write_data, png_flush);

    png_set_filter(png, 0, PNG_FILTER_PAETH);
    if (!uncompressed) {
        png_set_compression_level(png, 6);
    } else {
        png_set_compression_level(png, 0);
    }
    png_set_compression_strategy(png, 3);

    int width = getWidth ();
    int height = getHeight ();

    if (bps < 0) {
        bps = getBPS ();
    }
    if (bps > 16) {
        bps = 16;
    }

    png_set_IHDR(png, info, width, height, bps, PNG_COLOR_TYPE_RGB,
                 PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_BASE);

    if (profileData) {
#if PNG_LIBPNG_VER < 10500
        png_charp profdata = reinterpret_cast<png_charp>(profileData);
#else
        png_bytep profdata = reinterpret_cast<png_bytep>(profileData);
#endif
        png_set_iCCP(png, info, const_cast<png_charp>("icc"), 0, profdata, profileLength);
    }

    int rowlen = width * 3 * bps / 8;
    unsigned char *row = new unsigned char [rowlen];

    png_write_info(png, info);

    for (int i = 0; i < height; i++) {
        getScanline (i, row, bps);

        if (bps == 16) {
            // convert to network byte order
#if __BYTE_ORDER__==__ORDER_LITTLE_ENDIAN__
            for (int j = 0; j < width * 6; j += 2) {
                unsigned char tmp = row[j];
                row[j] = row[j + 1];
                row[j + 1] = tmp;
            }

#endif
        }

        png_write_row (png, (png_byte*)row);

        if (pl && !(i % 100)) {
            pl->setProgress ((double)(i + 1) / height);
        }
    }

    png_write_end(png, info);
    png_destroy_write_struct(&png, &info);

    delete [] row;
    fclose (file);

    if (!saveMetadata(fname)) {
        g_remove(fname.c_str());
        return IMIO_CANNOTWRITEFILE;
    }

    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_READY");
        pl->setProgress (1.0);
    }

    return IMIO_SUCCESS;
}



// Quality 0..100, subsampling: 1=low quality, 2=medium, 3=high
int ImageIO::saveJPEG (const Glib::ustring &fname, int quality, int subSamp) const
{
    if (getWidth() < 1 || getHeight() < 1) {
        return IMIO_HEADERERROR;
    }

    FILE* const file = g_fopen_withBinaryAndLock (fname);

    if (!file) {
        return IMIO_CANNOTWRITEFILE;
    }

    jpeg_compress_struct cinfo;
    /* We use our private extension JPEG error handler.
       Note that this struct must live as long as the main JPEG parameter
       struct, to avoid dangling-pointer problems.
    */
    //my_error_mgr jerr;
    /* We set up the normal JPEG error routines, then override error_exit. */
    //cinfo.err = jpeg_std_error(&jerr.pub);
    //jerr.pub.error_exit = my_error_exit;
    rt_jpeg_error_mgr jerr;
    cinfo.err = rt_jpeg_std_error(&jerr, fname.c_str(), pl);

//     /* Establish the setjmp return context for my_error_exit to use. */
// #if defined( WIN32 ) && defined( __x86_64__ ) && !defined(__clang__)

//     if (__builtin_setjmp(jerr.setjmp_buffer)) {
// #else

//     if (setjmp(jerr.setjmp_buffer)) {
// #endif
//         /* If we get here, the JPEG code has signaled an error.
//            We need to clean up the JPEG object, close the file, remove the already saved part of the file and return.
//         */
//         jpeg_destroy_compress(&cinfo);
//         fclose(file);
//         g_remove (fname.c_str());
//         return IMIO_CANNOTWRITEFILE;
//     }

    jpeg_create_compress (&cinfo);

    try {
        if (pl) {
            pl->setProgressStr ("PROGRESSBAR_SAVEJPEG");
            pl->setProgress (0.0);
        }

        jpeg_stdio_dest (&cinfo, file);

        int width = getWidth ();
        int height = getHeight ();

        cinfo.image_width  = width;
        cinfo.image_height = height;
        cinfo.in_color_space = JCS_RGB;
        cinfo.input_components = 3;
        jpeg_set_defaults (&cinfo);
        cinfo.write_JFIF_header = FALSE;

        // compute optimal Huffman coding tables for the image. Bit slower to generate, but size of result image is a bit less (default was FALSE)
        cinfo.optimize_coding = TRUE;

        // Since math coprocessors are common these days, FLOAT should be a bit more accurate AND fast (default is ISLOW)
        // (machine dependency is not really an issue, since we all run on x86 and having exactly the same file is not a requirement)
        cinfo.dct_method = JDCT_FLOAT;

        if (quality >= 0 && quality <= 100) {
            jpeg_set_quality (&cinfo, quality, true);
        }

        cinfo.comp_info[1].h_samp_factor = cinfo.comp_info[1].v_samp_factor = 1;
        cinfo.comp_info[2].h_samp_factor = cinfo.comp_info[2].v_samp_factor = 1;

        if (subSamp == 1) {
            // Best compression, default of the JPEG library:  2x2, 1x1, 1x1 (4:2:0)
            cinfo.comp_info[0].h_samp_factor = cinfo.comp_info[0].v_samp_factor = 2;
        } else if (subSamp == 2) {
            // Widely used normal ratio 2x1, 1x1, 1x1 (4:2:2)
            cinfo.comp_info[0].h_samp_factor = 2;
            cinfo.comp_info[0].v_samp_factor = 1;
        } else if (subSamp == 3) {
            // Best quality 1x1 1x1 1x1 (4:4:4)
            cinfo.comp_info[0].h_samp_factor = cinfo.comp_info[0].v_samp_factor = 1;
        }

        jpeg_start_compress(&cinfo, TRUE);

        // write icc profile to the output
        if (profileData) {
            write_icc_profile (&cinfo, (JOCTET*)profileData, profileLength);
        }

        // write image data
        int rowlen = width * 3;
        std::vector<unsigned char> vrow(rowlen);
        unsigned char *row = &(vrow[0]);//new unsigned char [rowlen];

//         /* To avoid memory leaks we establish a new setjmp return context for my_error_exit to use. */
// #if defined( WIN32 ) && defined( __x86_64__ ) && !defined(__clang__)

//         if (__builtin_setjmp(jerr.setjmp_buffer)) {
// #else

//             if (setjmp(jerr.setjmp_buffer)) {
// #endif
//                 /* If we get here, the JPEG code has signaled an error.
//                    We need to clean up the JPEG object, close the file, remove the already saved part of the file and return.
//                 */
//                 //delete [] row;
//                 jpeg_destroy_compress(&cinfo);
//                 fclose(file);
//                 g_remove (fname.c_str());
//                 return IMIO_CANNOTWRITEFILE;
//             }

        while (cinfo.next_scanline < cinfo.image_height) {

            getScanline (cinfo.next_scanline, row, 8);

            if (jpeg_write_scanlines (&cinfo, &row, 1) < 1) {
                jpeg_destroy_compress (&cinfo);
                //delete [] row;
                fclose (file);
                g_remove (fname.c_str());
                return IMIO_CANNOTWRITEFILE;
            }

            if (pl && !(cinfo.next_scanline % 100)) {
                pl->setProgress ((double)(cinfo.next_scanline) / cinfo.image_height);
            }
        }

        jpeg_finish_compress (&cinfo);
        jpeg_destroy_compress (&cinfo);

        //delete [] row;

        fclose (file);
    } catch (rt_jpeg_error &e) {
        jpeg_destroy_compress(&cinfo);
        fclose(file);
        g_remove(fname.c_str());
        return IMIO_CANNOTWRITEFILE;
    }

    if (!saveMetadata(fname)) {
        g_remove(fname.c_str());
        return IMIO_CANNOTWRITEFILE;
    }

    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_READY");
        pl->setProgress (1.0);
    }

    return IMIO_SUCCESS;
}


int ImageIO::saveTIFF (const Glib::ustring &fname, int bps, bool isFloat, bool uncompressed) const
{
    if (getWidth() < 1 || getHeight() < 1) {
        return IMIO_HEADERERROR;
    }

    bool writeOk = true;
    int width = getWidth ();
    int height = getHeight ();

    if (bps < 0) {
        bps = getBPS ();
    }

    int lineWidth = width * 3 * bps / 8;
    unsigned char* linebuffer = new unsigned char[lineWidth];

    // little hack to get libTiff to use proper byte order (see TIFFClienOpen()):
    const char *mode = "w";
#ifdef WIN32
    FILE *file = g_fopen_withBinaryAndLock (fname);
    int fileno = _fileno(file);
    int osfileno = _get_osfhandle(fileno);
    TIFF* out = TIFFFdOpen (osfileno, fname.c_str(), mode);
#else
    TIFF* out = TIFFOpen(fname.c_str(), mode);
    // int fileno = TIFFFileno (out);
#endif

    if (!out) {
        delete [] linebuffer;
        return IMIO_CANNOTWRITEFILE;
    }

    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_SAVETIFF");
        pl->setProgress (0.0);
    }

    bool needsReverse = false;

    TIFFSetField (out, TIFFTAG_SOFTWARE, RTNAME " " RTVERSION);
    TIFFSetField (out, TIFFTAG_IMAGEWIDTH, width);
    TIFFSetField (out, TIFFTAG_IMAGELENGTH, height);
    TIFFSetField (out, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
    TIFFSetField (out, TIFFTAG_SAMPLESPERPIXEL, 3);
    TIFFSetField (out, TIFFTAG_ROWSPERSTRIP, height);
    TIFFSetField (out, TIFFTAG_BITSPERSAMPLE, bps);
    TIFFSetField (out, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
    TIFFSetField (out, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);
    TIFFSetField (out, TIFFTAG_COMPRESSION, uncompressed ? COMPRESSION_NONE : COMPRESSION_ADOBE_DEFLATE);
    TIFFSetField (out, TIFFTAG_SAMPLEFORMAT, (bps == 16 || bps == 32) && isFloat ? SAMPLEFORMAT_IEEEFP : SAMPLEFORMAT_UINT);

    // somehow Exiv2 (tested with 0.27.3) doesn't seem to be able to update
    // XResolution and YResolution, so we do it ourselves here....
    constexpr float default_resolution = 300.f;
    float x_res = default_resolution;
    float y_res = default_resolution;
    int res_unit = RESUNIT_INCH;
    if (!metadataInfo.filename().empty()) {
        auto exif = metadataInfo.getOutputExifData();
        auto it = exif.findKey(Exiv2::ExifKey("Exif.Image.XResolution"));
        if (it != exif.end()) {
            x_res = it->toFloat();
        }
        it = exif.findKey(Exiv2::ExifKey("Exif.Image.YResolution"));
        if (it != exif.end()) {
            y_res = it->toFloat();
        }
        it = exif.findKey(Exiv2::ExifKey("Exif.Image.ResolutionUnit"));
        if (it != exif.end()) {
            res_unit = exiv2_to_long(*it);
        }
    }
    TIFFSetField(out, TIFFTAG_XRESOLUTION, x_res);
    TIFFSetField(out, TIFFTAG_YRESOLUTION, y_res);
    TIFFSetField(out, TIFFTAG_RESOLUTIONUNIT, res_unit);

    if (!uncompressed) {
        TIFFSetField (out, TIFFTAG_PREDICTOR, (bps == 16 || bps == 32) && isFloat ? PREDICTOR_FLOATINGPOINT : PREDICTOR_HORIZONTAL);
    }
    if (profileData) {
        TIFFSetField (out, TIFFTAG_ICCPROFILE, profileLength, profileData);
    }

    for (int row = 0; row < height; row++) {
        getScanline (row, linebuffer, bps, isFloat);

        if (bps == 16) {
            if(needsReverse && !uncompressed && isFloat) {
                for(int i = 0; i < lineWidth; i += 2) {
                    char temp = linebuffer[i];
                    linebuffer[i] = linebuffer[i + 1];
                    linebuffer[i + 1] = temp;
                }
            }
        } else if (bps == 32) {
            if(needsReverse && !uncompressed) {
                for(int i = 0; i < lineWidth; i += 4) {
                    char temp = linebuffer[i];
                    linebuffer[i] = linebuffer[i + 3];
                    linebuffer[i + 3] = temp;
                    temp = linebuffer[i + 1];
                    linebuffer[i + 1] = linebuffer[i + 2];
                    linebuffer[i + 2] = temp;
                }
            }
        }

        if (TIFFWriteScanline (out, linebuffer, row, 0) < 0) {
            TIFFClose (out);
            delete [] linebuffer;
            return IMIO_CANNOTWRITEFILE;
        }

        if (pl && !(row % 100)) {
            pl->setProgress ((double)(row + 1) / height);
        }
    }

    if (TIFFFlush(out) != 1) {
        writeOk = false;
    }

    TIFFClose (out);
#ifdef WIN32
    fclose (file);
#endif

    delete [] linebuffer;

    if (!saveMetadata(fname)) {
        writeOk = false;
    }

    if (pl) {
        pl->setProgressStr ("PROGRESSBAR_READY");
        pl->setProgress (1.0);
    }

    if(writeOk) {
        return IMIO_SUCCESS;
    } else {
        g_remove (fname.c_str());
        return IMIO_CANNOTWRITEFILE;
    }
}

// PNG read and write routines:

void png_read_data(png_structp png_ptr, png_bytep data, png_size_t length)
{
    png_size_t check;

    /* fread() returns 0 on error, so it is OK to store this in a png_size_t
     * instead of an int, which is what fread() actually returns.
     */
    check = (png_size_t)fread(data, (png_size_t)1, length, (FILE *)png_get_io_ptr(png_ptr));

    if (check != length) {
        png_error(png_ptr, "Read Error");
    }
}

void png_write_data(png_structp png_ptr, png_bytep data, png_size_t length)
{
    png_uint_32 check;

    check = fwrite(data, 1, length, (FILE *)png_get_io_ptr(png_ptr));

    if (check != length) {
        png_error(png_ptr, "Write Error");
    }
}

void png_flush(png_structp png_ptr)
{
    FILE *io_ptr;
    io_ptr = (FILE *)(png_get_io_ptr(png_ptr));

    if (io_ptr != nullptr) {
        fflush(io_ptr);
    }
}


int ImageIO::load(const Glib::ustring &fname, int maxw_hint, int maxh_hint)
{

    if (hasPngExtension(fname)) {
        return loadPNG(fname);
    } else if (hasJpegExtension(fname)) {
        return loadJPEG(fname, maxw_hint, maxh_hint);
    } else if (hasTiffExtension(fname)) {
        return loadTIFF(fname);
    } else {
        return IMIO_FILETYPENOTSUPPORTED;
    }
}


int ImageIO::save(const Glib::ustring &fname) const
{
    if (hasPngExtension(fname)) {
        return savePNG (fname);
    } else if (hasJpegExtension(fname)) {
        return saveJPEG (fname);
    } else if (hasTiffExtension(fname)) {
        return saveTIFF (fname);
    } else {
        return IMIO_FILETYPENOTSUPPORTED;
    }
}

void ImageIO::setProgressListener (ProgressListener* l)
{
    pl = l;
}

void ImageIO::setSampleFormat(IIOSampleFormat sFormat)
{
    sampleFormat = sFormat;
}

IIOSampleFormat ImageIO::getSampleFormat() const
{
    return sampleFormat;
}

void ImageIO::setSampleArrangement(IIOSampleArrangement sArrangement)
{
    sampleArrangement = sArrangement;
}

IIOSampleArrangement ImageIO::getSampleArrangement() const
{
    return sampleArrangement;
}

cmsHPROFILE ImageIO::getEmbeddedProfile () const
{
    return embProfile;
}

void ImageIO::getEmbeddedProfileData (int& length, unsigned char*& pdata) const
{
    length = loadedProfileLength;
    pdata = (unsigned char*)loadedProfileData;
}

MyMutex& ImageIO::mutex ()
{
    return imutex;
}

void ImageIO::deleteLoadedProfileData( )
{
    if(loadedProfileData) {
        if(loadedProfileDataJpg) {
            free(loadedProfileData);
        } else {
            delete[] loadedProfileData;
        }
    }

    loadedProfileData = nullptr;
}


bool ImageIO::saveMetadata(const Glib::ustring &fname) const
{
    if (metadataInfo.filename().empty()) {
        return true;
    }

    bool has_meta = true;
    try {
        metadataInfo.load();
    } catch (std::exception &exc) {
        // if (settings->verbose) {
        //     std::cout << "EXIF LOAD ERROR: " << exc.what() << std::endl;
        // }
        if (pl) {
            pl->error(Glib::ustring::compose(M("METADATA_LOAD_ERROR"), metadataInfo.filename(), exc.what()));
        }
        has_meta = false;
    }

    if (has_meta) {
        try {
            metadataInfo.saveToImage(pl, fname, false);
            if (!profileData) {
                Exiv2Metadata outmd(fname);
                outmd.exifData()["Exif.Photo.ColorSpace"] = 1;
                outmd.saveToImage(nullptr, fname, true);
            }
        } catch (std::exception &exc) {
            //std::cout << "EXIF ERROR: " << exc.what() << std::endl;
            //return false;
            if (pl) {
                pl->error(Glib::ustring::compose(M("METADATA_SAVE_ERROR"), fname, exc.what()));
            }
        }
    }

    return true;
}
