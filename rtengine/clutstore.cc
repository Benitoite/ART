#include <algorithm>
#include <unordered_map>

#include "clutstore.h"

#include "iccstore.h"
#include "imagefloat.h"
#include "opthelper.h"
#include "rt_math.h"
#include "stdimagesource.h"
#include "linalgebra.h"
#include "settings.h"

#include <giomm.h>
#include <sstream>
#include <iostream>
#include <fstream>

#include "../rtgui/options.h"
#include "../rtgui/multilangmgr.h"
#include "cJSON.h"

#ifdef _OPENMP
# include <omp.h>
#endif

namespace rtengine {

extern const Settings *settings;

namespace {

bool loadFile(
    const Glib::ustring& filename,
    const Glib::ustring& working_color_space,
    AlignedBuffer<std::uint16_t>& clut_image,
    unsigned int& clut_level
)
{
    rtengine::StdImageSource img_src;

    if (!Glib::file_test(filename, Glib::FILE_TEST_EXISTS) || img_src.load(filename)) {
        return false;
    }

    int fw, fh;
    img_src.getFullSize(fw, fh, TR_NONE);

    bool res = false;

    if (fw == fh) {
        int level = 1;

        while (level * level * level < fw) {
            ++level;
        }

        if (level * level * level == fw && level > 1) {
            clut_level = level;
            res = true;
        }
    }

    if (res) {
        rtengine::ColorTemp curr_wb = img_src.getWB();
        std::unique_ptr<rtengine::Imagefloat> img_float = std::unique_ptr<rtengine::Imagefloat>(new rtengine::Imagefloat(fw, fh));
        const PreviewProps pp(0, 0, fw, fh, 1);

        rtengine::procparams::ColorManagementParams icm;
        icm.workingProfile = working_color_space;

        img_src.getImage(curr_wb, TR_NONE, img_float.get(), pp, rtengine::procparams::ExposureParams(), rtengine::procparams::RAWParams());

        if (!working_color_space.empty()) {
            img_src.convertColorSpace(img_float.get(), icm, curr_wb);
        }

        AlignedBuffer<std::uint16_t> image(fw * fh * 4 + 4); // getClutValues() loads one pixel in advance

        std::size_t index = 0;

        for (int y = 0; y < fh; ++y) {
            for (int x = 0; x < fw; ++x) {
                image.data[index] = img_float->r(y, x);
                ++index;
                image.data[index] = img_float->g(y, x);
                ++index;
                image.data[index] = img_float->b(y, x);
                index += 2;
            }
        }

        clut_image.swap(image);
    }

    return res;
}

#ifdef __SSE2__
vfloat2 getClutValues(const AlignedBuffer<std::uint16_t>& clut_image, size_t index)
{
    const vint v_values = _mm_loadu_si128(reinterpret_cast<const vint*>(clut_image.data + index));
#ifdef __SSE4_1__
    return {
        _mm_cvtepi32_ps(_mm_cvtepu16_epi32(v_values)),
        _mm_cvtepi32_ps(_mm_cvtepu16_epi32(_mm_srli_si128(v_values, 8)))
    };
#else
    const vint v_mask = _mm_set1_epi32(0x0000FFFF);

    vint v_low = _mm_shuffle_epi32(v_values, _MM_SHUFFLE(1, 0, 1, 0));
    vint v_high = _mm_shuffle_epi32(v_values, _MM_SHUFFLE(3, 2, 3, 2));
    v_low = _mm_shufflelo_epi16(v_low, _MM_SHUFFLE(1, 1, 0, 0));
    v_high = _mm_shufflelo_epi16(v_high, _MM_SHUFFLE(1, 1, 0, 0));
    v_low = _mm_shufflehi_epi16(v_low, _MM_SHUFFLE(3, 3, 2, 2));
    v_high = _mm_shufflehi_epi16(v_high, _MM_SHUFFLE(3, 3, 2, 2));
    v_low = vandm(v_low, v_mask);
    v_high = vandm(v_high, v_mask);

    return {
        _mm_cvtepi32_ps(v_low),
        _mm_cvtepi32_ps(v_high)
    };
#endif
}
#endif

constexpr int TS = 112;

} // namespace

rtengine::HaldCLUT::HaldCLUT() :
    clut_level(0),
    flevel_minus_one(0.0f),
    flevel_minus_two(0.0f),
    clut_profile("sRGB")
{
}

rtengine::HaldCLUT::~HaldCLUT()
{
}

bool rtengine::HaldCLUT::load(const Glib::ustring& filename)
{
    if (loadFile(filename, "", clut_image, clut_level)) {
        Glib::ustring name, ext;
        rtengine::CLUTStore::splitClutFilename(filename, name, ext, clut_profile);

        clut_filename = filename;
        clut_level *= clut_level;
        flevel_minus_one = static_cast<float>(clut_level - 1) / 65535.0f;
        flevel_minus_two = static_cast<float>(clut_level - 2);
        return true;
    }

    return false;
}

rtengine::HaldCLUT::operator bool() const
{
    return !clut_image.isEmpty();
}

Glib::ustring rtengine::HaldCLUT::getFilename() const
{
    return clut_filename;
}

Glib::ustring rtengine::HaldCLUT::getProfile() const
{
    return clut_profile;
}

void rtengine::HaldCLUT::getRGB(
    float strength,
    std::size_t line_size,
    const float* r,
    const float* g,
    const float* b,
    float* out_rgbx
) const
{
    const unsigned int level = clut_level; // This is important

    const unsigned int level_square = level * level;

#ifdef __SSE2__
    const vfloat v_strength = F2V(strength);
#endif

    for (std::size_t column = 0; column < line_size; ++column, ++r, ++g, ++b, out_rgbx += 4) {
        const unsigned int red = std::min(flevel_minus_two, *r * flevel_minus_one);
        const unsigned int green = std::min(flevel_minus_two, *g * flevel_minus_one);
        const unsigned int blue = std::min(flevel_minus_two, *b * flevel_minus_one);

        const unsigned int color = red + green * level + blue * level_square;

#ifndef __SSE2__
        const float re = *r * flevel_minus_one - red;
        const float gr = *g * flevel_minus_one - green;
        const float bl = *b * flevel_minus_one - blue;

        size_t index = color * 4;

        float tmp1[4] ALIGNED16;
        tmp1[0] = intp<float>(re, clut_image.data[index + 4], clut_image.data[index]);
        tmp1[1] = intp<float>(re, clut_image.data[index + 5], clut_image.data[index + 1]);
        tmp1[2] = intp<float>(re, clut_image.data[index + 6], clut_image.data[index + 2]);

        index = (color + level) * 4;

        float tmp2[4] ALIGNED16;
        tmp2[0] = intp<float>(re, clut_image.data[index + 4], clut_image.data[index]);
        tmp2[1] = intp<float>(re, clut_image.data[index + 5], clut_image.data[index + 1]);
        tmp2[2] = intp<float>(re, clut_image.data[index + 6], clut_image.data[index + 2]);

        out_rgbx[0] = intp<float>(gr, tmp2[0], tmp1[0]);
        out_rgbx[1] = intp<float>(gr, tmp2[1], tmp1[1]);
        out_rgbx[2] = intp<float>(gr, tmp2[2], tmp1[2]);

        index = (color + level_square) * 4;

        tmp1[0] = intp<float>(re, clut_image.data[index + 4], clut_image.data[index]);
        tmp1[1] = intp<float>(re, clut_image.data[index + 5], clut_image.data[index + 1]);
        tmp1[2] = intp<float>(re, clut_image.data[index + 6], clut_image.data[index + 2]);

        index = (color + level + level_square) * 4;

        tmp2[0] = intp<float>(re, clut_image.data[index + 4], clut_image.data[index]);
        tmp2[1] = intp<float>(re, clut_image.data[index + 5], clut_image.data[index + 1]);
        tmp2[2] = intp<float>(re, clut_image.data[index + 6], clut_image.data[index + 2]);

        tmp1[0] = intp<float>(gr, tmp2[0], tmp1[0]);
        tmp1[1] = intp<float>(gr, tmp2[1], tmp1[1]);
        tmp1[2] = intp<float>(gr, tmp2[2], tmp1[2]);

        out_rgbx[0] = intp<float>(bl, tmp1[0], out_rgbx[0]);
        out_rgbx[1] = intp<float>(bl, tmp1[1], out_rgbx[1]);
        out_rgbx[2] = intp<float>(bl, tmp1[2], out_rgbx[2]);

        out_rgbx[0] = intp<float>(strength, out_rgbx[0], *r);
        out_rgbx[1] = intp<float>(strength, out_rgbx[1], *g);
        out_rgbx[2] = intp<float>(strength, out_rgbx[2], *b);
#else
        const vfloat v_in = _mm_set_ps(0.0f, *b, *g, *r);
        const vfloat v_tmp = v_in * F2V(flevel_minus_one);
        const vfloat v_rgb = v_tmp - _mm_cvtepi32_ps(_mm_cvttps_epi32(vminf(v_tmp, F2V(flevel_minus_two))));

        size_t index = color * 4;

        const vfloat v_r = PERMUTEPS(v_rgb, _MM_SHUFFLE(0, 0, 0, 0));

        vfloat2 v_clut_values = getClutValues(clut_image, index);
        vfloat v_tmp1 = vintpf(v_r, v_clut_values.y, v_clut_values.x);

        index = (color + level) * 4;

        v_clut_values = getClutValues(clut_image, index);
        vfloat v_tmp2 = vintpf(v_r, v_clut_values.y, v_clut_values.x);

        const vfloat v_g = PERMUTEPS(v_rgb, _MM_SHUFFLE(1, 1, 1, 1));

        vfloat v_out = vintpf(v_g, v_tmp2, v_tmp1);

        index = (color + level_square) * 4;

        v_clut_values = getClutValues(clut_image, index);
        v_tmp1 = vintpf(v_r, v_clut_values.y, v_clut_values.x);

        index = (color + level + level_square) * 4;

        v_clut_values = getClutValues(clut_image, index);
        v_tmp2 = vintpf(v_r, v_clut_values.y, v_clut_values.x);

        v_tmp1 = vintpf(v_g, v_tmp2, v_tmp1);

        const vfloat v_b = PERMUTEPS(v_rgb, _MM_SHUFFLE(2, 2, 2, 2));

        v_out = vintpf(v_b, v_tmp1, v_out);

        STVF(*out_rgbx, vintpf(v_strength, v_out, v_in));
#endif
    }
}


Glib::ustring rtengine::CLUTStore::getClutDisplayName(const Glib::ustring &filename)
{
    Glib::ustring name;
    
#ifdef ART_USE_CTL
    if (getFileExtension(filename) == "ctl") {
        const Glib::ustring full_filename =
            !Glib::path_is_absolute(filename)
            ? Glib::ustring(Glib::build_filename(options.clutsDir, filename))
            : filename;
        if (Glib::file_test(filename, Glib::FILE_TEST_EXISTS)) {
            auto fn = Glib::filename_from_utf8(filename);
            std::ifstream src(fn.c_str());
            std::string line;
            bool found = false;
            while (src && std::getline(src, line) && !found) {
                size_t s = 0;
                while (s < line.size() && std::isspace(line[s])) {
                    ++s;
                }
                if (s+1 < line.size() && line[s] == '/' && line[s+1] == '/') {
                    s += 2;
                }
                while (s < line.size() && std::isspace(line[s])) {
                    ++s;
                }
                if (line.find("@ART-label:", s) == s) {
                    line = line.substr(11+s);
                    cJSON *root = cJSON_Parse(line.c_str());
                    if (root) {
                        if (cJSON_IsString(root)) {
                            name = cJSON_GetStringValue(root);
                            if (!name.empty() && name[0] == '$') {
                                name = M(name.c_str()+1);
                            }
                            found = !name.empty();
                        }
                        cJSON_Delete(root);
                    }
                    break;
                }
            }
            if (found) {
                return name;
            }
        }
    }
#endif
    
    Glib::ustring dummy;
    splitClutFilename(filename, name, dummy, dummy);
    return name;
}


void rtengine::CLUTStore::splitClutFilename(
    const Glib::ustring& filename,
    Glib::ustring& name,
    Glib::ustring& extension,
    Glib::ustring& profile_name
)
{
    Glib::ustring basename = Glib::path_get_basename(filename);

    const Glib::ustring::size_type last_dot_pos = basename.rfind('.');

    if (last_dot_pos != Glib::ustring::npos) {
        name.assign(basename, 0, last_dot_pos);
        extension.assign(basename, last_dot_pos + 1, Glib::ustring::npos);
    } else {
        name = basename;
    }

    profile_name = "sRGB";

    bool search_profile_name = true;
#ifdef ART_USE_OCIO
    search_profile_name = !(extension.casefold().find("clf") == 0);
#endif // ART_USE_OCIO
#ifdef ART_USE_CTL
    search_profile_name = search_profile_name && !(extension.casefold().find("ctl") == 0);
#endif // ART_USE_CTL
    if (search_profile_name && !name.empty()) {
        for (const auto& working_profile : rtengine::ICCStore::getInstance()->getWorkingProfiles()) {
            if (
                !working_profile.empty() // This isn't strictly needed, but an empty wp name should be skipped anyway
                && std::search(name.rbegin(), name.rend(), working_profile.rbegin(), working_profile.rend()) == name.rbegin()
            ) {
                profile_name = working_profile;
                name.erase(name.size() - working_profile.size());
                break;
            }
        }
    } else if (!search_profile_name) {
        profile_name = "";
    }
}

rtengine::CLUTStore& rtengine::CLUTStore::getInstance()
{
    static CLUTStore instance;
    return instance;
}

std::shared_ptr<rtengine::HaldCLUT> rtengine::CLUTStore::getHaldClut(const Glib::ustring& filename) const
{
    MyMutex::MyLock lock(mutex_);
    std::shared_ptr<rtengine::HaldCLUT> result;

    const Glib::ustring full_filename =
        !Glib::path_is_absolute(filename)
            ? Glib::ustring(Glib::build_filename(options.clutsDir, filename))
            : filename;

    if (!cache.get(full_filename, result)) {
        std::unique_ptr<rtengine::HaldCLUT> clut(new rtengine::HaldCLUT);

        if (clut->load(full_filename)) {
            result = std::move(clut);
            cache.insert(full_filename, result);
        }
    }

    return result;
}


#ifdef ART_USE_OCIO

namespace {

std::string decompress_to_temp(const Glib::ustring &fname)
{
    std::string templ = Glib::build_filename(Glib::get_tmp_dir(), Glib::ustring::compose("ART-ocio-clf-%1-XXXXXX", Glib::path_get_basename(fname)));
    int fd = Glib::mkstemp(templ);
    if (fd < 0) {
        throw "error";
    }

    close(fd);
    bool err = false;

    auto f = Gio::File::create_for_path(templ);
    auto s = f->append_to();
    auto c = Gio::ZlibDecompressor::create(Gio::ZLIB_COMPRESSOR_FORMAT_GZIP);
    {
        auto stream = Gio::ConverterOutputStream::create(s, c);
        stream->set_close_base_stream(true);
        constexpr gsize chunk = 512;
        char buffer[chunk];
        auto src_fname = Glib::filename_from_utf8(fname);
        FILE *src = g_fopen(src_fname.c_str(), "rb");
        if (!src) {
            err = true; 
        } else {
            try {
                size_t n = 0;
                while ((n = fread(buffer, 1, chunk, src)) > 0) {
                    stream->write(buffer, n);
                }
            } catch (...) {
                err = true;
            }
            fclose(src);
        }
    }

    if (err) {
        g_remove(templ.c_str());
        throw "error";
    }
    return templ;
}


std::string copy_to_temp(const Glib::ustring &fname)
{
    std::string templ = Glib::build_filename(Glib::get_tmp_dir(), Glib::ustring::compose("ART-ocio-clf-%1-XXXXXX", Glib::path_get_basename(fname)));
    int fd = Glib::mkstemp(templ);
    if (fd < 0) {
        throw "error";
    }

    close(fd);
    bool err = false;

    auto f = Gio::File::create_for_path(templ);
    {
        auto stream = f->append_to();
        constexpr gsize chunk = 512;
        char buffer[chunk];
        auto src_fname = Glib::filename_from_utf8(fname);
        FILE *src = g_fopen(src_fname.c_str(), "rb");
        if (!src) {
            err = true; 
        } else {
            try {
                size_t n = 0;
                while ((n = fread(buffer, 1, chunk, src)) > 0) {
                    stream->write(buffer, n);
                }
            } catch (...) {
                err = true;
            }
            fclose(src);
        }
    }

    if (err) {
        g_remove(templ.c_str());
        throw "error";
    }
    return templ;
}

} // namespace

OCIO::ConstProcessorRcPtr rtengine::CLUTStore::getOCIOLut(const Glib::ustring& filename) const
{
    MyMutex::MyLock lock(mutex_);

    
    OCIOCacheEntry result;
    OCIO::ConstProcessorRcPtr retval;

    const Glib::ustring full_filename =
        !Glib::path_is_absolute(filename)
            ? Glib::ustring(Glib::build_filename(options.clutsDir, filename))
            : filename;

    auto ext = getFileExtension(full_filename);
    if (ext != "clf" && ext != "clfz") {
        return retval;
    }
    
    const auto md5 = getMD5(full_filename, true);

    bool found = ocio_cache_.get(full_filename, result);
    if (!found || result.second != md5) {
        std::string fn = "";
        bool del_fn = false;
        try {
            OCIO::ConstConfigRcPtr config = OCIO::Config::CreateRaw();
            OCIO::FileTransformRcPtr t = OCIO::FileTransform::Create();
            if (ext == "clfz") {
                fn = decompress_to_temp(full_filename);
                del_fn = true;
            } else {
                fn = copy_to_temp(full_filename);
                del_fn = true;
            }
            t->setSrc(fn.c_str());
            t->setInterpolation(OCIO::INTERP_BEST);
            retval = config->getProcessor(t);
            result = std::make_pair(retval, md5);
            ocio_cache_.set(full_filename, result);
        } catch (...) {
        }
        if (del_fn && !fn.empty()) {
            g_remove(fn.c_str());
        }
    } else {
        retval = result.first;
    }

    return retval;
}

#endif // ART_USE_OCIO

#ifdef ART_USE_CTL

namespace {

bool fill_from_json(std::unordered_map<std::string, int> &name2pos, std::vector<CLUTParamDescriptor> &params, cJSON *root)
{
    if (!cJSON_IsArray(root)) {
        return false;
    }
    
    auto sz = cJSON_GetArraySize(root);
    if (sz < 2) {
        return false;
    }
    
    auto n = cJSON_GetArrayItem(root, 0);
    if (!cJSON_IsString(n)) {
        return false;
    }
    std::string name = cJSON_GetStringValue(n);
    auto it = name2pos.find(name);
    if (it == name2pos.end()) {
        return false;
    }
    auto &desc = params[it->second];

    name2pos.erase(it);

    n = cJSON_GetArrayItem(root, 1);
    if (!cJSON_IsString(n)) {
        return false;
    }
    desc.gui_name = cJSON_GetStringValue(n);
    desc.gui_group = "";
    desc.gui_step = 1;

    const auto set_group =
        [&](int i) -> bool
        {
            auto n = cJSON_GetArrayItem(root, i);
            if (!cJSON_IsString(n)) {
                return false;
            }
            desc.gui_group = cJSON_GetStringValue(n);
            return true;
        };

    switch (desc.type) {
    case CLUTParamType::PT_BOOL:
        if (sz == 2) {
            return true;
        } else if (sz == 3 || sz == 4) {
            n = cJSON_GetArrayItem(root, 2);
            if (cJSON_IsBool(n)) {
                desc.value_default = cJSON_IsTrue(n);
            }
            if (sz == 4) {
                return set_group(3);
            } else {
                return true;
            }
        }
        break;
    case CLUTParamType::PT_FLOAT:
        if (sz >= 4 && sz <= 7) {
            n = cJSON_GetArrayItem(root, 2);
            if (cJSON_IsNumber(n)) {
                desc.value_min = n->valuedouble;
            } else {
                return false;
            }
            n = cJSON_GetArrayItem(root, 3);
            if (cJSON_IsNumber(n)) {
                desc.value_max = n->valuedouble;
            } else {
                return false;
            }
            if (sz >= 5) {
                n = cJSON_GetArrayItem(root, 4);
                if (cJSON_IsNumber(n)) {
                    desc.value_default = n->valuedouble;
                } else {
                    return false;
                }
                if (sz >= 6) {
                    n = cJSON_GetArrayItem(root, 5);
                    if (cJSON_IsNumber(n)) {
                        desc.gui_step = n->valuedouble;
                    } else {
                        return false;
                    }
                } else {
                    desc.gui_step = (desc.value_max - desc.value_min) / 100.0;
                }
                if (sz == 7) {
                    return set_group(6);
                }
            }
            return true;
        }
        break;
    case CLUTParamType::PT_INT:
        if (sz >= 3 && sz <= 6) {
            n = cJSON_GetArrayItem(root, 2);
            if (cJSON_IsArray(n)) {
                for (int i = 0, k = cJSON_GetArraySize(n); i < k; ++i) {
                    auto v = cJSON_GetArrayItem(n, i);
                    if (!cJSON_IsString(v)) {
                        return false;
                    }
                    desc.choices.push_back(cJSON_GetStringValue(v));
                }
                desc.type = CLUTParamType::PT_CHOICE;
                if (sz == 4) {
                    return set_group(3);
                } else {
                    return (sz == 3);
                }
            } else if (sz >= 4) {
                if (cJSON_IsNumber(n)) {
                    desc.value_min = n->valuedouble;
                } else {
                    return false;
                }
                n = cJSON_GetArrayItem(root, 3);
                if (cJSON_IsNumber(n)) {
                    desc.value_max = n->valuedouble;
                } else {
                    return false;
                }
                if (sz >= 5) {
                    n = cJSON_GetArrayItem(root, 4);
                    if (cJSON_IsNumber(n)) {
                        desc.value_default = n->valuedouble;
                    } else {
                        return false;
                    }
                    if (sz == 6) {
                        return set_group(5);
                    }
                }
                return true;
            } else {
                return false;
            }
        }
        break;
    default:
        return false;
    }

    return false;
}

/**
 * ART-compatible CTL scripts can contain parameters as additional uniform
 * input parameters to ART_main. Currently, only parameters of type "bool",
 * "int" or "float" are supported. Each such parameter must come with an
 * associated ART parameter definition in the CTL script. ART parameter
 * definitions are special comment lines of the following form:
 *
 *  // @ART-param: <param-def>
 *
 * where <param-def> is an array in JSON format, whose content depends on the
 * parameter type. The array must be at least of size 2; the first element is
 * a string containing the name of the parameter (which must mach the name
 * used in ART_main), and the second element is its GUI label. The rest of the
 * array has the following structure:
 *
 * - for "bool" parameters, the 3rd optional element specifies the default
 *   value; the 4th optional element instead is a "group name" for the GUI: if
 *   set, this will cause the control to appear under a collapsible panel with
 *   the given name in the GUI;
 *
 * - for "float" parameters, the array size must be at least 4 and at most 7.
 *   The 3rd and 4th elements are the minimum and maximum values for the
 *   GUI slider. The optional 5th element is the default value, the optional
 *   6th element the precision to use in the GUI (e.g. 0.01 will use 2 decimal
 *   digits in the GUI), and the optional last element is the GUI group name;
 *
 * - for "int" parameters, the array size must be at least 3 and at most 6.
 *   If the 3rd parameter is an array of strings, it is interpreted as a list
 *   of options in a choice menu, with values corresponding to their index in
 *   the array (i.e. the 1st option will give a value of 0, the 2nd a value of
 *   1, etc.). In this case, the array can contain at most 2 other elements,
 *   which are respectively the default value and the optional GUI group name.
 *   If the 3rd parameter is not an array of strings, then the array size must
 *   be at least 4, with the 3rd and 4th elements corresponding to the minimum
 *   and maximum values for the GUI slider. The optional 5th element is the
 *   default value, and the optional last element the GUI group name.
 *
 * If default values are not given in the ART parameter definition, they are
 * taken from the definition of the ART_main function. If no default is given,
 * zero is used.
 *
 * Example:
 *
 * // @ART-param: ["param_float", "A float slider", -1.0, 1.0, 0.5, 0.1]
 * // @ART-param: ["param_int", "An int slider", -10, 10]
 * // @ART-param: ["param_bool", "A checkbox", true]
 * // @ART-param: ["param_choice", "A combo box", ["Option A", "Option B"], 1, "Select between option A (value 0) and option B (value 1)"]
 * 
 * void ART_main(varying float r, varying float g, varying float b,
 *               output varying float or, output varying float og, output varying float ob,
 *               float param_float,
 *               int param_int,
 *               bool param_bool,
 *               int param_choice)
 * {
 *    // ...
 * }
 */ 
bool get_CTL_params(const Glib::ustring &filename, std::shared_ptr<Ctl::Interpreter> intp, Ctl::FunctionCallPtr func, std::vector<CLUTParamDescriptor> &out, Glib::ustring &colorspace)
{
    out.clear();
    std::unordered_map<std::string, int> name2pos;

    colorspace = "ACESp0";

    static std::unordered_map<std::string, std::string> profilemap = {
        {"aces2065-1", "ACESp0"},
        {"acescg", "ACESp1"},
        {"rec2020", "Rec2020"},
        {"prophoto", "ProPhoto"},
        {"rec709", "sRGB"},
        {"srgb", "sRGB"},
        {"adobergb", "Adobe RGB"},
        {"adobe", "Adobe RGB"}
    };

    const auto err =
        [&](const std::string &msg) -> bool
        {
            if (settings->verbose) {
                std::cout << "Error in CTL script from " << filename << ": "
                          << msg << std::endl;
            }
            return false;
        };
    
    for (size_t i = 3, n = func->numInputArgs(); i < n; ++i) {
        auto a = func->inputArg(i);
        if (a->isVarying()) {
            return err("varying parameter " + a->name());
        }
        CLUTParamType tp = CLUTParamType::PT_INT;
        switch (a->type()->cDataType()) {
        case Ctl::BoolTypeEnum:
            tp = CLUTParamType::PT_BOOL;
            break;
        case Ctl::IntTypeEnum:
            tp = a->type().cast<Ctl::BoolType>() ? CLUTParamType::PT_BOOL : CLUTParamType::PT_INT;
            break;
        case Ctl::FloatTypeEnum:
            tp = CLUTParamType::PT_FLOAT;
            break;
        default: return err("parameter " + a->name() + " is of unsupported type");
        }

        std::string name = a->name();
        name2pos[name] = out.size();
        
        out.emplace_back();
        auto &desc = out.back();
        desc.name = name;
        desc.type = tp;

        desc.value_min = 0;
        desc.value_max = 1;
        desc.value_default = 0;

        if (a->hasDefaultValue()) {
            switch (tp) {
            case CLUTParamType::PT_BOOL:
                desc.value_default = *reinterpret_cast<bool *>(a->data());
                break;
            case CLUTParamType::PT_FLOAT:
                desc.value_default = *reinterpret_cast<float *>(a->data());
                break;
            case CLUTParamType::PT_INT:
                desc.value_default = *reinterpret_cast<int *>(a->data());
            default:
                break;
            }
        }
    }

    if (Glib::file_test(filename, Glib::FILE_TEST_EXISTS)) {
        auto fn = Glib::filename_from_utf8(filename);
        std::ifstream src(fn.c_str());
        std::string line;
        while (src && std::getline(src, line)) {
            size_t s = 0;
            while (s < line.size() && std::isspace(line[s])) {
                ++s;
            }
            if (s+1 < line.size() && line[s] == '/' && line[s+1] == '/') {
                s += 2;
            }
            while (s < line.size() && std::isspace(line[s])) {
                ++s;
            }
            if (line.find("@ART-param:", s) == s) {
                line = line.substr(11+s);
                cJSON *root = cJSON_Parse(line.c_str());
                if (!root) {
                    return err("bad parameter definition: " + line);
                }
                bool ok = fill_from_json(name2pos, out, root);
                cJSON_Delete(root);
                if (!ok) {
                    return err("bad parameter definition: " + line);
                }
            } else if (line.find("@ART-colorspace:", s) == s) {
                line = line.substr(16+s);
                cJSON *root = cJSON_Parse(line.c_str());
                if (!root || !cJSON_IsString(root)) {
                    return err("invalid colorspace definition: " + line);
                }
                std::string name = Glib::ustring(cJSON_GetStringValue(root)).casefold();
                cJSON_Delete(root);
                auto it = profilemap.find(name);
                if (it != profilemap.end()) {
                    colorspace = it->second;
                } else {
                    return err("invalid colorspace definition: " + line);
                }
            }            
        }
    } else {
        return err("file reading error");
    }

    if (!name2pos.empty() && !out.empty()) {
        std::string msg = "missing parameter definitions: ";
        const char *sep = "";
        for (auto &p : name2pos) {
            msg += sep + p.first;
            sep = ", ";
        }
        return err(msg);
    }
    
    return true;
}

} // namespace

std::vector<Ctl::FunctionCallPtr> rtengine::CLUTStore::getCTLLut(const Glib::ustring& filename, int num_threads, int &chunk_size, std::vector<CLUTParamDescriptor> &params, Glib::ustring &colorspace) const
{
    MyMutex::MyLock lock(mutex_);
    
    CTLCacheEntry result;
    std::vector<Ctl::FunctionCallPtr> retval;
    std::shared_ptr<Ctl::Interpreter> intp;
    
    const Glib::ustring full_filename =
        !Glib::path_is_absolute(filename)
            ? Glib::ustring(Glib::build_filename(options.clutsDir, filename))
            : filename;
    if (!Glib::file_test(full_filename, Glib::FILE_TEST_IS_REGULAR) || getFileExtension(full_filename) != "ctl") {
        return retval;
    }
    const auto md5 = getMD5(full_filename, true);

    const auto err =
        [&](const char *msg) -> std::vector<Ctl::FunctionCallPtr>
        {
            if (settings->verbose) {
                std::cout << "Error in CTL script from " << full_filename << ": "
                          << msg << std::endl;
            }
            retval.clear();
            return retval;
        };

    try {
        bool found = ctl_cache_.get(full_filename, result);
        if (!found || result.md5 != md5) {
            intp = std::make_shared<Ctl::SimdInterpreter>();
            intp->setModulePaths({ Glib::path_get_dirname(full_filename) });
            intp->loadFile(full_filename);

            auto f = intp->newFunctionCall("ART_main");
            if (f->numInputArgs() < 3) {
                return err("wrong number of input arguments to ART_main");
            } else {
                for (int i = 0; i < 3; ++i) {
                    auto a = f->inputArg(i);
                    if (!a->type().cast<Ctl::FloatType>() || !a->isVarying()) {
                        return err("bad input arg type");
                    }
                }
            }
            if (f->numOutputArgs() != 3) {
                return err("wrong number of output arguments");
            } else {
                for (int i = 0; i < 3; ++i) {
                    auto a = f->outputArg(i);
                    if (!a->type().cast<Ctl::FloatType>() || !a->isVarying()) {
                        return err("bad output arg type");
                    }
                }
            }

            if (!get_CTL_params(full_filename, intp, f, params, colorspace)) {
                params.clear();
                return err("error in parsing CTL parameters");
            }

            result.intp = intp;
            result.md5 = md5;
            result.params = params;
            result.colorspace = colorspace;
            ctl_cache_.set(full_filename, result);
        } else {
            intp = result.intp;
            params = result.params;
            colorspace = result.colorspace;
        }
        if (intp) {
            for (int i = 0; i < num_threads; ++i) {
                retval.push_back(intp->newFunctionCall("ART_main"));
            }
            chunk_size = intp->maxSamples();
        }
    } catch (std::exception &exc) {
        return err(exc.what());
    }

    return retval;
}

#endif // ART_USE_CTL


void rtengine::CLUTStore::clearCache()
{
    MyMutex::MyLock lock(mutex_);
    
    cache.clear();
#ifdef ART_USE_OCIO
    ocio_cache_.clear();
#endif // ART_USE_OCIO
#ifdef ART_USE_CTL
    ctl_cache_.clear();
#endif // ART_USE_CTL
}


namespace {

inline float CTL_shaper_pq(float a, bool inv)
{
    constexpr float m1 = 2610.0 / 16384.0;
    constexpr float m2 = 2523.0 / 32.0;
    constexpr float c1 = 107.0 / 128.0;
    constexpr float c2 = 2413.0 / 128.0;
    constexpr float c3 = 2392.0 / 128.0;
    constexpr float scale = 100.0;

    if (a <= 0.f) {
        return 0.f;
    }
    
    if (!inv) {
        a /= scale;
        float aa = pow_F(a, m1);
        return pow_F((c1 + c2 * aa)/(1.f + c3 * aa), m2);
    } else {
        float p = pow_F(a, 1.f/m2);
        float aa = std::max(p - c1, 0.f) / (c2 - c3 * p);
        return pow_F(aa, 1.f/m1) * scale;
    }
}

} // namespace


rtengine::CLUTStore::CLUTStore() :
    cache(options.clutCacheSize)
#ifdef ART_USE_OCIO
    , ocio_cache_(options.clutCacheSize)
#endif // ART_USE_OCIO
#ifdef ART_USE_CTL
    , ctl_cache_(options.clutCacheSize)
#endif // ART_USE_CTL
{
#ifdef ART_USE_CTL
    ctl_pq_(65536);
    ctl_pq_inv_(65536);
    for (int i = 0; i < 65536; ++i) {
        float x = float(i) / 65535.f;
        ctl_pq_[i] = CTL_shaper_pq(x, false);
        ctl_pq_inv_[i] = CTL_shaper_pq(x, true);
    }    
#endif
}


#ifdef ART_USE_CTL

float CLUTStore::CTL_shaper(float a, bool inv)
{
    if (a >= 0.f && a <= 1.f) {
        return inv ? ctl_pq_inv_[a * 65535.f] : ctl_pq_[a * 65535.f];
    }
    return CTL_shaper_pq(a, inv);
}

#endif // ART_USE_CTL


//-----------------------------------------------------------------------------
// CLUTApplication
//-----------------------------------------------------------------------------

CLUTApplication::CLUTApplication(const Glib::ustring &clut_filename, const Glib::ustring &working_profile, float strength, int num_threads):
    clut_filename_(clut_filename),
    working_profile_(working_profile),
    ok_(false),
    clut_and_working_profiles_are_same_(false),
    multiThread_(num_threads > 1),
    strength_(strength)
{
    init(num_threads);
}


void CLUTApplication::init(int num_threads)
{
    hald_clut_ = CLUTStore::getInstance().getHaldClut(clut_filename_);
    if (!hald_clut_) {
#ifdef ART_USE_OCIO
        if (!OCIO_init())
#endif // ART_USE_OCIO
#ifdef ART_USE_CTL
        if (!CTL_init(num_threads))
#endif // ART_USE_CTL
            ok_ = false;
        return;
    }

    clut_and_working_profiles_are_same_ = hald_clut_->getProfile() == working_profile_;

    if (!clut_and_working_profiles_are_same_) {
        wprof_ = ICCStore::getInstance()->workingSpaceMatrix(working_profile_);
        wiprof_ = ICCStore::getInstance()->workingSpaceInverseMatrix(working_profile_);
        
        xyz2clut_ = ICCStore::getInstance()->workingSpaceInverseMatrix(hald_clut_->getProfile());
        clut2xyz_ = ICCStore::getInstance()->workingSpaceMatrix(hald_clut_->getProfile());

#ifdef __SSE2__
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                v_work2xyz_[i][j] = F2V(wprof_[i][j]);
                v_xyz2clut_[i][j] = F2V(xyz2clut_[i][j]);
                v_xyz2work_[i][j] = F2V(wiprof_[i][j]);
                v_clut2xyz_[i][j] = F2V(clut2xyz_[i][j]);
            }
        }
#endif
    }

    ok_ = true;
}


#ifdef ART_USE_OCIO

bool CLUTApplication::OCIO_init()
{
    auto proc = CLUTStore::getInstance().getOCIOLut(clut_filename_);
    if (!proc) {
        ok_ = false;
        return false;
    }

    try {
        ok_ = true;

        ocio_processor_ = proc->getOptimizedCPUProcessor(OCIO::BIT_DEPTH_F32, 
                                                         OCIO::BIT_DEPTH_F32,
                                                         OCIO::OPTIMIZATION_DEFAULT);
        init_matrices("ACESp0");
        return true;
    } catch (...) {
        ok_ = false;
        return false;
    }
}
#endif // ART_USE_OCIO


#ifdef ART_USE_CTL

bool CLUTApplication::CTL_init(int num_threads)
{
    try {
        Glib::ustring colorspace = "ACESp0";
        Glib::ustring lbl;
        auto func = CLUTStore::getInstance().getCTLLut(clut_filename_, num_threads, ctl_chunk_size_, ctl_params_, colorspace);
        if (func.empty()) {
            ok_ = false;
            return false;
        } else {
            ctl_func_ = std::move(func);
            init_matrices(colorspace);
            ok_ = true;
            return true;
        }
    } catch (...) {
        ok_ = false;
        return false;
    }
}


bool CLUTApplication::CTL_set_params(const std::vector<double> &values, Quality q)
{
    try {
        if (values.size() != ctl_params_.size() && !values.empty()) {
            if (settings->verbose) {
                std::cout << "Error in setting parameters for LUT " << clut_filename_ << ": " << (values.size() < ctl_params_.size() ? "not enough values" : "too many values") << std::endl;
            }
            return false;
        }
        for (size_t i = 0; i < ctl_params_.size(); ++i) {
            auto &desc = ctl_params_[i];
            auto v = i < values.size() ? values[i] : desc.value_default;
            switch (desc.type) {
            case CLUTParamType::PT_BOOL:
                for (auto f : ctl_func_) {
                    *reinterpret_cast<bool *>(f->inputArg(i+3)->data()) = v;
                }
                break;
            case CLUTParamType::PT_FLOAT:
                for (auto f : ctl_func_) {
                    *reinterpret_cast<float *>(f->inputArg(i+3)->data()) = v;
                }
                break;
            case CLUTParamType::PT_INT:
            case CLUTParamType::PT_CHOICE:
            default:
                for (auto f : ctl_func_) {
                    *reinterpret_cast<int *>(f->inputArg(i+3)->data()) = v;
                }
                break;
            }
        }
    } catch (std::exception &exc) {
        if (settings->verbose) {
            std::cout << "Error in setting parameters for LUT " << clut_filename_ << ": " << exc.what() << std::endl;
        }
        return false;
    } catch (...) {
        return false;
    }

    if (settings->ctl_scripts_fast_preview) {
        switch (q) {
        case Quality::LOW:
            CTL_init_lut(16);
            break;
        case Quality::MEDIUM:
            CTL_init_lut(50);
            break;
        default:
            break;
        }
    }

    return true;
}

namespace {

inline float CTL_shaper(float a, bool inv)
{
    return CLUTStore::getInstance().CTL_shaper(a, inv);
}


class CTLLutInitializer: public LUT3D::initializer {
public:
    CTLLutInitializer(const std::vector<float> *rgb):
        rgb_(rgb), i_(0) {}
    
    void operator()(float &r, float &g, float &b) override
    {
        r = rgb_[0][i_];
        g = rgb_[1][i_];
        b = rgb_[2][i_];
        ++i_;
    }

private:
    const std::vector<float> *rgb_;
    int i_;
};

} // namespace


void CLUTApplication::CTL_init_lut(int dim)
{
    std::vector<float> rgb[3];

    int sz = SQR(dim) * dim;
    for (int i = 0; i < 3; ++i) {
        rgb[i].reserve(sz);
    }
    
    for (int i = 0; i < dim; ++i) {
        float r = float(i)/(dim-1);
        for (int j = 0; j < dim; ++j) {
            float g = float(j)/(dim-1);
            for (int k = 0; k < dim; ++k) {
                float b = float(k)/(dim-1);
                rgb[0].push_back(CTL_shaper(r, true));
                rgb[1].push_back(CTL_shaper(g, true));
                rgb[2].push_back(CTL_shaper(b, true));
            }
        }
    }

    auto func = ctl_func_[0];

    for (int x = 0; x < sz; x += ctl_chunk_size_) {
        const auto n = (x + ctl_chunk_size_ < sz ? ctl_chunk_size_ : sz - x);
        for (int i = 0; i < 3; ++i) {
            memcpy(func->inputArg(i)->data(), &(rgb[i][x]), sizeof(float) * n);
        }
        func->callFunction(n);
        for (int i = 0; i < 3; ++i) {
            memcpy(&(rgb[i][x]), func->outputArg(i)->data(), sizeof(float) * n);
        }
    }

    CTLLutInitializer f(rgb);
    ctl_lut_.init(dim, f);
}

#endif // ART_USE_CTL


std::vector<CLUTParamDescriptor> CLUTApplication::get_param_descriptors() const
{
#ifdef ART_USE_CTL
    if (!ctl_func_.empty()) {
        return ctl_params_;
    }
#endif // ART_USE_CTL
    return {};
}


bool CLUTApplication::set_param_values(const std::vector<double> &values, Quality q)
{
#ifdef ART_USE_CTL
    if (!ctl_func_.empty()) {
        return CTL_set_params(values, q);
    }
#endif // ART_USE_CTL
    return values.empty();
}


std::vector<CLUTParamDescriptor> CLUTApplication::get_param_descriptors(const Glib::ustring &filename)
{
#ifdef ART_USE_CTL
    try {
        std::vector<CLUTParamDescriptor> params;
        int n;
        Glib::ustring colorspace;
        auto func = CLUTStore::getInstance().getCTLLut(filename, 1, n, params, colorspace);
        return params;
    } catch (...) {}
#endif // ART_USE_CTL
    return {};
}


#if defined ART_USE_OCIO || defined ART_USE_CTL

void CLUTApplication::init_matrices(const Glib::ustring &lut_profile)
{
    wprof_ = ICCStore::getInstance()->workingSpaceMatrix(working_profile_);
    wiprof_ = ICCStore::getInstance()->workingSpaceInverseMatrix(working_profile_);
    auto lprof = ICCStore::getInstance()->workingSpaceMatrix(lut_profile);
    auto liprof = ICCStore::getInstance()->workingSpaceInverseMatrix(lut_profile);
    
    auto ws = dot_product(liprof, wprof_);
    auto iws = dot_product(wiprof_, lprof);

    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            conv_[i][j] = ws[i][j];
            iconv_[i][j] = iws[i][j] * 65535.f;
        }
    }        
}

#endif // ART_USE_OCIO || ART_USE_CTL


void CLUTApplication::operator()(Imagefloat *img)
{
    if (!ok_) {
        return;
    }

    const int W = img->getWidth();
    const int H = img->getHeight();

#ifdef _OPENMP
#   pragma omp parallel for if (multiThread_)
#endif
    for (int y = 0; y < H; ++y) {
#ifdef _OPENMP
        int thread_id = omp_get_thread_num();
#else
        int thread_id = 0;
#endif
        apply(thread_id, W, img->r(y), img->g(y), img->b(y));
    }
}


inline void CLUTApplication::do_apply(int W, float *r, float *g, float *b)
{
    AlignedBuffer<float> buf_out_rgbx(4 * W); // Line buffer for CLUT
    AlignedBuffer<float> buf_clutr(W);
    AlignedBuffer<float> buf_clutg(W);
    AlignedBuffer<float> buf_clutb(W);
    float *out_rgbx = buf_out_rgbx.data;
    float *clutr = buf_clutr.data;
    float *clutg = buf_clutg.data;
    float *clutb = buf_clutb.data;
    
    if (!clut_and_working_profiles_are_same_) {
        // Convert from working to clut profile
        int j = 0;

#ifdef __SSE2__
        for (; j < W - 3; j += 4) {
            vfloat sourceR = LVF(r[j]);
            vfloat sourceG = LVF(g[j]);
            vfloat sourceB = LVF(b[j]);

            vfloat x;
            vfloat y;
            vfloat z;
            Color::rgbxyz(sourceR, sourceG, sourceB, x, y, z, v_work2xyz_);
            Color::xyz2rgb(x, y, z, sourceR, sourceG, sourceB, v_xyz2clut_);

            STVF(clutr[j], sourceR);
            STVF(clutg[j], sourceG);
            STVF(clutb[j], sourceB);
        }

#endif

        for (; j < W; j++) {
            float sourceR = r[j];
            float sourceG = g[j];
            float sourceB = b[j];

            float x, y, z;
            Color::rgbxyz(sourceR, sourceG, sourceB, x, y, z, wprof_);
            Color::xyz2rgb(x, y, z, clutr[j], clutg[j], clutb[j], xyz2clut_);
        }
    } else {
        memcpy(clutr, r, sizeof(float) * W);
        memcpy(clutg, g, sizeof(float) * W);
        memcpy(clutb, b, sizeof(float) * W);
    }

    for (int j = 0; j < W; j++) {
        float &sourceR = clutr[j];
        float &sourceG = clutg[j];
        float &sourceB = clutb[j];

        // Apply gamma sRGB (default RT)
        sourceR = Color::gamma_srgbclipped(sourceR);
        sourceG = Color::gamma_srgbclipped(sourceG);
        sourceB = Color::gamma_srgbclipped(sourceB);
    }

    hald_clut_->getRGB(strength_, W, clutr, clutg, clutb, out_rgbx);

    for (int j = 0; j < W; j++) {
        float &sourceR = clutr[j];
        float &sourceG = clutg[j];
        float &sourceB = clutb[j];

        // Apply inverse gamma sRGB
        sourceR = Color::igamma_srgb(out_rgbx[j * 4 + 0]);
        sourceG = Color::igamma_srgb(out_rgbx[j * 4 + 1]);
        sourceB = Color::igamma_srgb(out_rgbx[j * 4 + 2]);
    }

    if (!clut_and_working_profiles_are_same_) {
        // Convert from clut to working profile
        int j = 0;

#ifdef __SSE2__

        for (; j < W - 3; j += 4) {
            vfloat sourceR = LVF(clutr[j]);
            vfloat sourceG = LVF(clutg[j]);
            vfloat sourceB = LVF(clutb[j]);

            vfloat x;
            vfloat y;
            vfloat z;
            Color::rgbxyz(sourceR, sourceG, sourceB, x, y, z, v_clut2xyz_);
            Color::xyz2rgb(x, y, z, sourceR, sourceG, sourceB, v_xyz2work_);

            STVF(clutr[j], sourceR);
            STVF(clutg[j], sourceG);
            STVF(clutb[j], sourceB);
        }

#endif

        for (; j < W; j++) {
            float &sourceR = clutr[j];
            float &sourceG = clutg[j];
            float &sourceB = clutb[j];

            float x, y, z;
            Color::rgbxyz(sourceR, sourceG, sourceB, x, y, z, clut2xyz_);
            Color::xyz2rgb(x, y, z, sourceR, sourceG, sourceB, wiprof_);
        }
    }

    for (int j = 0; j < W; j++) {
        r[j] = clutr[j];
        g[j] = clutg[j];
        b[j] = clutb[j];
    }
}


#ifdef ART_USE_OCIO

inline void CLUTApplication::OCIO_apply(int W, float *r, float *g, float *b)
{
    const bool blend = strength_ < 1.f;

    Vec3<float> v;
    std::vector<float> data(W * 3);
    for (int x = 0, i = 0; x < W; ++x) {
        v[0] = r[x] / 65535.f;
        v[1] = g[x] / 65535.f;
        v[2] = b[x] / 65535.f;
        v = dot_product(conv_, v);
        data[i++] = v[0];
        data[i++] = v[1];
        data[i++] = v[2];
    }

    OCIO::PackedImageDesc pd(&data[0], W, 1, 3);
    ocio_processor_->apply(pd);
            
    for (int x = 0, i = 0; x < W; ++x) {
        v[0] = data[i++];
        v[1] = data[i++];
        v[2] = data[i++];
        v = dot_product(iconv_, v);
        // no need to renormalize to 65535 as this is already done in iconv_
        if (blend) {
            r[x] = intp(strength_, v[0], r[x]);
            g[x] = intp(strength_, v[1], g[x]);
            b[x] = intp(strength_, v[2], b[x]);
        } else {
            r[x] = v[0];
            g[x] = v[1];
            b[x] = v[2];
        }
    }
}

#endif // ART_USE_OCIO


#ifdef ART_USE_CTL

inline void CLUTApplication::CTL_apply(int thread_id, int W, float *r, float *g, float *b)
{
    if (!ctl_func_.empty()) {
        auto func = ctl_func_[thread_id];
        Vec3<float> v;
        std::vector<float> rgb[3];
        for (int i = 0; i < 3; ++i) {
            rgb[i].resize(W);
        }
        
        for (int x = 0; x < W; ++x) {
            v[0] = r[x] / 65535.f;
            v[1] = g[x] / 65535.f;
            v[2] = b[x] / 65535.f;
            v = dot_product(conv_, v);
            rgb[0][x] = v[0];
            rgb[1][x] = v[1];
            rgb[2][x] = v[2];
        }

        if (ctl_lut_) {
            for (int x = 0; x < W; ++x) {
                float r = CTL_shaper(rgb[0][x], false);
                float g = CTL_shaper(rgb[1][x], false);
                float b = CTL_shaper(rgb[2][x], false);
                ctl_lut_(r, g, b);
                rgb[0][x] = r;
                rgb[1][x] = g;
                rgb[2][x] = b;
            }            
        } else {
            for (int x = 0; x < W; x += ctl_chunk_size_) {
                const auto n = (x + ctl_chunk_size_ < W ? ctl_chunk_size_ : W - x);
                for (int i = 0; i < 3; ++i) {
                    memcpy(func->inputArg(i)->data(), &(rgb[i][x]), sizeof(float) * n);
                }
                func->callFunction(n);
                for (int i = 0; i < 3; ++i) {
                    memcpy(&(rgb[i][x]), func->outputArg(i)->data(), sizeof(float) * n);
                }
            }
        }
        
        const bool blend = strength_ < 1.f;
        
        for (int x = 0; x < W; ++x) {
            v[0] = rgb[0][x];
            v[1] = rgb[1][x];
            v[2] = rgb[2][x];
            v = dot_product(iconv_, v);
            // no need to renormalize to 65535 as this is already done in iconv_
            if (blend) {
                r[x] = intp(strength_, v[0], r[x]);
                g[x] = intp(strength_, v[1], g[x]);
                b[x] = intp(strength_, v[2], b[x]);
            } else {
                r[x] = v[0];
                g[x] = v[1];
                b[x] = v[2];
            }
        }

    }
}

#endif // ART_USE_CTL


void CLUTApplication::apply(int thread_id, int W, float *r, float *g, float *b)
{
    if (!ok_) {
        return;
    }

#ifdef ART_USE_OCIO
    if (ocio_processor_) {
        OCIO_apply(W, r, g, b);
        return;
    }
#endif
#ifdef ART_USE_CTL
    if (!ctl_func_.empty()) {
        CTL_apply(thread_id, W, r, g, b);
        return;
    }
#endif

    do_apply(W, r, g, b);
}

} // namespace rtengine
