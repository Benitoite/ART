#include <algorithm>

#include "clutstore.h"

#include "iccstore.h"
#include "imagefloat.h"
#include "opthelper.h"
#include "rt_math.h"
#include "stdimagesource.h"

#include "../rtgui/options.h"

namespace
{

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

        img_src.getImage(curr_wb, TR_NONE, img_float.get(), pp, rtengine::procparams::ToneCurveParams(), rtengine::procparams::RAWParams());

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

}

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
        splitClutFilename(filename, name, ext, clut_profile);

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

void rtengine::HaldCLUT::splitClutFilename(
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

    if (!name.empty()) {
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
    }
}

rtengine::CLUTStore& rtengine::CLUTStore::getInstance()
{
    static CLUTStore instance;
    return instance;
}

std::shared_ptr<rtengine::HaldCLUT> rtengine::CLUTStore::getClut(const Glib::ustring& filename) const
{
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

void rtengine::CLUTStore::clearCache()
{
    cache.clear();
}

rtengine::CLUTStore::CLUTStore() :
    cache(options.clutCacheSize)
{
}


//-----------------------------------------------------------------------------
// HaldCLUTApplication
//-----------------------------------------------------------------------------

namespace rtengine {

HaldCLUTApplication::HaldCLUTApplication(const Glib::ustring &clut_filename, const Glib::ustring &working_profile):
    clut_filename_(clut_filename),
    working_profile_(working_profile),
    ok_(false),
    clut_and_working_profiles_are_same_(false),
    TS_(0),
    strength_(0)
{
}


void HaldCLUTApplication::init(float strength, int tile_size)
{
    hald_clut_ = CLUTStore::getInstance().getClut(clut_filename_);
    if (!hald_clut_) {
        ok_ = false;
        return;
    }

    strength_ = strength;
    TS_ = tile_size;

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


void HaldCLUTApplication::operator()(float *r, float *g, float *b, int istart, int jstart, int tW, int tH)
{
    if (!ok_) {
        return;
    }

    float out_rgbx[4 * TS_] ALIGNED16; // Line buffer for CLUT
    float clutr[TS_] ALIGNED16;
    float clutg[TS_] ALIGNED16;
    float clutb[TS_] ALIGNED16;
    
    for (int i = istart, ti = 0; i < tH; i++, ti++) {
        std::cout << "HERE: " << i << std::endl;
        
        if (!clut_and_working_profiles_are_same_) {
            // Convert from working to clut profile
            int j = jstart;
            int tj = 0;

#ifdef __SSE2__
            for (; j < tW - 3; j += 4, tj += 4) {
                vfloat sourceR = LVF(r[ti * TS_ + tj]);
                vfloat sourceG = LVF(g[ti * TS_ + tj]);
                vfloat sourceB = LVF(b[ti * TS_ + tj]);

                vfloat x;
                vfloat y;
                vfloat z;
                Color::rgbxyz(sourceR, sourceG, sourceB, x, y, z, v_work2xyz_);
                Color::xyz2rgb(x, y, z, sourceR, sourceG, sourceB, v_xyz2clut_);

                STVF(clutr[tj], sourceR);
                STVF(clutg[tj], sourceG);
                STVF(clutb[tj], sourceB);
            }

#endif

            for (; j < tW; j++, tj++) {
                float sourceR = r[ti * TS_ + tj];
                float sourceG = g[ti * TS_ + tj];
                float sourceB = b[ti * TS_ + tj];

                float x, y, z;
                Color::rgbxyz(sourceR, sourceG, sourceB, x, y, z, wprof_);
                Color::xyz2rgb(x, y, z, clutr[tj], clutg[tj], clutb[tj], xyz2clut_);
            }
        } else {
            memcpy(clutr, &r[ti * TS_], sizeof(float) * TS_);
            memcpy(clutg, &g[ti * TS_], sizeof(float) * TS_);
            memcpy(clutb, &b[ti * TS_], sizeof(float) * TS_);
        }

        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
            float &sourceR = clutr[tj];
            float &sourceG = clutg[tj];
            float &sourceB = clutb[tj];

            // Apply gamma sRGB (default RT)
            sourceR = Color::gamma_srgbclipped(sourceR);
            sourceG = Color::gamma_srgbclipped(sourceG);
            sourceB = Color::gamma_srgbclipped(sourceB);
        }

        hald_clut_->getRGB(strength_, std::min(TS_, tW - jstart), clutr, clutg, clutb, out_rgbx);

        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
            float &sourceR = clutr[tj];
            float &sourceG = clutg[tj];
            float &sourceB = clutb[tj];

            // Apply inverse gamma sRGB
            sourceR = Color::igamma_srgb(out_rgbx[tj * 4 + 0]);
            sourceG = Color::igamma_srgb(out_rgbx[tj * 4 + 1]);
            sourceB = Color::igamma_srgb(out_rgbx[tj * 4 + 2]);
        }

        if (!clut_and_working_profiles_are_same_) {
            // Convert from clut to working profile
            int j = jstart;
            int tj = 0;

#ifdef __SSE2__

            for (; j < tW - 3; j += 4, tj += 4) {
                vfloat sourceR = LVF(clutr[tj]);
                vfloat sourceG = LVF(clutg[tj]);
                vfloat sourceB = LVF(clutb[tj]);

                vfloat x;
                vfloat y;
                vfloat z;
                Color::rgbxyz(sourceR, sourceG, sourceB, x, y, z, v_clut2xyz_);
                Color::xyz2rgb(x, y, z, sourceR, sourceG, sourceB, v_xyz2work_);

                STVF(clutr[tj], sourceR);
                STVF(clutg[tj], sourceG);
                STVF(clutb[tj], sourceB);
            }

#endif

            for (; j < tW; j++, tj++) {
                float &sourceR = clutr[tj];
                float &sourceG = clutg[tj];
                float &sourceB = clutb[tj];

                float x, y, z;
                Color::rgbxyz(sourceR, sourceG, sourceB, x, y, z, clut2xyz_);
                Color::xyz2rgb(x, y, z, sourceR, sourceG, sourceB, wiprof_);
            }
        }

        for (int j = jstart, tj = 0; j < tW; j++, tj++) {
            setUnlessOOG(r[ti * TS_ + tj], g[ti * TS_ + tj], b[ti * TS_ + tj], clutr[tj], clutg[tj], clutb[tj]);
        }
    }
}

} // namespace rtengine
