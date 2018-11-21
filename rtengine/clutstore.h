#pragma once

#include <memory>
#include <cstdint>

#include <gtkmm.h>

#include "cache.h"
#include "alignedbuffer.h"
#include "noncopyable.h"
#include "iccstore.h"

namespace rtengine
{

class HaldCLUT final :
    public NonCopyable
{
public:
    HaldCLUT();
    ~HaldCLUT();

    bool load(const Glib::ustring& filename);

    explicit operator bool() const;

    Glib::ustring getFilename() const;
    Glib::ustring getProfile() const;

    void getRGB(
        float strength,
        std::size_t line_size,
        const float* r,
        const float* g,
        const float* b,
        float* out_rgbx
    ) const;

    static void splitClutFilename(
        const Glib::ustring& filename,
        Glib::ustring& name,
        Glib::ustring& extension,
        Glib::ustring& profile_name
    );

private:
    AlignedBuffer<std::uint16_t> clut_image;
    unsigned int clut_level;
    float flevel_minus_one;
    float flevel_minus_two;
    Glib::ustring clut_filename;
    Glib::ustring clut_profile;
};


class CLUTStore final :
    public NonCopyable
{
public:
    static CLUTStore& getInstance();

    std::shared_ptr<HaldCLUT> getClut(const Glib::ustring& filename) const;

    void clearCache();

private:
    CLUTStore();

    mutable Cache<Glib::ustring, std::shared_ptr<HaldCLUT>> cache;
};


class HaldCLUTApplication {
public:
    HaldCLUTApplication(const Glib::ustring &clut_filename, const Glib::ustring &working_profile);
    void init(float strength, int tile_size);
    void operator()(float *r, float *g, float *b, int istart, int jstart, int tW, int tH);
    operator bool() const { return ok_; }

private:
    Glib::ustring clut_filename_;
    Glib::ustring working_profile_;
    bool ok_;
    bool clut_and_working_profiles_are_same_;
    int TS_;
    float strength_;
    std::shared_ptr<HaldCLUT> hald_clut_;
    TMatrix wprof_;
    TMatrix wiprof_;
    TMatrix xyz2clut_;
    TMatrix clut2xyz_;
#ifdef __SSE2__
    vfloat v_work2xyz_[3][3] ALIGNED16;
    vfloat v_xyz2clut_[3][3] ALIGNED16;
    vfloat v_clut2xyz_[3][3] ALIGNED16;
    vfloat v_xyz2work_[3][3] ALIGNED16;
#endif
};


}
