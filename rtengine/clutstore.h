/** -*- C++ -*-
 *  
 *  This file is part of ART
 *
 *  ART is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ART is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ART.  If not, see <http://www.gnu.org/licenses/>.
 */

#pragma once

#include <memory>
#include <cstdint>

#include <gtkmm.h>

#include "cache.h"
#include "alignedbuffer.h"
#include "noncopyable.h"
#include "iccstore.h"
#include "imagefloat.h"
#include "clutparams.h"

#ifdef ART_USE_OCIO
#  include <OpenColorIO/OpenColorIO.h>
namespace OCIO = OCIO_NAMESPACE;
#endif // ART_USE_OCIO

#ifdef ART_USE_CTL
#  include <CtlSimdInterpreter.h>
#  include "LUT3D.h"
#endif


namespace rtengine {

class HaldCLUT final: public NonCopyable {
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

private:
    AlignedBuffer<std::uint16_t> clut_image;
    unsigned int clut_level;
    float flevel_minus_one;
    float flevel_minus_two;
    Glib::ustring clut_filename;
    Glib::ustring clut_profile;
};


class CLUTStore final: public NonCopyable {
public:
    static CLUTStore& getInstance();

    std::shared_ptr<HaldCLUT> getHaldClut(const Glib::ustring& filename) const;
#ifdef ART_USE_OCIO
    OCIO::ConstProcessorRcPtr getOCIOLut(const Glib::ustring &filename) const;
#endif // ART_USE_OCIO
#ifdef ART_USE_CTL
    std::vector<Ctl::FunctionCallPtr> getCTLLut(const Glib::ustring &filename, int num_threads, int &chunk_size, std::vector<CLUTParamDescriptor> &params, Glib::ustring &colorspace) const;
    float CTL_shaper(float a, bool inv);
#endif // ART_USE_CTL

    void clearCache();

    static Glib::ustring getClutDisplayName(const Glib::ustring &filename);
    
    static void splitClutFilename(
        const Glib::ustring& filename,
        Glib::ustring& name,
        Glib::ustring& extension,
        Glib::ustring& profile_name
    );

private:
    CLUTStore();

    mutable Cache<Glib::ustring, std::shared_ptr<HaldCLUT>> cache;
#ifdef ART_USE_OCIO
    typedef std::pair<OCIO::ConstProcessorRcPtr, std::string> OCIOCacheEntry;
    mutable Cache<Glib::ustring, OCIOCacheEntry> ocio_cache_;
#endif // ART_USE_OCIO
#ifdef ART_USE_CTL
    struct CTLCacheEntry {
        std::shared_ptr<Ctl::Interpreter> intp;
        std::string md5;
        std::vector<CLUTParamDescriptor> params;
        Glib::ustring colorspace;
    };
    mutable Cache<Glib::ustring, CTLCacheEntry> ctl_cache_;
    LUTf ctl_pq_;
    LUTf ctl_pq_inv_;
#endif // ART_USE_CTL
    mutable MyMutex mutex_;
};


class CLUTApplication {
public:
    enum class Quality {
        LOW,
        MEDIUM,
        HIGH
    };
    CLUTApplication(const Glib::ustring &clut_filename, const Glib::ustring &working_profile="", float strength=1.f, int num_threads=1);
    void operator()(Imagefloat *img);
    void apply(int thread_id, int W, float *r, float *g, float *b);
    operator bool() const { return ok_; }

    std::vector<CLUTParamDescriptor> get_param_descriptors() const;
    bool set_param_values(const std::vector<double> &values, Quality q=Quality::HIGH);

    static std::vector<CLUTParamDescriptor> get_param_descriptors(const Glib::ustring &filename);

private:
    void init(int num_threads);
    void do_apply(int W, float *r, float *g, float *b);
    Glib::ustring clut_filename_;
    Glib::ustring working_profile_;
    bool ok_;
    bool clut_and_working_profiles_are_same_;
    bool multiThread_;
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
#endif // __SSE2__

#ifdef ART_USE_OCIO
    bool OCIO_init();
    void OCIO_apply(int W, float *r, float *g, float *b);

    OCIO::ConstCPUProcessorRcPtr ocio_processor_;
#endif // ART_USE_OCIO

#ifdef ART_USE_CTL
    bool CTL_init(int num_threads);
    void CTL_apply(int thread_id, int W, float *r, float *g, float *b);
    bool CTL_set_params(const std::vector<double> &values, Quality q);
    void CTL_init_lut(int dim);
    std::vector<Ctl::FunctionCallPtr> ctl_func_;
    int ctl_chunk_size_;
    std::vector<CLUTParamDescriptor> ctl_params_;
    LUT3D ctl_lut_;
#endif // ART_USE_CTL

#if defined ART_USE_OCIO || defined ART_USE_CTL
    void init_matrices(const Glib::ustring &lut_profile);
    
    float conv_[3][3];
    float iconv_[3][3];
#endif
};

} // namespace rtengine
