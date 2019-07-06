/** -*- C++ -*-
 *  
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 Alberto Griggio <alberto.griggio@gmail.com>
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
#include "fastexport.h"
#include "options.h"

namespace {

void adjust_fast_params(rtengine::procparams::ProcParams &params)
{
    if (!options.fastexport_use_fast_pipeline) {
        if (options.fastexport_bypass_sharpening) {
            params.sharpening.enabled = false;
        }

        if (options.fastexport_bypass_sharpenMicro) {
            params.sharpenMicro.enabled = false;
        }

        //if (options.fastexport_bypass_lumaDenoise) params.lumaDenoise.enabled = false;
        //if (options.fastexport_bypass_colorDenoise) params.colorDenoise.enabled = false;
        if (options.fastexport_bypass_defringe) {
            params.defringe.enabled = false;
        }

        if (options.fastexport_bypass_dirpyrDenoise) {
            params.denoise.enabled = false;
        }

        if (options.fastexport_bypass_dirpyrequalizer) {
            params.dirpyrequalizer.enabled = false;
        }

        //if (options.fastexport_bypass_raw_bayer_all_enhance) params.raw.bayersensor.all_enhance = false;
        if (options.fastexport_bypass_raw_bayer_dcb_iterations) {
            params.raw.bayersensor.dcb_iterations = 0;
        }

        if (options.fastexport_bypass_raw_bayer_dcb_enhance) {
            params.raw.bayersensor.dcb_enhance = false;
        }

        if (options.fastexport_bypass_raw_bayer_lmmse_iterations) {
            params.raw.bayersensor.lmmse_iterations = 0;
        }

        if (options.fastexport_bypass_raw_bayer_linenoise) {
            params.raw.bayersensor.linenoise = 0;
        }

        if (options.fastexport_bypass_raw_bayer_greenthresh) {
            params.raw.bayersensor.greenthresh = 0;
        }

        if (options.fastexport_bypass_raw_ccSteps) {
            params.raw.bayersensor.ccSteps = params.raw.xtranssensor.ccSteps = 0;
        }

        if (options.fastexport_bypass_raw_ca) {
            params.raw.ca_autocorrect = false;
            params.raw.cared = 0;
            params.raw.cablue = 0;
        }

        if (options.fastexport_bypass_raw_df) {
            params.raw.df_autoselect = false;
            params.raw.dark_frame = "";
        }

        if (options.fastexport_bypass_raw_ff) {
            params.raw.ff_AutoSelect = false;
            params.raw.ff_file = "";
        }

        params.raw.bayersensor.method = options.fastexport_raw_bayer_method;
        params.raw.xtranssensor.method = options.fastexport_raw_xtrans_method;
        params.icm.inputProfile = options.fastexport_icm_input_profile;
        params.icm.workingProfile = options.fastexport_icm_working_profile;
        params.icm.outputProfile = options.fastexport_icm_output_profile;
        params.icm.outputIntent = options.fastexport_icm_outputIntent;
        params.icm.outputBPC = options.fastexport_icm_outputBPC;
    }

    if (params.resize.enabled) {
        params.resize.width = rtengine::min(params.resize.width, options.fastexport_resize_width);
        params.resize.height = rtengine::min(params.resize.height, options.fastexport_resize_height);
    } else {
        params.resize.width = options.fastexport_resize_width;
        params.resize.height = options.fastexport_resize_height;
    }

    params.resize.enabled = options.fastexport_resize_enabled;
    params.resize.scale = options.fastexport_resize_scale;
    params.resize.appliesTo = options.fastexport_resize_appliesTo;
    params.resize.method = options.fastexport_resize_method;
    params.resize.dataspec = options.fastexport_resize_dataspec;
    params.resize.allowUpscaling = false;
}

} // namespace

rtengine::ProcessingJob *create_processing_job(const Glib::ustring &fname, bool is_raw, rtengine::procparams::ProcParams params, bool fast)
{
    if (fast) {
        adjust_fast_params(params);
    }

    auto ret = rtengine::ProcessingJob::create(fname, is_raw, params, fast && options.fastexport_use_fast_pipeline);
    return ret;
}


rtengine::ProcessingJob *create_processing_job(rtengine::InitialImage *initialImage, rtengine::procparams::ProcParams params, bool fast)
{
    if (fast) {
        adjust_fast_params(params);
    }

    auto ret = rtengine::ProcessingJob::create(initialImage, params, fast && options.fastexport_use_fast_pipeline);
    return ret;
}
