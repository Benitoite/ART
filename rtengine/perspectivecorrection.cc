/* -*- C++ -*-
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

// taken from darktable (src/iop/ashift.c)
/*
  This file is part of darktable,
  copyright (c) 2016 Ulrich Pegelow.

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
// Inspiration to this module comes from the program ShiftN (http://www.shiftn.de) by
// Marcus Hebel.

// Thanks to Marcus for his support when implementing part of the ShiftN functionality
// to darktable.


#include "perspectivecorrection.h"
#include "improcfun.h"
#include "rt_math.h"
#include <string.h>
#include <math.h>
#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "../rtgui/threadutils.h"
#include "settings.h"

namespace rtengine { extern const Settings *settings; }

#define _(msg) (msg)
#define dt_control_log(msg) \
    if (settings->verbose) { \
        printf("%s\n", msg);       \
        fflush(stdout);            \
    }


namespace rtengine {

namespace {

#define ROTATION_RANGE 10                   // allowed min/max default range for rotation parameter
#define ROTATION_RANGE_SOFT 20              // allowed min/max range for rotation parameter with manual adjustment
#define LENSSHIFT_RANGE 0.5                 // allowed min/max default range for lensshift parameters
#define LENSSHIFT_RANGE_SOFT 1              // allowed min/max range for lensshift parameters with manual adjustment
#define SHEAR_RANGE 0.2                     // allowed min/max range for shear parameter
#define SHEAR_RANGE_SOFT 0.5                // allowed min/max range for shear parameter with manual adjustment
#define MIN_LINE_LENGTH 5                   // the minimum length of a line in pixels to be regarded as relevant
#define MAX_TANGENTIAL_DEVIATION 30         // by how many degrees a line may deviate from the +/-180 and +/-90 to be regarded as relevant
#define LSD_SCALE 0.99                      // LSD: scaling factor for line detection
#define LSD_SIGMA_SCALE 0.6                 // LSD: sigma for Gaussian filter is computed as sigma = sigma_scale/scale
#define LSD_QUANT 2.0                       // LSD: bound to the quantization error on the gradient norm
#define LSD_ANG_TH 22.5                     // LSD: gradient angle tolerance in degrees
#define LSD_LOG_EPS 0.0                     // LSD: detection threshold: -log10(NFA) > log_eps
#define LSD_DENSITY_TH 0.7                  // LSD: minimal density of region points in rectangle
#define LSD_N_BINS 1024                     // LSD: number of bins in pseudo-ordering of gradient modulus
#define LSD_GAMMA 0.45                      // gamma correction to apply on raw images prior to line detection
#define RANSAC_RUNS 400                     // how many iterations to run in ransac
#define RANSAC_EPSILON 2                    // starting value for ransac epsilon (in -log10 units)
#define RANSAC_EPSILON_STEP 1               // step size of epsilon optimization (log10 units)
#define RANSAC_ELIMINATION_RATIO 60         // percentage of lines we try to eliminate as outliers
#define RANSAC_OPTIMIZATION_STEPS 5         // home many steps to optimize epsilon
#define RANSAC_OPTIMIZATION_DRY_RUNS 50     // how man runs per optimization steps
#define RANSAC_HURDLE 5                     // hurdle rate: the number of lines below which we do a complete permutation instead of random sampling
#define MINIMUM_FITLINES 4                  // minimum number of lines needed for automatic parameter fit
#define NMS_EPSILON 1e-3                    // break criterion for Nelder-Mead simplex
#define NMS_SCALE 1.0                       // scaling factor for Nelder-Mead simplex
#define NMS_ITERATIONS 400                  // number of iterations for Nelder-Mead simplex
#define NMS_CROP_EPSILON 100.0              // break criterion for Nelder-Mead simplex on crop fitting
#define NMS_CROP_SCALE 0.5                  // scaling factor for Nelder-Mead simplex on crop fitting
#define NMS_CROP_ITERATIONS 100             // number of iterations for Nelder-Mead simplex on crop fitting
#define NMS_ALPHA 1.0                       // reflection coefficient for Nelder-Mead simplex
#define NMS_BETA 0.5                        // contraction coefficient for Nelder-Mead simplex
#define NMS_GAMMA 2.0                       // expansion coefficient for Nelder-Mead simplex
#define DEFAULT_F_LENGTH 28.0               // focal length we assume if no exif data are available


inline int mat3inv(float *const dst, const float *const src)
{
    std::array<std::array<float, 3>, 3> tmpsrc;
    std::array<std::array<float, 3>, 3> tmpdst;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            tmpsrc[i][j] = src[3 * i + j];
        }
    }
    if (invertMatrix(tmpsrc, tmpdst)) {
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                dst[3 * i + j] = tmpdst[i][j];
            }
        }
        return 0;
    } else {
        return 1;
    }
}


// the darktable ashift iop (adapted to RT), which does most of the work
#include "ashift_dt.c"

} // namespace


PerspectiveCorrection::PerspectiveCorrection():
    ok_(false),
    scale_(1.0),
    offx_(0.0),
    offy_(0.0)
{
}


void PerspectiveCorrection::init(int width, int height, const procparams::PerspectiveParams &params, bool fill)
{
    homography((float *)ihomograph_, params.angle, params.vertical / 100.0, -params.horizontal / 100.0, params.shear / 100.0, DEFAULT_F_LENGTH, 0.f, 1.f, width, height, ASHIFT_HOMOGRAPH_INVERTED);

    ok_ = true;
    calc_scale(width, height, params, fill);
}


inline void PerspectiveCorrection::correct(double &x, double &y, double scale, double offx, double offy)
{
    if (ok_) {       
        float pin[3], pout[3];
        pout[0] = x;
        pout[1] = y;
        pout[0] *= scale;//_;
        pout[1] *= scale;//_;
        pout[0] += offx;//_;
        pout[1] += offy;//_;
        pout[2] = 1.f;
        mat3mulv(pin, (float *)ihomograph_, pout);
        pin[0] /= pin[2];
        pin[1] /= pin[2];
        x = pin[0];
        y = pin[1];
    }
}


void PerspectiveCorrection::operator()(double &x, double &y)
{
    correct(x, y, scale_, offx_, offy_);
}


namespace {

std::vector<Coord2D> get_corners(int w, int h)
{
    int x1 = 0, y1 = 0;
    int x2 = x1 + w - 1;
    int y2 = y1 + h - 1;

    std::vector<Coord2D> corners = {
        Coord2D(x1, y1),
        Coord2D(x1, y2),
        Coord2D(x2, y2),
        Coord2D(x2, y1)
    };
    return corners;
}

void init_dt_structures(dt_iop_ashift_params_t *p, dt_iop_ashift_gui_data_t *g,
                        const procparams::PerspectiveParams *params)
{
    dt_iop_ashift_params_t dp = {
        0.0f,
        0.0f,
        0.0f,
        0.0f,
        DEFAULT_F_LENGTH,
        1.0f,
        100.0f,
        1.0f,
        ASHIFT_MODE_GENERIC,
        0,
        ASHIFT_CROP_OFF,
        0.0f,
        1.0f,
        0.0f,
        1.0f
    };
    *p = dp;

    g->buf = NULL;
    g->buf_width = 0;
    g->buf_height = 0;
    g->buf_x_off = 0;
    g->buf_y_off = 0;
    g->buf_scale = 1.0f;
    g->buf_hash = 0;
    g->isflipped = 0;
    g->lastfit = ASHIFT_FIT_NONE;
    g->fitting = 0;
    g->lines = NULL;
    g->lines_count =0;
    g->horizontal_count = 0;
    g->vertical_count = 0;
    g->grid_hash = 0;
    g->lines_hash = 0;
    g->rotation_range = ROTATION_RANGE_SOFT;
    g->lensshift_v_range = LENSSHIFT_RANGE_SOFT;
    g->lensshift_h_range = LENSSHIFT_RANGE_SOFT;
    g->shear_range = SHEAR_RANGE_SOFT;
    g->lines_suppressed = 0;
    g->lines_version = 0;
    g->show_guides = 0;
    g->isselecting = 0;
    g->isdeselecting = 0;
    g->isbounding = ASHIFT_BOUNDING_OFF;
    g->near_delta = 0;
    g->selecting_lines_version = 0;
    g->points = NULL;
    g->points_idx = NULL;
    g->points_lines_count = 0;
    g->points_version = 0;
    g->jobcode = ASHIFT_JOBCODE_NONE;
    g->jobparams = 0;
    g->adjust_crop = FALSE;
    g->lastx = g->lasty = -1.0f;
    g->crop_cx = g->crop_cy = 1.0f;

    if (params) {
        p->rotation = params->angle;
        p->lensshift_v = params->vertical / 100.0;
        p->lensshift_h = -params->horizontal / 100.0;
        p->shear = params->shear / 100.0;
    }
}

} // namespace

void PerspectiveCorrection::calc_scale(int w, int h, const procparams::PerspectiveParams &params, bool fill)
{
    double cw, ch;
    get_view_size(w, h, params, cw, ch);

    if (!fill) {
        scale_ = max(cw / double(w), ch / double(h));
        offx_ = (cw - w * scale_) * 0.5;
        offy_ = (ch - h * scale_) * 0.5;
    } else {
        dt_iop_ashift_params_t p;
        dt_iop_ashift_gui_data_t g;
        init_dt_structures(&p, &g, &params);
        dt_iop_module_t module = { &g, false };
        g.buf_width = w;
        g.buf_height = h;
        p.cropmode = ASHIFT_CROP_ASPECT;
        do_crop(&module, &p);
        offx_ = p.cl * cw;
        offy_ = p.ct * ch;
        scale_ = (p.cr - p.cl) * cw/double(w);
        scale_ *= 0.99; // add a safety factor to avoid border issues -- TODO
    }
}


void PerspectiveCorrection::get_view_size(int w, int h, const procparams::PerspectiveParams &params, double &cw, double &ch)
{
    double min_x = RT_INFINITY, max_x = -RT_INFINITY;
    double min_y = RT_INFINITY, max_y = -RT_INFINITY;

    auto corners = get_corners(w, h);

    float homo[3][3];
    homography((float *)homo, params.angle, params.vertical / 100.0, -params.horizontal / 100.0, params.shear / 100.0, DEFAULT_F_LENGTH, 0.f, 1.f, w, h, ASHIFT_HOMOGRAPH_FORWARD);
    
    for (auto &c : corners) {
        float pin[3] = { float(c.x), float(c.y), 1.f };
        float pout[3];
        mat3mulv(pout, (float *)homo, pin);
        double x = pout[0] / pout[2];
        double y = pout[1] / pout[2];
        min_x = min(min_x, x);
        max_x = max(max_x, x);
        min_y = min(min_y, y);
        max_y = max(max_y, y);
    }

    cw = floor(max_x - min_x + 1);
    ch = floor(max_y - min_y + 1);
}    


procparams::PerspectiveParams PerspectiveCorrection::autocompute(ImageSource *src, Direction dir, const procparams::ProcParams *pparams)
{
    dt_iop_ashift_params_t p;
    dt_iop_ashift_gui_data_t g;
    init_dt_structures(&p, &g, nullptr);
    dt_iop_module_t module;
    module.gui_data = &g;
    module.is_raw = src->isRAW();

    int tr = getCoarseBitMask(pparams->coarse);
    int fw, fh;
    src->getFullSize(fw, fh, tr);
    int skip = max(float(max(fw, fh)) / 900.f + 0.5f, 1.f);
    PreviewProps pp(0, 0, fw, fh, skip);
    int w, h;
    src->getSize(pp, w, h);
    std::unique_ptr<Imagefloat> img(new Imagefloat(w, h));

    ProcParams neutral;
    neutral.raw.bayersensor.method = RAWParams::BayerSensor::getMethodString(RAWParams::BayerSensor::Method::FAST);
    neutral.raw.xtranssensor.method = RAWParams::XTransSensor::getMethodString(RAWParams::XTransSensor::Method::FAST);
    neutral.icm.outputProfile = ColorManagementParams::NoICMString;    
    src->getImage(src->getWB(), tr, img.get(), pp, neutral.toneCurve, neutral.raw);
    src->convertColorSpace(img.get(), pparams->icm, src->getWB());

    neutral.rotate = pparams->rotate;
    neutral.distortion = pparams->distortion;
    neutral.lensProf = pparams->lensProf;
    ImProcFunctions ipf(&neutral, true);
    if (ipf.needsTransform()) {
        Imagefloat *tmp = new Imagefloat(w, h);
        ipf.transform(img.get(), tmp, 0, 0, 0, 0, w, h, w, h,
                      src->getMetaData(), src->getRotateDegree(), false);
        img.reset(tmp);
    }

    // allocate the gui buffer
    g.buf = static_cast<float *>(malloc(sizeof(float) * w * h * 4));
    g.buf_width = w;
    g.buf_height = h;

    img->normalizeFloatTo1();
    
#ifdef _OPENMP
#   pragma omp parallel for
#endif
    for (int y = 0; y < h; ++y) {
        for (int x = 0; x < w; ++x) {
            int i = (y * w + x) * 4;
            g.buf[i] = img->r(y, x);
            g.buf[i+1] = img->g(y, x);
            g.buf[i+2] = img->b(y, x);
            g.buf[i+3] = 1.f;
        }
    }

    dt_iop_ashift_fitaxis_t fitaxis = ASHIFT_FIT_NONE;
    switch (dir) {
    case HORIZONTAL:
        fitaxis = ASHIFT_FIT_HORIZONTALLY;
        break;  
    case VERTICAL:
        fitaxis = ASHIFT_FIT_VERTICALLY;
        break;
    default:
        fitaxis = ASHIFT_FIT_BOTH_SHEAR;
        break;
    }
    auto res = do_get_structure(&module, &p, ASHIFT_ENHANCE_EDGES) && do_fit(&module, &p, fitaxis);
    procparams::PerspectiveParams retval;

    // cleanup the gui
    if (g.lines) free(g.lines);
    if (g.points) free(g.points);
    if (g.points_idx) free(g.points_idx);
    free(g.buf);

    if (res) {
        retval.horizontal = -p.lensshift_h * 100;
        retval.vertical = p.lensshift_v * 100;
        retval.angle = p.rotation;
        retval.shear = p.shear * 100;
    }
    return retval;
}


void PerspectiveCorrection::autocrop(int width, int height, const procparams::PerspectiveParams &params, int &x, int &y, int &w, int &h)
{
    double cw, ch;
    get_view_size(width, height, params, cw, ch);
    double s = min(double(width)/cw, double(height)/ch);
    dt_iop_ashift_params_t p;
    dt_iop_ashift_gui_data_t g;
    init_dt_structures(&p, &g, &params);
    dt_iop_module_t module = { &g, false };
    g.buf_width = width;
    g.buf_height = height;
    p.cropmode = ASHIFT_CROP_ASPECT;
    do_crop(&module, &p);
    s *= 0.99; // make the crop a bit smaller just to stay safe
    cw *= s;
    ch *= s;
    double ox = p.cl * cw;
    double oy = p.ct * ch;
    x = ceil(ox - (cw - width)/s * 0.5);
    y = ceil(oy - (ch - height)/s * 0.5);
    w = floor((p.cr - p.cl) * cw);
    h = floor((p.cb - p.ct) * ch);
}

} // namespace rtengine
