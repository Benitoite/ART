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
#ifndef _PARAMEDITED_H_
#define _PARAMEDITED_H_

#include <glibmm.h>
#include <vector>
//#include "../rtengine/procparams.h"
//#include "../rtengine/rtengine.h"

class ParamsEdited {

public:
    bool general;
    bool exposure;
    bool brightContrSat;
    bool toneCurve;
    bool labCurve;
    bool localContrast;
    bool rgbCurves;
    bool sharpening;
    bool prsharpening;
    bool wb;
    bool defringe;
    bool denoise;
    bool textureBoost;
    bool fattal;
    bool logenc;
    bool impulseDenoise;
    bool sh;
    bool toneEqualizer;
    bool crop;
    bool coarse;
    bool commonTrans;
    bool rotate;
    bool distortion;
    bool lensProf;
    bool perspective;
    bool gradient;
    bool pcvignette;
    bool cacorrection;
    bool vignetting;
    bool chmixer;
    bool blackwhite;
    bool resize;
    bool icm;

    //bool raw;
    bool demosaic;
    bool rawBlack;
    bool rawWhite;
    bool rawPreprocessing;
    bool hotDeadPixelFilter;
    bool darkframe;
    bool flatfield;
    bool rawCA;
    
    bool dirpyrequalizer;
    bool filmSimulation;
    bool softlight;
    bool dehaze;
    bool grain;
    bool smoothing;
    bool colorcorrection;
    bool metadata;
    bool exif;
    bool iptc;

    explicit ParamsEdited(bool value = false);

    void set(bool v);
    // void initFrom(const std::vector<rtengine::procparams::ProcParams>& src);
    // void combine(rtengine::procparams::ProcParams& toEdit, const rtengine::procparams::ProcParams& mods, bool forceSet);
};

#endif
