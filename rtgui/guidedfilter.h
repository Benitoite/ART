/** -*- C++ -*-
 *  
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2018 Alberto Griggio <alberto.griggio@gmail.com>
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
#pragma once

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "mycurve.h"

class GuidedFilter: public ToolParamBlock, public AdjusterListener, public FoldableToolPanel, public CurveListener
{
private:
    Adjuster *smoothingRadius;
    Adjuster *smoothingEpsilon;
    Adjuster *smoothingIterations;
    Adjuster *smoothingLumaBlend;
    Adjuster *smoothingChromaBlend;

    Adjuster *decompRadius;
    Adjuster *decompEpsilon;
    Adjuster *decompDetailBoost;
    
    DiagonalCurveEditor *decompBaseCurve1;
    DiagonalCurveEditor *decompBaseCurve2;

    rtengine::ProcEvent EvEnabled;
    rtengine::ProcEvent EvSmoothingRadius;
    rtengine::ProcEvent EvSmoothingEpsilon;
    rtengine::ProcEvent EvSmoothingIterations;
    rtengine::ProcEvent EvSmoothingLumaBlend;
    rtengine::ProcEvent EvSmoothingChromaBlend;
    rtengine::ProcEvent EvDecompRadius;
    rtengine::ProcEvent EvDecompEpsilon;
    rtengine::ProcEvent EvDecompDetailBoost;
    rtengine::ProcEvent EvDecompBaseCurve1;
    rtengine::ProcEvent EvDecompBaseCurve2;
    
public:

    GuidedFilter();

    void read(const rtengine::procparams::ProcParams *pp, const ParamsEdited *pedited=nullptr);
    void write(rtengine::procparams::ProcParams *pp, ParamsEdited *pedited=nullptr);
    void setDefaults(const rtengine::procparams::ProcParams *defParams, const ParamsEdited *pedited=nullptr);
    void setBatchMode(bool batchMode);

    void adjusterChanged(Adjuster *a, double newval);
    void enabledChanged();
    void adjusterAutoToggled(Adjuster *a, bool newval) {}
    void curveChanged(CurveEditor *ce);
    void autoOpenCurve();
};

