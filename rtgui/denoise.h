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
#pragma once

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "colorprovider.h"
#include "guiutils.h"
#include "options.h"

class Denoise:
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel,
    public rtengine::AutoChromaListener,
    public CurveListener
{
public:
    Denoise();
    ~Denoise() override;

    void read(const rtengine::procparams::ProcParams* pp) override;
    void write(rtengine::procparams::ProcParams* pp) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams) override;
    void adjusterChanged(Adjuster* a, double newval) override;
    void adjusterAutoToggled(Adjuster* a, bool newval) override;
    void enabledChanged() override;

    void chromaChanged(double autchroma, double autred, double autblue) override;
    bool chromaComputed(double chroma, double red, double blue);
    
    void noiseChanged(double nresid, double highresid) override {}
    void noiseTilePrev(int tileX, int tileY, int prevX, int prevY, int sizeT, int sizeP) override {}

    void trimValues(rtengine::procparams::ProcParams* pp) override;

    void curveChanged() override;
    void autoOpenCurve() override;

private:
    void colorSpaceChanged();
    void aggressiveChanged();
    void chrominanceMethodChanged();
    void medianTypeChanged();
    void medianMethodChanged();
    void smoothingEnabledToggled();
    void smoothingMethodChanged();
    
    rtengine::ProcEvent EvSmoothingMethod;
    rtengine::ProcEvent EvGuidedLumaRadius;
    rtengine::ProcEvent EvGuidedLumaStrength;
    rtengine::ProcEvent EvGuidedChromaRadius;
    rtengine::ProcEvent EvGuidedChromaStrength;
    rtengine::ProcEvent EvChrominanceAutoFactor;
    rtengine::ProcEvent EvLuminanceDetailThreshold;
    rtengine::ProcEvent EvLuminanceCurve;

    MyComboBoxText *colorSpace;
    MyComboBoxText *aggressive;
    Adjuster *gamma;
    Adjuster *luminance;
    Adjuster *luminanceDetail;
    Adjuster *luminanceDetailThreshold;
    MyComboBoxText *chrominanceMethod;
    Adjuster *chrominanceAutoFactor;
    Adjuster *chrominance;
    Adjuster *chrominanceRedGreen;
    Adjuster *chrominanceBlueYellow;
    MyExpander *smoothingEnabled;
    MyComboBoxText *smoothingMethod;
    Gtk::VBox *medianBox;
    MyComboBoxText *medianType;
    MyComboBoxText *medianMethod;
    Adjuster *medianIterations;
    Gtk::VBox *guidedBox;
    Adjuster *guidedLumaRadius;
    Adjuster *guidedLumaStrength;
    Adjuster *guidedChromaRadius;
    Adjuster *guidedChromaStrength;
    CurveEditorGroup *cg;
    FlatCurveEditor *luminanceCurve;

    IdleRegister idle_register;
};

