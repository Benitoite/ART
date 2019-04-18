/** -*- C++ -*-
 *  
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2017 Alberto Griggio <alberto.griggio@gmail.com>
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

class LocalContrast: public ToolParamBlock, public AdjusterListener, public FoldableToolPanel, public CurveListener
{
private:
    MyComboBoxText *mode;
    Gtk::VBox *usm;
    Gtk::VBox *wavelets;
    Adjuster *radius;
    Adjuster *amount;
    Adjuster *darkness;
    Adjuster *lightness;
    Adjuster *contrast;
    FlatCurveEditor *curve;
    
    rtengine::ProcEvent EvLocalContrastEnabled;
    rtengine::ProcEvent EvLocalContrastMode;
    rtengine::ProcEvent EvLocalContrastRadius;
    rtengine::ProcEvent EvLocalContrastAmount;
    rtengine::ProcEvent EvLocalContrastDarkness;
    rtengine::ProcEvent EvLocalContrastLightness;
    rtengine::ProcEvent EvLocalContrastContrast;
    rtengine::ProcEvent EvLocalContrastCurve;
    
public:

    LocalContrast();

    void read(const rtengine::procparams::ProcParams *pp, const ParamsEdited *pedited=nullptr) override;
    void write(rtengine::procparams::ProcParams *pp, ParamsEdited *pedited=nullptr) override;
    void setDefaults(const rtengine::procparams::ProcParams *defParams, const ParamsEdited *pedited=nullptr) override;
    void setBatchMode(bool batchMode) override;

    void adjusterChanged(Adjuster *a, double newval) override;
    void adjusterAutoToggled(Adjuster* a, bool newval) override;
    void enabledChanged() override;
    void setAdjusterBehavior(bool radiusAdd, bool amountAdd, bool darknessAdd, bool lightnessAdd);
    void curveChanged() override;
    void autoOpenCurve() override;

private:
    void modeChanged();
};

