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
#include "labmaskspanel.h"

class Smoothing: public ToolParamBlock, public AdjusterListener, public FoldableToolPanel, public PParamsChangeListener
{
public:

    Smoothing();

    void read(const rtengine::procparams::ProcParams *pp) override;
    void write(rtengine::procparams::ProcParams *pp) override;
    void setDefaults(const rtengine::procparams::ProcParams *defParams) override;

    void adjusterChanged(Adjuster *a, double newval) override;
    void enabledChanged() override;
    void adjusterAutoToggled(Adjuster *a, bool newval) override {}

    void setEditProvider(EditDataProvider *provider) override;

    PParamsChangeListener *getPParamsChangeListener() override { return this; }
    void procParamsChanged(
        const rtengine::procparams::ProcParams* params,
        const rtengine::ProcEvent& ev,
        const Glib::ustring& descr,
        const ParamsEdited* paramsEdited = nullptr) override;
    void clearParamChanges() override {}

    void updateGeometry(int fullWidth, int fullHeight);
    void setAreaDrawListener(AreaDrawListener *l);
    void setDeltaEColorProvider(DeltaEColorProvider *p);

private:
    void regionGet(int idx);
    void regionShow(int idx);
    void channelChanged();
    
    rtengine::ProcEvent EvEnabled;
    rtengine::ProcEvent EvChannel;
    rtengine::ProcEvent EvRadius;
    rtengine::ProcEvent EvEpsilon;

    rtengine::ProcEvent EvList;
    rtengine::ProcEvent EvHueMask;
    rtengine::ProcEvent EvChromaticityMask;
    rtengine::ProcEvent EvLightnessMask;
    rtengine::ProcEvent EvMaskBlur;
    rtengine::ProcEvent EvShowMask;
    rtengine::ProcEvent EvAreaMask;
    rtengine::ProcEvent EvDeltaEMask;

    std::vector<rtengine::procparams::GuidedSmoothingParams::Region> data;

    friend class SmoothingMasksContentProvider;
    std::unique_ptr<LabMasksContentProvider> labMasksContentProvider;
    LabMasksPanel *labMasks;
    
    MyComboBoxText *channel;
    Adjuster *radius;
    Adjuster *epsilon;
    Gtk::VBox *box;
};

