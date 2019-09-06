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
#include "labgrid.h"

class ColorCorrection: public ToolParamBlock, public AdjusterListener, public FoldableToolPanel, public PParamsChangeListener
{
public:

    ColorCorrection();

    void read(const rtengine::procparams::ProcParams *pp);
    void write(rtengine::procparams::ProcParams *pp);
    void setDefaults(const rtengine::procparams::ProcParams *defParams);

    void adjusterChanged(Adjuster *a, double newval);
    void enabledChanged();
    void adjusterAutoToggled(Adjuster *a, bool newval) {}

    void setEditProvider(EditDataProvider *provider) override;
    void setListener(ToolPanelListener *tpl) override;

    PParamsChangeListener *getPParamsChangeListener() override { return this; }
    void procParamsChanged(
        const rtengine::procparams::ProcParams* params,
        const rtengine::ProcEvent& ev,
        const Glib::ustring& descr,
        const ParamsEdited* paramsEdited = nullptr) override;
    void clearParamChanges() override {}

    void updateGeometry(int fullWidth, int fullHeight);

private:
    void regionGet(int idx);
    void regionShow(int idx);
    void channelChanged();
    
    rtengine::ProcEvent EvEnabled;
    rtengine::ProcEvent EvAB;
    rtengine::ProcEvent EvSaturation;
    rtengine::ProcEvent EvLightness;
    rtengine::ProcEvent EvSlope;
    rtengine::ProcEvent EvOffset;
    rtengine::ProcEvent EvPower;    
    rtengine::ProcEvent EvPivot;    
    rtengine::ProcEvent EvChannel;

    rtengine::ProcEvent EvList;
    rtengine::ProcEvent EvHueMask;
    rtengine::ProcEvent EvChromaticityMask;
    rtengine::ProcEvent EvLightnessMask;
    rtengine::ProcEvent EvMaskBlur;
    rtengine::ProcEvent EvShowMask;
    rtengine::ProcEvent EvAreaMask;

    std::vector<rtengine::procparams::ColorCorrectionParams::Region> data;
    int showMaskIdx;

    friend class ColorCorrectionMasksContentProvider;
    std::unique_ptr<LabMasksContentProvider> labMasksContentProvider;
    LabMasksPanel *labMasks;
    
    Gtk::VBox *box;
    LabGrid *gridAB;
    Adjuster *saturation;
    Adjuster *slope;
    Adjuster *offset;
    Adjuster *power;
    Adjuster *pivot;
    MyComboBoxText *channel;
};

