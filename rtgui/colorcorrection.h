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
#include "colorprovider.h"
#include "thresholdadjuster.h"

class ColorCorrection:
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel,
    public PParamsChangeListener,
    public ThresholdAdjusterListener,
    public ColorProvider
{
public:

    ColorCorrection();

    void read(const rtengine::procparams::ProcParams *pp) override;
    void write(rtengine::procparams::ProcParams *pp) override;
    void setDefaults(const rtengine::procparams::ProcParams *defParams) override;

    void adjusterChanged(Adjuster *a, double newval) override;
    void enabledChanged() override;
    void adjusterAutoToggled(Adjuster *a, bool newval) override {}

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
    void setAreaDrawListener(AreaDrawListener *l);
    void setDeltaEColorProvider(DeltaEColorProvider *p);

    void adjusterChanged(ThresholdAdjuster *a, double newBottom, double newTop) override;
    void adjusterChanged(ThresholdAdjuster *a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight) override {}
    void adjusterChanged(ThresholdAdjuster *a, int newBottom, int newTop) override {}
    void adjusterChanged(ThresholdAdjuster *a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight) override {}
    void adjusterChanged2(ThresholdAdjuster *a, int newBottomL, int newTopL, int newBottomR, int newTopR) override {}

    void colorForValue(double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller *caller) override;
    
private:
    void regionGet(int idx);
    void regionShow(int idx);
    void modeChanged();
    void syncSlidersToggled();
    
    rtengine::ProcEvent EvEnabled;
    rtengine::ProcEvent EvAB;
    rtengine::ProcEvent EvSaturation;
    rtengine::ProcEvent EvLightness;
    rtengine::ProcEvent EvSlope;
    rtengine::ProcEvent EvOffset;
    rtengine::ProcEvent EvPower;    
    rtengine::ProcEvent EvPivot;    
    rtengine::ProcEvent EvMode;

    rtengine::ProcEvent EvList;
    rtengine::ProcEvent EvHueMask;
    rtengine::ProcEvent EvChromaticityMask;
    rtengine::ProcEvent EvLightnessMask;
    rtengine::ProcEvent EvMaskBlur;
    rtengine::ProcEvent EvShowMask;
    rtengine::ProcEvent EvAreaMask;
    rtengine::ProcEvent EvDeltaEMask;

    std::vector<rtengine::procparams::ColorCorrectionParams::Region> data;

    friend class ColorCorrectionMasksContentProvider;
    std::unique_ptr<LabMasksContentProvider> labMasksContentProvider;
    LabMasksPanel *labMasks;
    
    Gtk::VBox *box;
    MyComboBoxText *mode;

    Gtk::VBox *box_combined;
    Gtk::VBox *box_rgb;
    Gtk::VBox *box_hsl;
    
    LabGrid *gridAB;
    Adjuster *saturation;
    Adjuster *slope;
    Adjuster *offset;
    Adjuster *power;
    Adjuster *pivot;

    Adjuster *slope_rgb[3];
    Adjuster *offset_rgb[3];
    Adjuster *power_rgb[3];
    Adjuster *pivot_rgb[3];

    Gtk::CheckButton *sync_rgb_sliders;
    
    Adjuster *lfactor[3];
    ThresholdAdjuster *huesat[3];
};

