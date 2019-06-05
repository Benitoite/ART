/* -*- C++ -*-
 *  
 *  This file is part of RawTherapee.
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
 *
 *  (C) 2010 Emil Martinec <ejmartin@uchicago.edu>
 */

#ifndef DIRPYREQUALIZER_H_INCLUDED
#define DIRPYREQUALIZER_H_INCLUDED

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "thresholdadjuster.h"
#include "colorprovider.h"
#include "labmaskspanel.h"

class DirPyrEqualizer : public ToolParamBlock, public ThresholdAdjusterListener, public AdjusterListener, public FoldableToolPanel, public PParamsChangeListener

{

protected:
    void levelsGet(int idx);
    void levelsShow(int idx);
    
    rtengine::ProcEvent EvList;
    rtengine::ProcEvent EvHueMask;
    rtengine::ProcEvent EvChromaticityMask;
    rtengine::ProcEvent EvLightnessMask;
    rtengine::ProcEvent EvMaskBlur;
    rtengine::ProcEvent EvShowMask;
    rtengine::ProcEvent EvAreaMask;

    std::vector<rtengine::procparams::DirPyrEqualizerParams::Levels> levelsData;
    int showMaskIdx;

    friend class DirPyrEqMasksContentProvider;
    std::unique_ptr<LabMasksContentProvider> labMasksContentProvider;
    LabMasksPanel *labMasks;
    Gtk::VBox *box;
    Adjuster* multiplier[6];
    Adjuster* threshold;

public:

    DirPyrEqualizer();
    ~DirPyrEqualizer() override;

    void read(const rtengine::procparams::ProcParams* pp) override;
    void write(rtengine::procparams::ProcParams* pp) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams) override;
    void trimValues          (rtengine::procparams::ProcParams* pp) override;
    void adjusterChanged (Adjuster* a, double newval) override;
    void adjusterAutoToggled(Adjuster* a, bool newval) override;
    void enabledChanged() override;
    void lumaneutralPressed ();
    void lumacontrastPlusPressed ();
    void lumacontrastMinusPressed ();

    void setEditProvider(EditDataProvider *provider) override;

    PParamsChangeListener *getPParamsChangeListener() override { return this; }
    void procParamsChanged(
        const rtengine::procparams::ProcParams* params,
        const rtengine::ProcEvent& ev,
        const Glib::ustring& descr,
        const ParamsEdited* paramsEdited = nullptr) override;
    void clearParamChanges() override {}

    void updateGeometry(int fullWidth, int fullHeight);

    void adjusterChanged(ThresholdAdjuster* a, double newBottom, double newTop) override {}
    void adjusterChanged(ThresholdAdjuster* a, double newBottomLeft, double newTopLeft, double newBottomRight, double newTopRight) override {}
    void adjusterChanged(ThresholdAdjuster* a, int newBottom, int newTop) override {}
    void adjusterChanged(ThresholdAdjuster* a, int newBottomLeft, int newTopLeft, int newBottomRight, int newTopRight) override {}
    void adjusterChanged2(ThresholdAdjuster* a, int newBottomL, int newTopL, int newBottomR, int newTopR) override {}
};

#endif
