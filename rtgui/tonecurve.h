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
#include "mycurve.h"
#include "guiutils.h"

class ToneCurve: public ToolParamBlock, public FoldableToolPanel, public CurveListener, public AdjusterListener
{
private:
    IdleRegister idle_register;

protected:
    Adjuster *contrast;
    MyComboBoxText* toneCurveMode;
    MyComboBoxText* toneCurveMode2;
    Gtk::ToggleButton *histmatching;
    bool fromHistMatching;

    sigc::connection tcmodeconn, tcmode2conn;
    sigc::connection histmatchconn;
    CurveEditorGroup* curveEditorG;
    CurveEditorGroup* curveEditorG2;
    DiagonalCurveEditor* shape;
    DiagonalCurveEditor* shape2;
    CurveEditorGroup *satcurveG;
    FlatCurveEditor *satcurve;

    rtengine::ProcEvent EvHistMatching;
    rtengine::ProcEvent EvHistMatchingBatch;
    rtengine::ProcEvent EvSatCurve;

    // used temporarily in eventing
    rtengine::procparams::ToneCurveParams::TcMode nextToneCurveMode;
    std::vector<double> nextToneCurve;

    void setHistmatching(bool enabled);

public:
    ToneCurve();
    ~ToneCurve() override;

    void read(const rtengine::procparams::ProcParams* pp) override;
    void write(rtengine::procparams::ProcParams* pp) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams) override;
    void trimValues(rtengine::procparams::ProcParams* pp) override;
    void autoOpenCurve() override;
    void setEditProvider(EditDataProvider *provider) override;

    float blendPipetteValues(CurveEditor *ce, float chan1, float chan2, float chan3) override;

    void enableAll(bool yes=true);
    void curveChanged (CurveEditor* ce) override;
    void curveMode1Changed ();
    bool curveMode1Changed_ ();
    void curveMode2Changed ();
    bool curveMode2Changed_ ();
    void expandCurve (bool isExpanded);
    bool isCurveExpanded ();
    void updateCurveBackgroundHistogram(
        const LUTu& histToneCurve,
        const LUTu& histLCurve,
        const LUTu& histCCurve,
        const LUTu& histLCAM,
        const LUTu& histCCAM,
        const LUTu& histRed,
        const LUTu& histGreen,
        const LUTu& histBlue,
        const LUTu& histLuma,
        const LUTu& histLRETI
    );

    void histmatchingToggled();
    void autoMatchedToneCurveChanged(rtengine::procparams::ToneCurveParams::TcMode curveMode, const std::vector<double>& curve);
    void setRaw (bool raw);

    void adjusterChanged(Adjuster *a, double newval) override;
};
