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

class Exposure: public ToolParamBlock, public AdjusterListener, public FoldableToolPanel {
private:
    IdleRegister idle_register;

protected:
    Gtk::CheckButton *hrenabled;
    MyComboBoxText *method;
    sigc::connection methconn;
    sigc::connection enaconn;
    bool lasthrEnabled;

    Gtk::HBox *abox;
    Gtk::HBox *hlrbox;

    Gtk::ToggleButton *autolevels;
    Gtk::Label *lclip;
    MySpinButton *sclip;
    Gtk::Button *neutral;
    Adjuster *expcomp;
    Adjuster *black;
    Adjuster *hlcompr;
    Adjuster *hlcomprthresh;
    Adjuster *shcompr;
    Gtk::CheckButton *clampOOG;

    bool clipDirty, lastAuto;
    sigc::connection autoconn, neutralconn;

    rtengine::ProcEvent EvClampOOG;

    // used temporarily in eventing
    double nextExpcomp;
    int nextHlcompr;
    int nextHlcomprthresh;
    int nextBlack;
    bool nextHLRecons;

public:
    Exposure();
    ~Exposure() override;

    void read(const rtengine::procparams::ProcParams* pp) override;
    void write(rtengine::procparams::ProcParams* pp) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams) override;
    void trimValues(rtengine::procparams::ProcParams* pp) override;

    void adjusterChanged(Adjuster* a, double newval) override;
    void adjusterAutoToggled(Adjuster* a, bool newval) override;
    void neutral_pressed();
    void autolevels_toggled();
    void clip_changed();
    bool clip_changed_();
    void waitForAutoExp();
    void enableAll();

    void autoExpChanged(double expcomp, int bright, int contr, int black, int hlcompr, int hlcomprthresh, bool hlrecons);

    void setRaw(bool raw);

    void hrenabledChanged();
    void methodChanged();
    void clampOOGChanged();

    Gtk::ToggleButton *getAutoLevelsButton() { return autolevels; }
};

