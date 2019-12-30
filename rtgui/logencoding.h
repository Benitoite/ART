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

class LogEncoding: public ToolParamBlock, public AdjusterListener, public rtengine::AutoLogListener, public FoldableToolPanel
{
protected:
    Gtk::ToggleButton *autocompute;
    Adjuster *sourceGray;
    Adjuster *targetGray;
    Adjuster *blackEv;
    Adjuster *whiteEv;
    Adjuster *localContrast;

    rtengine::ProcEvent EvEnabled;
    rtengine::ProcEvent EvAuto;
    rtengine::ProcEvent EvAutoGrayOn;
    rtengine::ProcEvent EvAutoGrayOff;
    rtengine::ProcEvent EvAutoBatch;
    rtengine::ProcEvent EvSourceGray;
    rtengine::ProcEvent EvSourceGrayAuto;
    rtengine::ProcEvent EvTargetGray;
    rtengine::ProcEvent EvBlackEv;
    rtengine::ProcEvent EvWhiteEv;
    rtengine::ProcEvent EvLocalContrast;

    sigc::connection autoconn;
    
public:
    LogEncoding();

    void read(const rtengine::procparams::ProcParams *pp);
    void write(rtengine::procparams::ProcParams *pp);
    void setDefaults(const rtengine::procparams::ProcParams *defParams);

    void adjusterChanged(Adjuster* a, double newval);
    void adjusterAutoToggled(Adjuster* a, bool newval);
    void enabledChanged();

    void logEncodingChanged(const rtengine::LogEncodingParams &params);
    void autocomputeToggled();
};

