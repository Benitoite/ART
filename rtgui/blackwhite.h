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
#ifndef _BLACKWHITE_H_
#define _BLACKWHITE_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "guiutils.h"
#include "curveeditor.h"
#include "curveeditorgroup.h"
#include "mycurve.h"
#include "colorprovider.h"

class BlackWhite final :
    public ToolParamBlock,
    public AdjusterListener,
    public FoldableToolPanel,
    public ColorProvider
{
public:

    BlackWhite ();
    ~BlackWhite () override;

    void read(const rtengine::procparams::ProcParams* pp) override;
    void write(rtengine::procparams::ProcParams* pp) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams) override;
    void neutral_pressed ();

    void updateRGBLabel      ();
    void adjusterChanged     (Adjuster* a, double newval) override;
    void adjusterAutoToggled (Adjuster* a, bool newval) override;
    void trimValues          (rtengine::procparams::ProcParams* pp) override;
    void enabledChanged      () override;
    void methodChanged       ();
    void filterChanged       ();
    void settingChanged      ();

    Glib::ustring getSettingString ();
    Glib::ustring getFilterString  ();

private:
    void showFilter();
    void hideFilter();
    void showMixer(int nChannels, bool RGBIsSensitive = true);
    void hideMixer();
    void showGamma();
    void hideGamma();

    Gtk::Button*         neutral;
    Gtk::Label*          RGBLabels;

    Adjuster *mixerRed;
    Adjuster *mixerGreen;
    Adjuster *mixerBlue;
    Adjuster *gammaRed;
    Adjuster *gammaGreen;
    Adjuster *gammaBlue;
    MyComboBoxText*   method;
    sigc::connection  methodconn;
    Gtk::HBox*        filterHBox;
    Gtk::HSeparator*  filterSep, *filterSep2;
    MyComboBoxText*   filter;
    sigc::connection  filterconn;
    Gtk::HBox*        settingHBox;
    MyComboBoxText*   setting;
    sigc::connection  settingconn;
    Gtk::Frame* mixerFrame;
    Gtk::VBox * mixerVBox;
    Gtk::Frame* gammaFrame;

    Gtk::Image *imgIcon[11];

    Gtk::HSeparator* enabledccSep;
    Gtk::CheckButton* enabledcc;
    bool lastEnabledcc, lastAuto;
    sigc::connection enaccconn, tcmodeconn, tcmodeconn2, autoconn, neutralconn;

    double nextredbw;
    double nextgreenbw;
    double nextbluebw;
    int nextcount = 0;

    IdleRegister idle_register;
};

#endif
