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
#include "dehaze.h"
#include "eventmapper.h"
#include <iomanip>
#include <cmath>

using namespace rtengine;
using namespace rtengine::procparams;

Dehaze::Dehaze(): FoldableToolPanel(this, "dehaze", M("TP_DEHAZE_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    EvDehazeEnabled = m->newEvent(HDR, "HISTORY_MSG_DEHAZE_ENABLED");
    EvDehazeStrength = m->newEvent(HDR, "HISTORY_MSG_DEHAZE_STRENGTH");
    EvDehazeShowDepthMap = m->newEvent(HDR, "HISTORY_MSG_DEHAZE_SHOW_DEPTH_MAP");
    EvDehazeDepth = m->newEvent(HDR, "HISTORY_MSG_DEHAZE_DEPTH");
    EvDehazeLuminance = m->newEvent(HDR, "HISTORY_MSG_DEHAZE_LUMINANCE");
    
    strength = Gtk::manage(new Adjuster(M("TP_DEHAZE_STRENGTH"), 0., 100., 1., 50.));
    strength->setAdjusterListener(this);
    strength->show();

    depth = Gtk::manage(new Adjuster(M("TP_DEHAZE_DEPTH"), 0., 100., 1., 25.));
    depth->setAdjusterListener(this);
    depth->show();

    Gtk::HBox *hb = Gtk::manage (new Gtk::HBox ());
    hb->pack_start(*Gtk::manage(new Gtk::Label(M("TP_DEHAZE_MODE") + ": ")), Gtk::PACK_SHRINK);
    luminance = Gtk::manage(new MyComboBoxText());
    luminance->append(M("TP_DEHAZE_RGB"));
    luminance->append(M("TP_DEHAZE_LUMINANCE"));
    hb->pack_start(*luminance);
    pack_start(*hb);
    luminance->signal_changed().connect(sigc::mem_fun(*this, &Dehaze::luminanceChanged));
    hb->show();
    luminance->show();

    showDepthMap = Gtk::manage(new Gtk::CheckButton(M("TP_DEHAZE_SHOW_DEPTH_MAP")));
    showDepthMap->signal_toggled().connect(sigc::mem_fun(*this, &Dehaze::showDepthMapChanged));
    showDepthMap->show();

    pack_start(*strength);
    pack_start(*depth);
    //pack_start(*luminance);
    pack_start(*showDepthMap);
}


void Dehaze::read(const ProcParams *pp)
{
    disableListener();

    setEnabled(pp->dehaze.enabled);
    strength->setValue(pp->dehaze.strength);
    depth->setValue(pp->dehaze.depth);
    showDepthMap->set_active(pp->dehaze.showDepthMap);
    luminance->set_active(pp->dehaze.luminance ? 1 : 0);

    enableListener();
}


void Dehaze::write(ProcParams *pp)
{
    pp->dehaze.strength = strength->getValue();
    pp->dehaze.depth = depth->getValue();
    pp->dehaze.enabled = getEnabled();
    pp->dehaze.showDepthMap = showDepthMap->get_active();
    pp->dehaze.luminance = luminance->get_active_row_number() == 1;
}

void Dehaze::setDefaults(const ProcParams *defParams)
{
    strength->setDefault(defParams->dehaze.strength);
    depth->setDefault(defParams->dehaze.depth);
}


void Dehaze::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        if (a == strength) {
            listener->panelChanged(EvDehazeStrength, a->getTextValue());
        } else if (a == depth) {
            listener->panelChanged(EvDehazeDepth, a->getTextValue());
        }
    }
}


void Dehaze::enabledChanged ()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(EvDehazeEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(EvDehazeEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvDehazeEnabled, M("GENERAL_DISABLED"));
        }
    }
}


void Dehaze::showDepthMapChanged()
{
    if (listener) {
        listener->panelChanged(EvDehazeShowDepthMap, showDepthMap->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}


void Dehaze::luminanceChanged()
{
    if (listener) {
        listener->panelChanged(EvDehazeLuminance, luminance->get_active_row_number() == 1 ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}
