/*
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
#include "brightcontrsat.h"
#include "adjuster.h"
#include <sigc++/slot.h>
#include <iomanip>
#include "ppversion.h"
#include "edit.h"
#include "eventmapper.h"

using namespace rtengine;
using namespace rtengine::procparams;

BrightnessContrastSaturation::BrightnessContrastSaturation():
    FoldableToolPanel(this, "brightcontrsat", M("TP_BRIGHTCONTRSAT_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    EvVibrance = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_BRIGHTCONTRSAT_VIBRANCE");
    EvToolEnabled.set_action(LUMINANCECURVE);
    autolevels = nullptr;
    
    brightness = Gtk::manage (new Adjuster (M("TP_BRIGHTCONTRSAT_BRIGHTNESS"), -100, 100, 1, 0));
//    pack_start (*brightness);
    contrast   = Gtk::manage (new Adjuster (M("TP_BRIGHTCONTRSAT_CONTRAST"), -100, 100, 1, 0));
//    pack_start (*contrast);
    saturation = Gtk::manage (new Adjuster (M("TP_BRIGHTCONTRSAT_SATURATION"), -100, 100, 1, 0));
    pack_start (*saturation);
    vibrance = Gtk::manage (new Adjuster (M("TP_BRIGHTCONTRSAT_VIBRANCE"), -100, 100, 1, 0));
    pack_start(*vibrance);

    brightness->setLogScale(2, 0, true);
    contrast->setLogScale(2, 0, true);
    saturation->setLogScale(2, 0, true);
    vibrance->setLogScale(2, 0, true);

    brightness->setAdjusterListener(this);
    contrast->setAdjusterListener(this);
    saturation->setAdjusterListener(this);
    vibrance->setAdjusterListener(this);
}


BrightnessContrastSaturation::~BrightnessContrastSaturation ()
{
}


void BrightnessContrastSaturation::read(const ProcParams* pp)
{
    disableListener ();

    setEnabled(pp->brightContrSat.enabled);
    brightness->setValue(pp->brightContrSat.brightness);
    contrast->setValue(pp->brightContrSat.contrast);
    saturation->setValue(pp->brightContrSat.saturation);
    vibrance->setValue(pp->brightContrSat.vibrance);

    enableListener ();
}


void BrightnessContrastSaturation::write(ProcParams* pp)
{
    pp->brightContrSat.enabled = getEnabled();
    pp->brightContrSat.brightness = (int)brightness->getValue ();
    pp->brightContrSat.contrast = (int)contrast->getValue ();
    pp->brightContrSat.saturation = (int)saturation->getValue ();
    pp->brightContrSat.vibrance = (int)vibrance->getValue ();
}


void BrightnessContrastSaturation::setDefaults (const ProcParams* defParams)
{
    brightness->setDefault (defParams->brightContrSat.brightness);
    contrast->setDefault (defParams->brightContrSat.contrast);
    saturation->setDefault(defParams->brightContrSat.saturation);
    vibrance->setDefault(defParams->brightContrSat.vibrance);
}


void BrightnessContrastSaturation::adjusterChanged(Adjuster* a, double newval)
{
    // Switch off auto exposure if user changes sliders manually
    if (autolevels && autolevels->get_active() && (a == brightness || a == contrast)) {
        autolevels->set_active (false);
    }

    if (!listener || !getEnabled()) {
        return;
    }

    Glib::ustring costr = Glib::ustring::format ((int)a->getValue());

    if (a == brightness) {
        listener->panelChanged(EvBrightness, costr);
    } else if (a == contrast) {
        listener->panelChanged(EvContrast, costr);
    } else if (a == saturation) {
        listener->panelChanged(EvSaturation, costr);
    } else if (a == vibrance) {
        listener->panelChanged(EvVibrance, costr);
    }
}

void BrightnessContrastSaturation::adjusterAutoToggled(Adjuster* a, bool newval)
{
}


void BrightnessContrastSaturation::trimValues (rtengine::procparams::ProcParams* pp)
{
    brightness->trimValue(pp->brightContrSat.brightness);
    contrast->trimValue(pp->brightContrSat.contrast);
    saturation->trimValue(pp->brightContrSat.saturation);
    vibrance->trimValue(pp->brightContrSat.vibrance);
}


void BrightnessContrastSaturation::autoExpChanged(double expcomp, int bright, int contr, int black, int hlcompr, int hlcomprthresh, bool hlrecons)
{
    disableListener();
    setEnabled(true);
    brightness->setValue(bright);
    contrast->setValue(contr);
    enableListener();
}
