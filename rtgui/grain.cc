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
#include "grain.h"
#include "eventmapper.h"
#include <iomanip>
#include <cmath>

using namespace rtengine;
using namespace rtengine::procparams;

FilmGrain::FilmGrain(): FoldableToolPanel(this, "grain", M("TP_GRAIN_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    EvEnabled = m->newEvent(SHARPENING, "HISTORY_MSG_GRAIN_ENABLED");
    EvStrength = m->newEvent(SHARPENING, "HISTORY_MSG_GRAIN_STRENGTH");
    EvISO = m->newEvent(SHARPENING, "HISTORY_MSG_GRAIN_ISO");
    EvScale = m->newEvent(SHARPENING, "HISTORY_MSG_GRAIN_SCALE");
    
    iso = Gtk::manage(new Adjuster(M("TP_GRAIN_ISO"), 20., 6400., 1., 400.));
    iso->setAdjusterListener(this);
    iso->show();

    strength = Gtk::manage(new Adjuster(M("TP_GRAIN_STRENGTH"), 0., 100., 1., 25.));
    strength->setAdjusterListener(this);
    strength->show();

    scale = Gtk::manage(new Adjuster(M("TP_GRAIN_SCALE"), 0., 100., 1., 100.));
    scale->setAdjusterListener(this);
    scale->show();
    
    pack_start(*iso);
    pack_start(*strength);
    pack_start(*scale);
}


void FilmGrain::read(const ProcParams *pp, const ParamsEdited *pedited)
{
    disableListener();

    if (pedited) {
        iso->setEditedState(pedited->grain.iso ? Edited : UnEdited);
        strength->setEditedState(pedited->grain.strength ? Edited : UnEdited);
        scale->setEditedState(pedited->grain.scale ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->grain.enabled);
    }

    setEnabled(pp->grain.enabled);
    iso->setValue(pp->grain.iso);
    strength->setValue(pp->grain.strength);
    scale->setValue(pp->grain.scale);

    enableListener();
}


void FilmGrain::write(ProcParams *pp, ParamsEdited *pedited)
{
    pp->grain.enabled = getEnabled();
    pp->grain.iso = iso->getValue();
    pp->grain.strength = strength->getValue();
    pp->grain.scale = scale->getValue();

    if (pedited) {
        pedited->grain.enabled = !get_inconsistent();
        pedited->grain.iso = iso->getEditedState();
        pedited->grain.strength = strength->getEditedState();
        pedited->grain.scale = scale->getEditedState();
    }
}

void FilmGrain::setDefaults(const ProcParams *defParams, const ParamsEdited *pedited)
{
    iso->setDefault(defParams->grain.iso);
    strength->setDefault(defParams->grain.strength);
    scale->setDefault(defParams->grain.scale);

    if (pedited) {
        iso->setDefaultEditedState(pedited->grain.iso ? Edited : UnEdited);
        strength->setDefaultEditedState(pedited->grain.strength ? Edited : UnEdited);
        scale->setDefaultEditedState(pedited->grain.scale ? Edited : UnEdited);
    } else {
        iso->setDefaultEditedState(Irrelevant);
        strength->setDefaultEditedState(Irrelevant);
        scale->setDefaultEditedState(Irrelevant);
    }
}


void FilmGrain::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        if (a == strength) {
            listener->panelChanged(EvStrength, a->getTextValue());
        } else if (a == iso) {
            listener->panelChanged(EvISO, a->getTextValue());
        } else if (a == scale) {
            listener->panelChanged(EvScale, a->getTextValue());
        }
    }
}


void FilmGrain::enabledChanged ()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(EvEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(EvEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvEnabled, M("GENERAL_DISABLED"));
        }
    }
}


void FilmGrain::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);

    iso->showEditedCB();
    strength->showEditedCB();
    scale->showEditedCB();
}
