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
#include "shadowshighlights.h"
#include "eventmapper.h"

using namespace rtengine;
using namespace rtengine::procparams;

ShadowsHighlights::ShadowsHighlights () : FoldableToolPanel(this, "shadowshighlights", M("TP_SHADOWSHLIGHTS_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    EvLevels = m->newEvent(RGBCURVE, "HISTORY_MSG_SH_LEVELS");

    for (size_t i = 0; i < levels.size(); ++i) {
        levels[i] = Gtk::manage(new Adjuster(M("TP_SHADOWSHLIGHTS_LEVEL_" + std::to_string(i)), -100, 100, 1, 0));
        levels[i]->setAdjusterListener(this);
        pack_start(*levels[i]);
    }
    
    show_all_children ();
}

void ShadowsHighlights::read (const ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener ();

    if (pedited) {
        set_inconsistent(multiImage && !pedited->sh.enabled);
        for (size_t i = 0; i < levels.size(); ++i) {
            levels[i]->setEditedState(pedited->sh.levels ? Edited : UnEdited);
        }
    }

    setEnabled (pp->sh.enabled);

    for (size_t i = 0; i < levels.size(); ++i) {
        levels[i]->setValue(pp->sh.levels[i]);
    }
    
    enableListener ();
}

void ShadowsHighlights::write (ProcParams* pp, ParamsEdited* pedited)
{
    for (size_t i = 0; i < levels.size(); ++i) {
        pp->sh.levels[i] = levels[i]->getValue();
    }
    pp->sh.enabled = getEnabled();

    if (pedited) {
        pedited->sh.enabled = !get_inconsistent();
        pedited->sh.levels = false;
        for (size_t i = 0; i < levels.size(); ++i) {
            pedited->sh.levels = pedited->sh.levels || levels[i]->getEditedState();
        }
    }
}

void ShadowsHighlights::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{
    for (size_t i = 0; i < levels.size(); ++i) {
        levels[i]->setDefault(defParams->sh.levels[i]);
    }

    if (pedited) {
        for (size_t i = 0; i < levels.size(); ++i) {
            levels[i]->setDefaultEditedState(pedited->sh.levels ? Edited : UnEdited);
        }
    } else {
        for (size_t i = 0; i < levels.size(); ++i) {
            levels[i]->setDefaultEditedState(Irrelevant);
        }
    }
}

void ShadowsHighlights::adjusterChanged (Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        Glib::ustring s;
        for (size_t i = 0; i < levels.size(); ++i) {
            s += Glib::ustring::format((int)a->getValue()) + " ";
        }
        listener->panelChanged(EvLevels, s);
    }
}

void ShadowsHighlights::adjusterAutoToggled(Adjuster* a, bool newval)
{
}

void ShadowsHighlights::enabledChanged ()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvSHEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvSHEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvSHEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void ShadowsHighlights::setBatchMode (bool batchMode)
{

    ToolPanel::setBatchMode (batchMode);
    for (size_t i = 0; i < levels.size(); ++i) {
        levels[i]->showEditedCB();
    }
}

void ShadowsHighlights::setAdjusterBehavior (bool hadd, bool sadd)
{

    // highlights->setAddMode(hadd);
    // shadows->setAddMode(sadd);
}

void ShadowsHighlights::trimValues (rtengine::procparams::ProcParams* pp)
{
    for (size_t i = 0; i < levels.size(); ++i) {
        levels[i]->trimValue(pp->sh.levels[i]);
    }
}
