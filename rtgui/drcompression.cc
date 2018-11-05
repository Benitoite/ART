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
#include "drcompression.h"
#include "eventmapper.h"
#include <iomanip>
#include <cmath>

using namespace rtengine;
using namespace rtengine::procparams;

DRCompression::DRCompression(): FoldableToolPanel(this, "fattal", M("TP_TM_FATTAL_LABEL"), true, true)
{
    auto m = ProcEventMapper::getInstance();
    EvTMFattalAnchor = m->newEvent(HDR, "HISTORY_MSG_TM_FATTAL_ANCHOR");

    amount = Gtk::manage(new Adjuster (M("TP_TM_FATTAL_AMOUNT"), 1., 100., 1., 20.));
    threshold = Gtk::manage(new Adjuster (M("TP_TM_FATTAL_THRESHOLD"), -100., 300., 1., 30.0));
    Gtk::Image *al = Gtk::manage(new RTImage("circle-black-small.png"));
    Gtk::Image *ar = Gtk::manage(new RTImage("circle-white-small.png"));
    anchor = Gtk::manage(new Adjuster(M("TP_TM_FATTAL_ANCHOR"), 1, 100, 1, 50, al, ar));

    amount->setAdjusterListener(this);
    threshold->setAdjusterListener(this);
    anchor->setAdjusterListener(this);

    threshold->setLogScale(10, 0);

    amount->show();
    threshold->show();
    anchor->show();

    pack_start(*amount);
    pack_start(*threshold);
    pack_start(*anchor);
}

void DRCompression::read(const ProcParams *pp, const ParamsEdited *pedited)
{
    disableListener();

    if (pedited) {
        threshold->setEditedState(pedited->drcomp.threshold ? Edited : UnEdited);
        amount->setEditedState(pedited->drcomp.amount ? Edited : UnEdited);
        anchor->setEditedState(pedited->drcomp.anchor ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->drcomp.enabled);
    }

    setEnabled(pp->drcomp.enabled);
    threshold->setValue(pp->drcomp.threshold);
    amount->setValue(pp->drcomp.amount);
    anchor->setValue(pp->drcomp.anchor);

    enableListener();
}

void DRCompression::write(ProcParams *pp, ParamsEdited *pedited)
{
    pp->drcomp.threshold = threshold->getValue();
    pp->drcomp.amount = amount->getValue();
    pp->drcomp.anchor = anchor->getValue();
    pp->drcomp.enabled = getEnabled();

    if (pedited) {
        pedited->drcomp.threshold = threshold->getEditedState();
        pedited->drcomp.amount = amount->getEditedState();
        pedited->drcomp.anchor = anchor->getEditedState();
        pedited->drcomp.enabled = !get_inconsistent();
    }
}

void DRCompression::setDefaults(const ProcParams *defParams, const ParamsEdited *pedited)
{
    threshold->setDefault(defParams->drcomp.threshold);
    amount->setDefault(defParams->drcomp.amount);
    anchor->setDefault(defParams->drcomp.anchor);

    if (pedited) {
        threshold->setDefaultEditedState(pedited->drcomp.threshold ? Edited : UnEdited);
        amount->setDefaultEditedState(pedited->drcomp.amount ? Edited : UnEdited);
        anchor->setDefaultEditedState(pedited->drcomp.anchor ? Edited : UnEdited);
    } else {
        threshold->setDefaultEditedState(Irrelevant);
        amount->setDefaultEditedState(Irrelevant);
        anchor->setDefaultEditedState(Irrelevant);
    }
}

void DRCompression::adjusterChanged(Adjuster* a, double newval)
{
    if(listener && getEnabled()) {
        if(a == threshold) {
            listener->panelChanged(EvTMFattalThreshold, a->getTextValue());
        } else if(a == amount) {
            listener->panelChanged(EvTMFattalAmount, a->getTextValue());
        } else if(a == anchor) {
            listener->panelChanged(EvTMFattalAnchor, a->getTextValue());
        }
    }
}

void DRCompression::adjusterAutoToggled(Adjuster* a, bool newval)
{
}

void DRCompression::enabledChanged ()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvTMFattalEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvTMFattalEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvTMFattalEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void DRCompression::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);

    threshold->showEditedCB();
    amount->showEditedCB();
    anchor->showEditedCB();
}

void DRCompression::setAdjusterBehavior(bool amountAdd, bool thresholdAdd, bool anchorAdd)
{
    amount->setAddMode(amountAdd);
    threshold->setAddMode(thresholdAdd);
    anchor->setAddMode(anchorAdd);
}
