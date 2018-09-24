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
    EvDRCompMethod = m->newEvent(HDR, "HISTORY_MSG_DR_COMP_METHOD");
    EvDRCompPower = m->newEvent(HDR, "HISTORY_MSG_DR_COMP_POWER");
    EvDRCompSlope = m->newEvent(HDR, "HISTORY_MSG_DR_COMP_SLOPE");
    EvDRCompOffset = m->newEvent(HDR, "HISTORY_MSG_DR_COMP_OFFSET");

    Gtk::HBox *b = Gtk::manage(new Gtk::HBox());
    b->pack_start(*Gtk::manage(new Gtk::Label(M("TP_SHARPENING_METHOD") + ":")), Gtk::PACK_SHRINK, 4);
    method = Gtk::manage(new MyComboBoxText());
    method->append(M("TP_DR_COMP_METHOD_FATTAL"));
    method->append(M("TP_DR_COMP_METHOD_GAMMA"));

    method->show();
    b->pack_start(*method);
    pack_start(*b);

    fattalbox = Gtk::manage(new Gtk::VBox());
    amount = Gtk::manage(new Adjuster (M("TP_TM_FATTAL_AMOUNT"), 1., 100., 1., 30.));
    threshold = Gtk::manage(new Adjuster (M("TP_TM_FATTAL_THRESHOLD"), -100., 100., 1., 0.0));
    Gtk::Image *al = Gtk::manage(new RTImage("circle-black-small.png"));
    Gtk::Image *ar = Gtk::manage(new RTImage("circle-white-small.png"));
    anchor = Gtk::manage(new Adjuster(M("TP_TM_FATTAL_ANCHOR"), 1, 100, 1, 50, al, ar));

    amount->setAdjusterListener(this);
    threshold->setAdjusterListener(this);
    anchor->setAdjusterListener(this);

    fattalbox->show();
    amount->show();
    threshold->show();
    anchor->show();

    fattalbox->pack_start(*amount);
    fattalbox->pack_start(*threshold);
    fattalbox->pack_start(*anchor);

    pack_start(*fattalbox);

    gammabox = Gtk::manage(new Gtk::VBox());
    power = Gtk::manage(new Adjuster (M("TP_DR_COMP_POWER"), 0.0, 2.0, 0.005, 1.));
    slope = Gtk::manage(new Adjuster (M("TP_DR_COMP_SLOPE"), 0.001, 2.0, 0.005, 1.0));

    power->setAdjusterListener(this);
    slope->setAdjusterListener(this);

    gammabox->show();
    power->show();
    slope->show();

    gammabox->pack_start(*power);
    gammabox->pack_start(*slope);

    pack_start(*gammabox);

    method->signal_changed().connect(sigc::mem_fun(*this, &DRCompression::method_changed));
}

void DRCompression::read(const ProcParams *pp, const ParamsEdited *pedited)
{
    disableListener();

    if (pedited) {
        threshold->setEditedState(pedited->drcomp.threshold ? Edited : UnEdited);
        amount->setEditedState(pedited->drcomp.amount ? Edited : UnEdited);
        anchor->setEditedState(pedited->drcomp.anchor ? Edited : UnEdited);
        power->setEditedState(pedited->drcomp.power ? Edited : UnEdited);
        slope->setEditedState(pedited->drcomp.slope ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->drcomp.enabled);
    }

    setEnabled(pp->drcomp.enabled);
    threshold->setValue(pp->drcomp.threshold);
    amount->setValue(pp->drcomp.amount);
    anchor->setValue(pp->drcomp.anchor);

    power->setValue(pp->drcomp.power);
    slope->setValue(pp->drcomp.slope);

    if (pedited && !pedited->drcomp.method) {
        method->set_active(2);
    } else {
        method->set_active(pp->drcomp.method);
        if (pp->drcomp.method == rtengine::procparams::DRCompressionParams::DR_COMP_FATTAL) {
            fattalbox->show();
            gammabox->hide();
        } else {
            fattalbox->hide();
            gammabox->show();
        }
    }
    
    enableListener();
}

void DRCompression::write(ProcParams *pp, ParamsEdited *pedited)
{
    pp->drcomp.threshold = threshold->getValue();
    pp->drcomp.amount = amount->getValue();
    pp->drcomp.anchor = anchor->getValue();
    pp->drcomp.enabled = getEnabled();
    pp->drcomp.power = power->getValue();
    pp->drcomp.slope = slope->getValue();
    if (method->get_active_row_number() != 2) {
        pp->drcomp.method = method->get_active_row_number();
    }

    if(pedited) {
        pedited->drcomp.threshold = threshold->getEditedState();
        pedited->drcomp.amount = amount->getEditedState();
        pedited->drcomp.anchor = anchor->getEditedState();
        pedited->drcomp.enabled = !get_inconsistent();
        pedited->drcomp.power = power->getEditedState();
        pedited->drcomp.slope = slope->getEditedState();
        pedited->drcomp.method = method->get_active_row_number() != 2;
    }
}

void DRCompression::setDefaults(const ProcParams *defParams, const ParamsEdited *pedited)
{
    threshold->setDefault(defParams->drcomp.threshold);
    amount->setDefault(defParams->drcomp.amount);
    anchor->setDefault(defParams->drcomp.anchor);

    power->setDefault(defParams->drcomp.power);
    slope->setDefault(defParams->drcomp.slope);
    
    if(pedited) {
        threshold->setDefaultEditedState(pedited->drcomp.threshold ? Edited : UnEdited);
        amount->setDefaultEditedState(pedited->drcomp.amount ? Edited : UnEdited);
        anchor->setDefaultEditedState(pedited->drcomp.anchor ? Edited : UnEdited);

        power->setDefaultEditedState(pedited->drcomp.power ? Edited : UnEdited);
        slope->setDefaultEditedState(pedited->drcomp.slope ? Edited : UnEdited);
    } else {
        threshold->setDefaultEditedState(Irrelevant);
        amount->setDefaultEditedState(Irrelevant);
        anchor->setDefaultEditedState(Irrelevant);

        power->setDefaultEditedState(Irrelevant);
        slope->setDefaultEditedState(Irrelevant);
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
        } else if (a == power) {
            listener->panelChanged(EvDRCompPower, a->getTextValue());
        } else if (a == slope) {
            listener->panelChanged(EvDRCompSlope, a->getTextValue());
        }
    }
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

    power->showEditedCB();
    slope->showEditedCB();
    
    method->append(M("GENERAL_UNCHANGED"));
}

void DRCompression::setAdjusterBehavior(bool amountAdd, bool thresholdAdd, bool anchorAdd)
{
    amount->setAddMode(amountAdd);
    threshold->setAddMode(thresholdAdd);
    anchor->setAddMode(anchorAdd);
}


void DRCompression::method_changed()
{
    if (!batchMode) {
        fattalbox->show();
        gammabox->show();

        if (method->get_active_row_number() == 0) {
            gammabox->hide();
        } else if (method->get_active_row_number() == 1) {
            fattalbox->hide();
        }
    }

    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvDRCompMethod, method->get_active_text());
    }
}
