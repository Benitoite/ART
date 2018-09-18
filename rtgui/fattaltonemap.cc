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
#include "fattaltonemap.h"
#include "eventmapper.h"
#include <iomanip>
#include <cmath>

using namespace rtengine;
using namespace rtengine::procparams;

FattalToneMapping::FattalToneMapping(): FoldableToolPanel(this, "fattal", M("TP_TM_FATTAL_LABEL"), true, true)
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
    offset = Gtk::manage(new Adjuster (M("TP_DR_COMP_OFFSET"), -0.002, 0.002, 0.0001, 0.0));

    power->setAdjusterListener(this);
    slope->setAdjusterListener(this);
    offset->setAdjusterListener(this);

    gammabox->show();
    power->show();
    slope->show();
    offset->show();

    gammabox->pack_start(*power);
    gammabox->pack_start(*slope);
    gammabox->pack_start(*offset);

    pack_start(*gammabox);

    method->signal_changed().connect(sigc::mem_fun(*this, &FattalToneMapping::method_changed));
}

void FattalToneMapping::read(const ProcParams *pp, const ParamsEdited *pedited)
{
    disableListener();

    if (pedited) {
        threshold->setEditedState(pedited->fattal.threshold ? Edited : UnEdited);
        amount->setEditedState(pedited->fattal.amount ? Edited : UnEdited);
        anchor->setEditedState(pedited->fattal.anchor ? Edited : UnEdited);
        power->setEditedState(pedited->fattal.power ? Edited : UnEdited);
        slope->setEditedState(pedited->fattal.slope ? Edited : UnEdited);
        offset->setEditedState(pedited->fattal.offset ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->fattal.enabled);
    }

    setEnabled(pp->fattal.enabled);
    threshold->setValue(pp->fattal.threshold);
    amount->setValue(pp->fattal.amount);
    anchor->setValue(pp->fattal.anchor);

    power->setValue(pp->fattal.power);
    slope->setValue(pp->fattal.slope);
    offset->setValue(pp->fattal.offset);

    if (pedited && !pedited->fattal.method) {
        method->set_active(2);
    } else {
        method->set_active(pp->fattal.method);
        if (pp->fattal.method == rtengine::procparams::FattalToneMappingParams::DR_COMP_FATTAL) {
            fattalbox->show();
            gammabox->hide();
        } else {
            fattalbox->hide();
            gammabox->show();
        }
    }
    
    enableListener();
}

void FattalToneMapping::write(ProcParams *pp, ParamsEdited *pedited)
{
    pp->fattal.threshold = threshold->getValue();
    pp->fattal.amount = amount->getValue();
    pp->fattal.anchor = anchor->getValue();
    pp->fattal.enabled = getEnabled();
    pp->fattal.power = power->getValue();
    pp->fattal.slope = slope->getValue();
    pp->fattal.offset = offset->getValue();
    if (method->get_active_row_number() != 2) {
        pp->fattal.method = method->get_active_row_number();
    }

    if(pedited) {
        pedited->fattal.threshold = threshold->getEditedState();
        pedited->fattal.amount = amount->getEditedState();
        pedited->fattal.anchor = anchor->getEditedState();
        pedited->fattal.enabled = !get_inconsistent();
        pedited->fattal.power = power->getEditedState();
        pedited->fattal.slope = slope->getEditedState();
        pedited->fattal.offset = offset->getEditedState();
        pedited->fattal.method = method->get_active_row_number() != 2;
    }
}

void FattalToneMapping::setDefaults(const ProcParams *defParams, const ParamsEdited *pedited)
{
    threshold->setDefault(defParams->fattal.threshold);
    amount->setDefault(defParams->fattal.amount);
    anchor->setDefault(defParams->fattal.anchor);

    power->setDefault(defParams->fattal.power);
    slope->setDefault(defParams->fattal.slope);
    offset->setDefault(defParams->fattal.offset);
    
    if(pedited) {
        threshold->setDefaultEditedState(pedited->fattal.threshold ? Edited : UnEdited);
        amount->setDefaultEditedState(pedited->fattal.amount ? Edited : UnEdited);
        anchor->setDefaultEditedState(pedited->fattal.anchor ? Edited : UnEdited);

        power->setDefaultEditedState(pedited->fattal.power ? Edited : UnEdited);
        slope->setDefaultEditedState(pedited->fattal.slope ? Edited : UnEdited);
        offset->setDefaultEditedState(pedited->fattal.offset ? Edited : UnEdited);
    } else {
        threshold->setDefaultEditedState(Irrelevant);
        amount->setDefaultEditedState(Irrelevant);
        anchor->setDefaultEditedState(Irrelevant);

        power->setDefaultEditedState(Irrelevant);
        slope->setDefaultEditedState(Irrelevant);
        offset->setDefaultEditedState(Irrelevant);
    }
}

void FattalToneMapping::adjusterChanged(Adjuster* a, double newval)
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
        } else if (a == offset) {
            listener->panelChanged(EvDRCompOffset, a->getTextValue());
        }
    }
}

void FattalToneMapping::enabledChanged ()
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

void FattalToneMapping::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);

    threshold->showEditedCB();
    amount->showEditedCB();
    anchor->showEditedCB();

    power->showEditedCB();
    slope->showEditedCB();
    offset->showEditedCB();
    
    method->append(M("GENERAL_UNCHANGED"));
}

void FattalToneMapping::setAdjusterBehavior(bool amountAdd, bool thresholdAdd, bool anchorAdd)
{
    amount->setAddMode(amountAdd);
    threshold->setAddMode(thresholdAdd);
    anchor->setAddMode(anchorAdd);
}


void FattalToneMapping::method_changed()
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
