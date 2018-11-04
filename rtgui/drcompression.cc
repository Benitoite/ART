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
    EvDRLogDynamicRange = m->newEvent(HDR, "HISTORY_MSG_DR_COMP_LOG_DYNAMIC_RANGE");
    EvDRLogGrayPoint = m->newEvent(HDR, "HISTORY_MSG_DR_COMP_LOG_GRAY_POINT");
    EvDRLogShadowsRange = m->newEvent(HDR, "HISTORY_MSG_DR_COMP_LOG_SHADOWS_RANGE");

    Gtk::HBox *b = Gtk::manage(new Gtk::HBox());
    b->pack_start(*Gtk::manage(new Gtk::Label(M("TP_SHARPENING_METHOD") + ":")), Gtk::PACK_SHRINK, 4);
    method = Gtk::manage(new MyComboBoxText());
    method->append(M("TP_DR_COMP_METHOD_FATTAL"));
    method->append(M("TP_DR_COMP_METHOD_LOG"));

    method->show();
    b->pack_start(*method);
    pack_start(*b);

    fattalbox = Gtk::manage(new Gtk::VBox());
    amount = Gtk::manage(new Adjuster (M("TP_TM_FATTAL_AMOUNT"), 1., 100., 1., 30.));
    threshold = Gtk::manage(new Adjuster (M("TP_TM_FATTAL_THRESHOLD"), -100., 300., 1., 0.0));
    Gtk::Image *al = Gtk::manage(new RTImage("circle-black-small.png"));
    Gtk::Image *ar = Gtk::manage(new RTImage("circle-white-small.png"));
    anchor = Gtk::manage(new Adjuster(M("TP_TM_FATTAL_ANCHOR"), 1, 100, 1, 50, al, ar));

    amount->setAdjusterListener(this);
    threshold->setAdjusterListener(this);
    anchor->setAdjusterListener(this);

    threshold->setLogScale(10, 0);

    fattalbox->show();
    amount->show();
    threshold->show();
    anchor->show();

    fattalbox->pack_start(*amount);
    fattalbox->pack_start(*threshold);
    fattalbox->pack_start(*anchor);

    pack_start(*fattalbox);

    logbox = Gtk::manage(new Gtk::VBox());
    dynamicRange = Gtk::manage(new Adjuster(M("TP_DR_COMP_LOG_DYNAMIC_RANGE"), 1.0, 32.0, 0.1, 16.0));
    grayPoint = Gtk::manage(new Adjuster(M("TP_DR_COMP_LOG_GRAY_POINT"), 1.0, 100.0, 0.1, 18.0));
    shadowsRange = Gtk::manage(new Adjuster(M("TP_DR_COMP_LOG_SHADOWS_RANGE"), -16.0, 0.0, 0.1, -12.0));

    dynamicRange->setAdjusterListener(this);
    grayPoint->setAdjusterListener(this);
    shadowsRange->setAdjusterListener(this);

    logbox->show();
    dynamicRange->show();
    grayPoint->show();
    shadowsRange->show();

    logbox->pack_start(*grayPoint);
    logbox->pack_start(*shadowsRange);
    logbox->pack_start(*dynamicRange);

    pack_start(*logbox);

    method->signal_changed().connect(sigc::mem_fun(*this, &DRCompression::method_changed));
}

void DRCompression::read(const ProcParams *pp, const ParamsEdited *pedited)
{
    disableListener();

    if (pedited) {
        threshold->setEditedState(pedited->drcomp.threshold ? Edited : UnEdited);
        amount->setEditedState(pedited->drcomp.amount ? Edited : UnEdited);
        anchor->setEditedState(pedited->drcomp.anchor ? Edited : UnEdited);
        dynamicRange->setEditedState(pedited->drcomp.dynamicRange ? Edited : UnEdited);
        grayPoint->setEditedState(pedited->drcomp.grayPoint ? Edited : UnEdited);
        shadowsRange->setEditedState(pedited->drcomp.shadowsRange ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->drcomp.enabled);
    }

    setEnabled(pp->drcomp.enabled);
    threshold->setValue(pp->drcomp.threshold);
    amount->setValue(pp->drcomp.amount);
    anchor->setValue(pp->drcomp.anchor);

    dynamicRange->setValue(pp->drcomp.dynamicRange);
    grayPoint->setValue(pp->drcomp.grayPoint);
    shadowsRange->setValue(pp->drcomp.shadowsRange);

    if (pedited && !pedited->drcomp.method) {
        method->set_active(2);
    } else {
        method->set_active(pp->drcomp.method);
        if (pp->drcomp.method == rtengine::procparams::DRCompressionParams::DR_COMP_FATTAL) {
            fattalbox->show();
            logbox->hide();
        } else {
            fattalbox->hide();
            logbox->show();
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
    pp->drcomp.dynamicRange = dynamicRange->getValue();
    pp->drcomp.grayPoint = grayPoint->getValue();
    pp->drcomp.shadowsRange = shadowsRange->getValue();
    if (method->get_active_row_number() != 2) {
        pp->drcomp.method = method->get_active_row_number();
    }

    if (pedited) {
        pedited->drcomp.threshold = threshold->getEditedState();
        pedited->drcomp.amount = amount->getEditedState();
        pedited->drcomp.anchor = anchor->getEditedState();
        pedited->drcomp.enabled = !get_inconsistent();
        pedited->drcomp.dynamicRange = dynamicRange->getEditedState();
        pedited->drcomp.grayPoint = grayPoint->getEditedState();
        pedited->drcomp.shadowsRange = shadowsRange->getEditedState();
        pedited->drcomp.method = method->get_active_row_number() != 2;
    }
}

void DRCompression::setDefaults(const ProcParams *defParams, const ParamsEdited *pedited)
{
    threshold->setDefault(defParams->drcomp.threshold);
    amount->setDefault(defParams->drcomp.amount);
    anchor->setDefault(defParams->drcomp.anchor);

    dynamicRange->setDefault(defParams->drcomp.dynamicRange);
    grayPoint->setDefault(defParams->drcomp.grayPoint);
    shadowsRange->setDefault(defParams->drcomp.shadowsRange);
    
    if (pedited) {
        threshold->setDefaultEditedState(pedited->drcomp.threshold ? Edited : UnEdited);
        amount->setDefaultEditedState(pedited->drcomp.amount ? Edited : UnEdited);
        anchor->setDefaultEditedState(pedited->drcomp.anchor ? Edited : UnEdited);

        dynamicRange->setDefaultEditedState(pedited->drcomp.dynamicRange ? Edited : UnEdited);
        grayPoint->setDefaultEditedState(pedited->drcomp.grayPoint ? Edited : UnEdited);
        shadowsRange->setDefaultEditedState(pedited->drcomp.shadowsRange ? Edited : UnEdited);
    } else {
        threshold->setDefaultEditedState(Irrelevant);
        amount->setDefaultEditedState(Irrelevant);
        anchor->setDefaultEditedState(Irrelevant);

        dynamicRange->setDefaultEditedState(Irrelevant);
        grayPoint->setDefaultEditedState(Irrelevant);
        shadowsRange->setDefaultEditedState(Irrelevant);
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
        } else if (a == dynamicRange) {
            listener->panelChanged(EvDRLogDynamicRange, a->getTextValue());
        } else if (a == grayPoint) {
            listener->panelChanged(EvDRLogGrayPoint, a->getTextValue());
        } else if (a == shadowsRange) {
            listener->panelChanged(EvDRLogShadowsRange, a->getTextValue());
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

    dynamicRange->showEditedCB();
    grayPoint->showEditedCB();
    shadowsRange->showEditedCB();
    
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
        logbox->show();

        if (method->get_active_row_number() == 0) {
            logbox->hide();
        } else if (method->get_active_row_number() == 1) {
            fattalbox->hide();
        }
    }

    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged (EvDRCompMethod, method->get_active_text());
    }
}
