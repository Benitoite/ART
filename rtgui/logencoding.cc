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
#include "logencoding.h"
#include "eventmapper.h"
#include <iomanip>
#include <cmath>

using namespace rtengine;
using namespace rtengine::procparams;

LogEncoding::LogEncoding(): FoldableToolPanel(this, "log", M("TP_TM_LOG_LABEL"), true, true)
{
    auto m = ProcEventMapper::getInstance();
    EvEnabled = m->newEvent(RGBCURVE, "HISTORY_MSG_DR_COMP_LOG_ENABLED");
    EvDRLogDynamicRange = m->newEvent(RGBCURVE, "HISTORY_MSG_DR_COMP_LOG_DYNAMIC_RANGE");
    EvDRLogGrayPoint = m->newEvent(RGBCURVE, "HISTORY_MSG_DR_COMP_LOG_GRAY_POINT");
    EvDRLogShadowsRange = m->newEvent(RGBCURVE, "HISTORY_MSG_DR_COMP_LOG_SHADOWS_RANGE");

    dynamicRange = Gtk::manage(new Adjuster(M("TP_DR_COMP_LOG_DYNAMIC_RANGE"), 1.0, 32.0, 0.1, 10.0));
    grayPoint = Gtk::manage(new Adjuster(M("TP_DR_COMP_LOG_GRAY_POINT"), 1.0, 100.0, 0.1, 18.0));
    shadowsRange = Gtk::manage(new Adjuster(M("TP_DR_COMP_LOG_SHADOWS_RANGE"), -16.0, 0.0, 0.1, -5.0));

    dynamicRange->setAdjusterListener(this);
    grayPoint->setAdjusterListener(this);
    shadowsRange->setAdjusterListener(this);

    dynamicRange->show();
    grayPoint->show();
    shadowsRange->show();

    pack_start(*grayPoint);
    pack_start(*shadowsRange);
    pack_start(*dynamicRange);
}


void LogEncoding::read(const ProcParams *pp, const ParamsEdited *pedited)
{
    disableListener();

    if (pedited) {
        dynamicRange->setEditedState(pedited->logenc.dynamicRange ? Edited : UnEdited);
        grayPoint->setEditedState(pedited->logenc.grayPoint ? Edited : UnEdited);
        shadowsRange->setEditedState(pedited->logenc.shadowsRange ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->logenc.enabled);
    }

    setEnabled(pp->logenc.enabled);

    dynamicRange->setValue(pp->logenc.dynamicRange);
    grayPoint->setValue(pp->logenc.grayPoint);
    shadowsRange->setValue(pp->logenc.shadowsRange);

    enableListener();
}

void LogEncoding::write(ProcParams *pp, ParamsEdited *pedited)
{
    pp->logenc.enabled = getEnabled();
    pp->logenc.dynamicRange = dynamicRange->getValue();
    pp->logenc.grayPoint = grayPoint->getValue();
    pp->logenc.shadowsRange = shadowsRange->getValue();

    if (pedited) {
        pedited->logenc.enabled = !get_inconsistent();
        pedited->logenc.dynamicRange = dynamicRange->getEditedState();
        pedited->logenc.grayPoint = grayPoint->getEditedState();
        pedited->logenc.shadowsRange = shadowsRange->getEditedState();
    }
}

void LogEncoding::setDefaults(const ProcParams *defParams, const ParamsEdited *pedited)
{
    dynamicRange->setDefault(defParams->logenc.dynamicRange);
    grayPoint->setDefault(defParams->logenc.grayPoint);
    shadowsRange->setDefault(defParams->logenc.shadowsRange);
    
    if (pedited) {
        dynamicRange->setDefaultEditedState(pedited->logenc.dynamicRange ? Edited : UnEdited);
        grayPoint->setDefaultEditedState(pedited->logenc.grayPoint ? Edited : UnEdited);
        shadowsRange->setDefaultEditedState(pedited->logenc.shadowsRange ? Edited : UnEdited);
    } else {
        dynamicRange->setDefaultEditedState(Irrelevant);
        grayPoint->setDefaultEditedState(Irrelevant);
        shadowsRange->setDefaultEditedState(Irrelevant);
    }
}

void LogEncoding::adjusterChanged(Adjuster* a, double newval)
{
    if(listener && getEnabled()) {
        if (a == dynamicRange) {
            listener->panelChanged(EvDRLogDynamicRange, a->getTextValue());
        } else if (a == grayPoint) {
            listener->panelChanged(EvDRLogGrayPoint, a->getTextValue());
        } else if (a == shadowsRange) {
            listener->panelChanged(EvDRLogShadowsRange, a->getTextValue());
        }
    }
}

void LogEncoding::adjusterAutoToggled(Adjuster* a, bool newval)
{
}

void LogEncoding::enabledChanged ()
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

void LogEncoding::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);

    dynamicRange->showEditedCB();
    grayPoint->showEditedCB();
    shadowsRange->showEditedCB();
}
