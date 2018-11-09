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

LogEncoding::LogEncoding(): FoldableToolPanel(this, "log", M("TP_LOGENC_LABEL"), true, true)
{
    auto m = ProcEventMapper::getInstance();
    EvEnabled = m->newEvent(RGBCURVE | M_AUTOEXP, "HISTORY_MSG_LOGENC_ENABLED");
    EvAuto = m->newEvent(AUTOEXP, "HISTORY_MSG_LOGENC_AUTO");
    EvAutoBatch = m->newEvent(M_VOID, "HISTORY_MSG_LOGENC_AUTO");
    EvDynamicRange = m->newEvent(RGBCURVE, "HISTORY_MSG_LOGENC_DYNAMIC_RANGE");
    EvGrayPoint = m->newEvent(RGBCURVE, "HISTORY_MSG_LOGENC_GRAY_POINT");
    EvShadowsRange = m->newEvent(RGBCURVE, "HISTORY_MSG_LOGENC_SHADOWS_RANGE");

    autocompute = Gtk::manage(new Gtk::ToggleButton(M("TP_LOGENC_AUTO")));
    autoconn = autocompute->signal_toggled().connect(sigc::mem_fun(*this, &LogEncoding::autocomputeToggled));
    
    dynamicRange = Gtk::manage(new Adjuster(M("TP_LOGENC_DYNAMIC_RANGE"), 1.0, 32.0, 0.1, 10.0));
    grayPoint = Gtk::manage(new Adjuster(M("TP_LOGENC_GRAY_POINT"), 1.0, 100.0, 0.1, 18.0));
    shadowsRange = Gtk::manage(new Adjuster(M("TP_LOGENC_SHADOWS_RANGE"), -16.0, 0.0, 0.1, -5.0));

    dynamicRange->setAdjusterListener(this);
    grayPoint->setAdjusterListener(this);
    shadowsRange->setAdjusterListener(this);

    autocompute->show();
    dynamicRange->show();
    grayPoint->show();
    shadowsRange->show();

    pack_start(*autocompute);
    pack_start(*grayPoint);
    pack_start(*shadowsRange);
    pack_start(*dynamicRange);
}


void LogEncoding::read(const ProcParams *pp, const ParamsEdited *pedited)
{
    disableListener();
    ConnectionBlocker cbl(autoconn);

    if (pedited) {
        dynamicRange->setEditedState(pedited->logenc.dynamicRange ? Edited : UnEdited);
        grayPoint->setEditedState(pedited->logenc.grayPoint ? Edited : UnEdited);
        shadowsRange->setEditedState(pedited->logenc.shadowsRange ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->logenc.enabled);
        autocompute->set_inconsistent(!pedited->logenc.autocompute);
    }

    setEnabled(pp->logenc.enabled);

    autocompute->set_active(pp->logenc.autocompute);    
    dynamicRange->setValue(pp->logenc.dynamicRange);
    grayPoint->setValue(pp->logenc.grayPoint);
    shadowsRange->setValue(pp->logenc.shadowsRange);

    enableListener();
}

void LogEncoding::write(ProcParams *pp, ParamsEdited *pedited)
{
    pp->logenc.enabled = getEnabled();
    pp->logenc.autocompute = autocompute->get_active();
    pp->logenc.dynamicRange = dynamicRange->getValue();
    pp->logenc.grayPoint = grayPoint->getValue();
    pp->logenc.shadowsRange = shadowsRange->getValue();

    if (pedited) {
        pedited->logenc.enabled = !get_inconsistent();
        pedited->logenc.autocompute = !autocompute->get_inconsistent();
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
    ConnectionBlocker cbl(autoconn);
    autocompute->set_active(false);
    
    if (listener && getEnabled()) {
        if (a == dynamicRange) {
            listener->panelChanged(EvDynamicRange, a->getTextValue());
        } else if (a == grayPoint) {
            listener->panelChanged(EvGrayPoint, a->getTextValue());
        } else if (a == shadowsRange) {
            listener->panelChanged(EvShadowsRange, a->getTextValue());
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


void LogEncoding::autocomputeToggled()
{
    if (listener) {
        if (!batchMode) {
            if (autocompute->get_active()) {
                listener->panelChanged(EvAuto, M("GENERAL_ENABLED"));
                dynamicRange->setEnabled(false);
                grayPoint->setEnabled(false);
                shadowsRange->setEnabled(false);
            } else {
                listener->panelChanged(EvAuto, M("GENERAL_DISABLED"));
            }
        } else {
            listener->panelChanged(EvAutoBatch, autocompute->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
        }
    }
}


void LogEncoding::logEncodingChanged(const rtengine::LogEncodingParams &params)
{
    GThreadLock lock;
    
    disableListener();
    ConnectionBlocker cbl(autoconn);

    dynamicRange->setEnabled(true);
    grayPoint->setEnabled(true);
    shadowsRange->setEnabled(true);

    dynamicRange->setValue(params.dynamicRange);
    grayPoint->setValue(params.grayPoint);
    shadowsRange->setValue(params.shadowsRange);
    
    enableListener();
}
