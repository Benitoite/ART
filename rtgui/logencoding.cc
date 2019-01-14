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

LogEncoding::LogEncoding(): FoldableToolPanel(this, "log", M("TP_LOGENC_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    EvEnabled = m->newEvent(RGBCURVE | M_AUTOEXP, "HISTORY_MSG_LOGENC_ENABLED");
    EvAuto = m->newEvent(AUTOEXP, "HISTORY_MSG_LOGENC_AUTO");
    EvAutoGrayOn = m->newEvent(AUTOEXP, "HISTORY_MSG_LOGENC_AUTOGRAY");
    EvAutoGrayOff = m->newEvent(M_VOID, "HISTORY_MSG_LOGENC_AUTOGRAY");
    EvAutoBatch = m->newEvent(M_VOID, "HISTORY_MSG_LOGENC_AUTO");
    EvSourceGray = m->newEvent(M_LUMINANCE, "HISTORY_MSG_LOGENC_SOURCE_GRAY");
    EvSourceGrayAuto = m->newEvent(M_LUMINANCE | M_AUTOEXP, "HISTORY_MSG_LOGENC_GRAY_POINT");
    EvTargetGray = m->newEvent(M_LUMINANCE, "HISTORY_MSG_LOGENC_TARGET_GRAY");
    EvBlackEv = m->newEvent(M_LUMINANCE, "HISTORY_MSG_LOGENC_BLACK_EV");
    EvWhiteEv = m->newEvent(M_LUMINANCE, "HISTORY_MSG_LOGENC_WHITE_EV");

    autocompute = Gtk::manage(new Gtk::ToggleButton(M("TP_LOGENC_AUTO")));
    autoconn = autocompute->signal_toggled().connect(sigc::mem_fun(*this, &LogEncoding::autocomputeToggled));
    
    sourceGray = Gtk::manage(new Adjuster(M("TP_LOGENC_SOURCE_GRAY"), 1.0, 100.0, 0.1, 18.0));
    sourceGray->addAutoButton();
    targetGray = Gtk::manage(new Adjuster(M("TP_LOGENC_TARGET_GRAY"), 5.0, 80.0, 0.1, 18.0));
    blackEv = Gtk::manage(new Adjuster(M("TP_LOGENC_BLACK_EV"), -16.0, 0.0, 0.1, -5.0));
    whiteEv = Gtk::manage(new Adjuster(M("TP_LOGENC_WHITE_EV"), 0.0, 32.0, 0.1, 10.0));

    sourceGray->delay = options.adjusterMaxDelay;
    blackEv->delay = options.adjusterMaxDelay;
    whiteEv->delay = options.adjusterMaxDelay;
    targetGray->delay = options.adjusterMaxDelay;
    
    whiteEv->setAdjusterListener(this);
    sourceGray->setAdjusterListener(this);
    blackEv->setAdjusterListener(this);
    targetGray->setAdjusterListener(this);

    whiteEv->setLogScale(16, 0);
    blackEv->setLogScale(2, -8);

    sourceGray->show();
    autocompute->show();
    whiteEv->show();
    blackEv->show();
    targetGray->show();

    pack_start(*sourceGray);
    pack_start(*targetGray);
    pack_start(*autocompute);
    pack_start(*blackEv);
    pack_start(*whiteEv);
}


void LogEncoding::read(const ProcParams *pp, const ParamsEdited *pedited)
{
    disableListener();
    ConnectionBlocker cbl(autoconn);

    if (pedited) {
        sourceGray->setEditedState(pedited->logenc.sourceGray ? Edited : UnEdited);
        blackEv->setEditedState(pedited->logenc.blackEv ? Edited : UnEdited);
        whiteEv->setEditedState(pedited->logenc.whiteEv ? Edited : UnEdited);
        targetGray->setEditedState(pedited->logenc.targetGray ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->logenc.enabled);
        autocompute->set_inconsistent(!pedited->logenc.autocompute);
    }

    setEnabled(pp->logenc.enabled);

    autocompute->set_active(pp->logenc.autocompute);    
    sourceGray->setValue(pp->logenc.sourceGray);
    sourceGray->setAutoValue(pp->logenc.autogray);
    blackEv->setValue(pp->logenc.blackEv);
    whiteEv->setValue(pp->logenc.whiteEv);
    targetGray->setValue(pp->logenc.targetGray);

    enableListener();
}

void LogEncoding::write(ProcParams *pp, ParamsEdited *pedited)
{
    pp->logenc.enabled = getEnabled();
    pp->logenc.autocompute = autocompute->get_active();
    pp->logenc.autogray = sourceGray->getAutoValue();
    pp->logenc.sourceGray = sourceGray->getValue();
    pp->logenc.blackEv = blackEv->getValue();
    pp->logenc.whiteEv = whiteEv->getValue();
    pp->logenc.targetGray = targetGray->getValue();

    if (pedited) {
        pedited->logenc.enabled = !get_inconsistent();
        pedited->logenc.autocompute = !autocompute->get_inconsistent();
        pedited->logenc.sourceGray = sourceGray->getEditedState();
        pedited->logenc.blackEv = blackEv->getEditedState();
        pedited->logenc.whiteEv = whiteEv->getEditedState();
        pedited->logenc.targetGray = targetGray->getEditedState();
    }
}

void LogEncoding::setDefaults(const ProcParams *defParams, const ParamsEdited *pedited)
{
    sourceGray->setDefault(defParams->logenc.sourceGray);
    blackEv->setDefault(defParams->logenc.blackEv);
    whiteEv->setDefault(defParams->logenc.whiteEv);
    targetGray->setDefault(defParams->logenc.targetGray);
    
    if (pedited) {
        sourceGray->setDefaultEditedState(pedited->logenc.sourceGray ? Edited : UnEdited);
        blackEv->setDefaultEditedState(pedited->logenc.blackEv ? Edited : UnEdited);
        whiteEv->setDefaultEditedState(pedited->logenc.whiteEv ? Edited : UnEdited);
        targetGray->setDefaultEditedState(pedited->logenc.targetGray ? Edited : UnEdited);
    } else {
        sourceGray->setDefaultEditedState(Irrelevant);
        blackEv->setDefaultEditedState(Irrelevant);
        whiteEv->setDefaultEditedState(Irrelevant);
        targetGray->setDefaultEditedState(Irrelevant);
    }
}

void LogEncoding::adjusterChanged(Adjuster* a, double newval)
{
    ConnectionBlocker cbl(autoconn);
    if (a != sourceGray && a != targetGray) {
        autocompute->set_active(false);
    }
    
    if (listener && getEnabled()) {
        if (a == sourceGray) {
            listener->panelChanged(autocompute->get_active() ? EvSourceGrayAuto : EvSourceGray, a->getTextValue());
        } else if (a == blackEv) {
            listener->panelChanged(EvBlackEv, a->getTextValue());
        } else if (a == whiteEv) {
            listener->panelChanged(EvWhiteEv, a->getTextValue());
        } else if (a == targetGray) {
            listener->panelChanged(EvTargetGray, a->getTextValue());
        }
    }
}

void LogEncoding::adjusterAutoToggled(Adjuster* a, bool newval)
{
    if (listener) {
        if (a == sourceGray) {
            auto e = (batchMode || !newval) ? EvAutoGrayOff : EvAutoGrayOn;
            listener->panelChanged(e, newval ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
        }
    }
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

    sourceGray->showEditedCB();
    blackEv->showEditedCB();
    whiteEv->showEditedCB();
    targetGray->showEditedCB();
}


void LogEncoding::autocomputeToggled()
{
    if (listener) {
        if (!batchMode) {
            if (autocompute->get_active()) {
                listener->panelChanged(EvAuto, M("GENERAL_ENABLED"));
                blackEv->setEnabled(false);
                whiteEv->setEnabled(false);
                //targetGray->setEnabled(false);
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

    blackEv->setEnabled(true);
    whiteEv->setEnabled(true);
//    targetGray->setEnabled(true);

    sourceGray->setValue(params.sourceGray);
    blackEv->setValue(params.blackEv);
    whiteEv->setValue(params.whiteEv);
//    targetGray->setValue(params.targetGray);
    
    enableListener();
}
