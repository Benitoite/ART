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
    EvAutoBatch = m->newEvent(M_VOID, "HISTORY_MSG_LOGENC_AUTO");
    EvGrayPoint = m->newEvent(M_LUMINANCE, "HISTORY_MSG_LOGENC_GRAY_POINT");
    EvBlackEv = m->newEvent(M_LUMINANCE, "HISTORY_MSG_LOGENC_BLACK_EV");
    EvWhiteEv = m->newEvent(M_LUMINANCE, "HISTORY_MSG_LOGENC_WHITE_EV");
    EvBase = m->newEvent(M_LUMINANCE, "HISTORY_MSG_LOGENC_BASE");

    autocompute = Gtk::manage(new Gtk::ToggleButton(M("TP_LOGENC_AUTO")));
    autoconn = autocompute->signal_toggled().connect(sigc::mem_fun(*this, &LogEncoding::autocomputeToggled));
    
    grayPoint = Gtk::manage(new Adjuster(M("TP_LOGENC_GRAY_POINT"), 1.0, 100.0, 0.1, 18.0));
    blackEv = Gtk::manage(new Adjuster(M("TP_LOGENC_BLACK_EV"), -16.0, 0.0, 0.1, -5.0));
    whiteEv = Gtk::manage(new Adjuster(M("TP_LOGENC_WHITE_EV"), 0.0, 32.0, 0.1, 10.0));
    base = Gtk::manage(new Adjuster(M("TP_LOGENC_BASE"), 0.0, 100.0, 0.1, 3.9));
    base->setLogScale(10, 0);

    whiteEv->setAdjusterListener(this);
    grayPoint->setAdjusterListener(this);
    blackEv->setAdjusterListener(this);
    base->setAdjusterListener(this);

    autocompute->show();
    whiteEv->show();
    grayPoint->show();
    blackEv->show();
    base->show();

    pack_start(*autocompute);
    pack_start(*grayPoint);
    pack_start(*blackEv);
    pack_start(*whiteEv);
    pack_start(*base);
}


void LogEncoding::read(const ProcParams *pp, const ParamsEdited *pedited)
{
    disableListener();
    ConnectionBlocker cbl(autoconn);

    if (pedited) {
        grayPoint->setEditedState(pedited->logenc.grayPoint ? Edited : UnEdited);
        blackEv->setEditedState(pedited->logenc.blackEv ? Edited : UnEdited);
        whiteEv->setEditedState(pedited->logenc.whiteEv ? Edited : UnEdited);
        base->setEditedState(pedited->logenc.base ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->logenc.enabled);
        autocompute->set_inconsistent(!pedited->logenc.autocompute);
    }

    setEnabled(pp->logenc.enabled);

    autocompute->set_active(pp->logenc.autocompute);    
    grayPoint->setValue(pp->logenc.grayPoint);
    blackEv->setValue(pp->logenc.blackEv);
    whiteEv->setValue(pp->logenc.whiteEv);
    base->setValue(pp->logenc.base);

    enableListener();
}

void LogEncoding::write(ProcParams *pp, ParamsEdited *pedited)
{
    pp->logenc.enabled = getEnabled();
    pp->logenc.autocompute = autocompute->get_active();
    pp->logenc.grayPoint = grayPoint->getValue();
    pp->logenc.blackEv = blackEv->getValue();
    pp->logenc.whiteEv = whiteEv->getValue();
    pp->logenc.base = base->getValue();

    if (pedited) {
        pedited->logenc.enabled = !get_inconsistent();
        pedited->logenc.autocompute = !autocompute->get_inconsistent();
        pedited->logenc.grayPoint = grayPoint->getEditedState();
        pedited->logenc.blackEv = blackEv->getEditedState();
        pedited->logenc.whiteEv = whiteEv->getEditedState();
        pedited->logenc.base = base->getEditedState();
    }
}

void LogEncoding::setDefaults(const ProcParams *defParams, const ParamsEdited *pedited)
{
    grayPoint->setDefault(defParams->logenc.grayPoint);
    blackEv->setDefault(defParams->logenc.blackEv);
    whiteEv->setDefault(defParams->logenc.whiteEv);
    base->setDefault(defParams->logenc.base);
    
    if (pedited) {
        grayPoint->setDefaultEditedState(pedited->logenc.grayPoint ? Edited : UnEdited);
        blackEv->setDefaultEditedState(pedited->logenc.blackEv ? Edited : UnEdited);
        whiteEv->setDefaultEditedState(pedited->logenc.whiteEv ? Edited : UnEdited);
        base->setDefaultEditedState(pedited->logenc.base ? Edited : UnEdited);
    } else {
        grayPoint->setDefaultEditedState(Irrelevant);
        blackEv->setDefaultEditedState(Irrelevant);
        whiteEv->setDefaultEditedState(Irrelevant);
        base->setDefaultEditedState(Irrelevant);
    }
}

void LogEncoding::adjusterChanged(Adjuster* a, double newval)
{
    ConnectionBlocker cbl(autoconn);
    autocompute->set_active(false);
    
    if (listener && getEnabled()) {
        if (a == grayPoint) {
            listener->panelChanged(EvGrayPoint, a->getTextValue());
        } else if (a == blackEv) {
            listener->panelChanged(EvBlackEv, a->getTextValue());
        } else if (a == whiteEv) {
            listener->panelChanged(EvWhiteEv, a->getTextValue());
        } else if (a == base) {
            listener->panelChanged(EvBase, a->getTextValue());
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

    grayPoint->showEditedCB();
    blackEv->showEditedCB();
    whiteEv->showEditedCB();
    base->showEditedCB();
}


void LogEncoding::autocomputeToggled()
{
    if (listener) {
        if (!batchMode) {
            if (autocompute->get_active()) {
                listener->panelChanged(EvAuto, M("GENERAL_ENABLED"));
                grayPoint->setEnabled(false);
                blackEv->setEnabled(false);
                whiteEv->setEnabled(false);
                base->setEnabled(false);
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

    grayPoint->setEnabled(true);
    blackEv->setEnabled(true);
    whiteEv->setEnabled(true);
    base->setEnabled(true);

    grayPoint->setValue(params.grayPoint);
    blackEv->setValue(params.blackEv);
    whiteEv->setValue(params.whiteEv);
    base->setValue(params.base);
    
    enableListener();
}
