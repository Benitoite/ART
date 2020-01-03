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
    const auto EVENT = LUMINANCECURVE;
    EvEnabled = m->newEvent(RGBCURVE | M_AUTOEXP, "HISTORY_MSG_LOGENC_ENABLED");
    EvAuto = m->newEvent(AUTOEXP, "HISTORY_MSG_LOGENC_AUTO");
    EvAutoGrayOn = m->newEvent(AUTOEXP, "HISTORY_MSG_LOGENC_AUTOGRAY");
    EvAutoGrayOff = m->newEvent(M_VOID, "HISTORY_MSG_LOGENC_AUTOGRAY");
    EvAutoBatch = m->newEvent(M_VOID, "HISTORY_MSG_LOGENC_AUTO");
    EvSourceGray = m->newEvent(EVENT, "HISTORY_MSG_LOGENC_SOURCE_GRAY");
    EvSourceGrayAuto = m->newEvent(AUTOEXP, "HISTORY_MSG_LOGENC_SOURCE_GRAY");
    EvTargetGray = m->newEvent(EVENT, "HISTORY_MSG_LOGENC_TARGET_GRAY");
    EvBlackEv = m->newEvent(EVENT, "HISTORY_MSG_LOGENC_BLACK_EV");
    EvWhiteEv = m->newEvent(EVENT, "HISTORY_MSG_LOGENC_WHITE_EV");
    EvRegularization = m->newEvent(EVENT, "HISTORY_MSG_LOGENC_REGULARIZATION");

    autocompute = Gtk::manage(new Gtk::ToggleButton(M("TP_LOGENC_AUTO")));
    autoconn = autocompute->signal_toggled().connect(sigc::mem_fun(*this, &LogEncoding::autocomputeToggled));
    
    sourceGray = Gtk::manage(new Adjuster(M("TP_LOGENC_SOURCE_GRAY"), 1.0, 100.0, 0.1, 18.0));
    sourceGray->addAutoButton();
    targetGray = Gtk::manage(new Adjuster(M("TP_LOGENC_TARGET_GRAY"), 5.0, 80.0, 0.1, 18.0));
    blackEv = Gtk::manage(new Adjuster(M("TP_LOGENC_BLACK_EV"), -16.0, 0.0, 0.1, -5.0));
    whiteEv = Gtk::manage(new Adjuster(M("TP_LOGENC_WHITE_EV"), 0.0, 32.0, 0.1, 10.0));
    regularization = Gtk::manage(new Adjuster(M("TP_LOGENC_REGULARIZATION"), 0, 100, 1, 65));

    Gtk::Frame *evFrame = Gtk::manage(new Gtk::Frame(M("TP_LOGENC_EV_LEVELS")));
    evFrame->set_label_align(0.025, 0.5);

    Gtk::VBox *evBox = Gtk::manage(new Gtk::VBox());
    evBox->set_border_width(2);
    evBox->set_spacing(2);

    sourceGray->delay = options.adjusterMaxDelay;
    blackEv->delay = options.adjusterMaxDelay;
    whiteEv->delay = options.adjusterMaxDelay;
    targetGray->delay = options.adjusterMaxDelay;
    
    whiteEv->setAdjusterListener(this);
    sourceGray->setAdjusterListener(this);
    blackEv->setAdjusterListener(this);
    targetGray->setAdjusterListener(this);
    regularization->setAdjusterListener(this);

    whiteEv->setLogScale(16, 0);
    blackEv->setLogScale(2, -8);

    sourceGray->show();
    autocompute->show();
    whiteEv->show();
    blackEv->show();
    targetGray->show();

    evBox->pack_start(*autocompute);
    evBox->pack_start(*blackEv);
    evBox->pack_start(*whiteEv);
    evFrame->add(*evBox);
    pack_start(*evFrame);
    pack_start(*sourceGray);
    pack_start(*targetGray);
    pack_start(*regularization);
}


void LogEncoding::read(const ProcParams *pp)
{
    disableListener();
    ConnectionBlocker cbl(autoconn);

    setEnabled(pp->logenc.enabled);

    autocompute->set_active(pp->logenc.autocompute);    
    sourceGray->setValue(pp->logenc.sourceGray);
    sourceGray->setAutoValue(pp->logenc.autogray);
    blackEv->setValue(pp->logenc.blackEv);
    whiteEv->setValue(pp->logenc.whiteEv);
    targetGray->setValue(pp->logenc.targetGray);
    regularization->setValue(pp->logenc.regularization);

    enableListener();
}

void LogEncoding::write(ProcParams *pp)
{
    pp->logenc.enabled = getEnabled();
    pp->logenc.autocompute = autocompute->get_active();
    pp->logenc.autogray = sourceGray->getAutoValue();
    pp->logenc.sourceGray = sourceGray->getValue();
    pp->logenc.blackEv = blackEv->getValue();
    pp->logenc.whiteEv = whiteEv->getValue();
    pp->logenc.targetGray = targetGray->getValue();
    pp->logenc.regularization = regularization->getValue();
}

void LogEncoding::setDefaults(const ProcParams *defParams)
{
    sourceGray->setDefault(defParams->logenc.sourceGray);
    blackEv->setDefault(defParams->logenc.blackEv);
    whiteEv->setDefault(defParams->logenc.whiteEv);
    targetGray->setDefault(defParams->logenc.targetGray);
    regularization->setDefault(defParams->logenc.regularization);
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
        } else if (a == regularization) {
            listener->panelChanged(EvRegularization, a->getTextValue());
        }
    }
}

void LogEncoding::adjusterAutoToggled(Adjuster* a, bool newval)
{
    if (listener) {
        if (a == sourceGray) {
            auto e = (!newval) ? EvAutoGrayOff : EvAutoGrayOn;
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


void LogEncoding::autocomputeToggled()
{
    if (listener) {
        if (autocompute->get_active()) {
            listener->panelChanged(EvAuto, M("GENERAL_ENABLED"));
            blackEv->setEnabled(false);
            whiteEv->setEnabled(false);
            //targetGray->setEnabled(false);
        } else {
            listener->panelChanged(EvAuto, M("GENERAL_DISABLED"));
            // blackEv->setEnabled(true);
            // whiteEv->setEnabled(true);
        }
    }
}


void LogEncoding::logEncodingChanged(const rtengine::LogEncodingParams &params)
{
    GThreadLock lock;
    
    disableListener();
    ConnectionBlocker cbl(autoconn);

    // blackEv->setEnabled(!params.autocompute);
    // whiteEv->setEnabled(!params.autocompute);
    blackEv->setEnabled(true);
    whiteEv->setEnabled(true);
//    targetGray->setEnabled(true);

    sourceGray->setValue(params.sourceGray);
    blackEv->setValue(params.blackEv);
    whiteEv->setValue(params.whiteEv);
//    targetGray->setValue(params.targetGray);
    
    enableListener();
}
