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
#include "localcontrast.h"
#include "eventmapper.h"
#include <iomanip>
#include <cmath>

using namespace rtengine;
using namespace rtengine::procparams;

LocalContrast::LocalContrast(): FoldableToolPanel(this, "localcontrast", M("TP_LOCALCONTRAST_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    EvLocalContrastEnabled = m->newEvent(DIRPYREQUALIZER, "HISTORY_MSG_LOCALCONTRAST_ENABLED");
    EvLocalContrastMode = m->newEvent(DIRPYREQUALIZER, "HISTORY_MSG_LOCALCONTRAST_MODE");
    EvLocalContrastRadius = m->newEvent(DIRPYREQUALIZER, "HISTORY_MSG_LOCALCONTRAST_RADIUS");
    EvLocalContrastAmount = m->newEvent(DIRPYREQUALIZER, "HISTORY_MSG_LOCALCONTRAST_AMOUNT");
    EvLocalContrastDarkness = m->newEvent(DIRPYREQUALIZER, "HISTORY_MSG_LOCALCONTRAST_DARKNESS");
    EvLocalContrastLightness = m->newEvent(DIRPYREQUALIZER, "HISTORY_MSG_LOCALCONTRAST_LIGHTNESS");
    EvLocalContrastContrast = m->newEvent(DIRPYREQUALIZER, "HISTORY_MSG_LOCALCONTRAST_CONTRAST");
    EvLocalContrastCurve = m->newEvent(DIRPYREQUALIZER, "HISTORY_MSG_LOCALCONTRAST_CURVE");

    usm = Gtk::manage(new Gtk::VBox());
    wavelets = Gtk::manage(new Gtk::VBox());

    Gtk::HBox *hb = Gtk::manage(new Gtk::HBox());
    mode = Gtk::manage(new MyComboBoxText ());
    mode->append(M("TP_LOCALCONTRAST_USM"));
    mode->append(M("TP_LOCALCONTRAST_WAVELETS"));
    hb->pack_start(*Gtk::manage(new Gtk::Label(M("TP_LOCALCONTRAST_MODE") + ":")), Gtk::PACK_SHRINK, 4);
    hb->pack_start(*mode);
    pack_start(*hb);
    
    radius = Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_RADIUS"), 20., 200., 1., 80.));
    amount = Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_AMOUNT"), 0., 1., 0.01, 0.2));
    darkness = Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_DARKNESS"), 0., 3., 0.01, 1.));
    lightness = Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_LIGHTNESS"), 0., 3., 0.01, 1.));
    contrast = Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_CONTRAST"), -100., 100., 0.1, 0.));

    const LocalContrastParams default_params;
    
    cg = Gtk::manage(new CurveEditorGroup(options.lastColorToningCurvesDir, M("TP_LOCALCONTRAST_CURVE"), 0.7));
    cg->setCurveListener(this);
    curve = static_cast<FlatCurveEditor *>(cg->addCurve(CT_Flat, "", nullptr, false, false));
    curve->setIdentityValue(0.);
    curve->setResetCurve(FlatCurveType(default_params.curve.at(0)), default_params.curve);
    cg->curveListComplete();
    cg->show();
    
    mode->signal_changed().connect(sigc::mem_fun(*this, &LocalContrast::modeChanged));
    
    radius->setAdjusterListener(this);
    amount->setAdjusterListener(this);
    darkness->setAdjusterListener(this);
    lightness->setAdjusterListener(this);
    contrast->setAdjusterListener(this);

    radius->show();
    amount->show();
    darkness->show();
    lightness->show();

    usm->pack_start(*radius);
    usm->pack_start(*amount);
    usm->pack_start(*darkness);
    usm->pack_start(*lightness);

    wavelets->pack_start(*cg);
    wavelets->pack_start(*contrast);

    pack_start(*usm);
    pack_start(*wavelets);
}


void LocalContrast::read(const ProcParams *pp)
{
    disableListener();

    setEnabled(pp->localContrast.enabled);
    if (pp->localContrast.mode == LocalContrastParams::USM) {
        mode->set_active(0);
    } else {
        mode->set_active(1);
    }
    radius->setValue(pp->localContrast.radius);
    amount->setValue(pp->localContrast.amount);
    darkness->setValue(pp->localContrast.darkness);
    lightness->setValue(pp->localContrast.lightness);
    contrast->setValue(pp->localContrast.contrast);
    curve->setCurve(pp->localContrast.curve);

    modeChanged();

    enableListener();
}


void LocalContrast::write(ProcParams *pp)
{
    pp->localContrast.mode = LocalContrastParams::Mode(min(mode->get_active_row_number(), 1));
    pp->localContrast.radius = radius->getValue();
    pp->localContrast.amount = amount->getValue();
    pp->localContrast.darkness = darkness->getValue();
    pp->localContrast.lightness = lightness->getValue();
    pp->localContrast.enabled = getEnabled();
    pp->localContrast.contrast = contrast->getValue();
    pp->localContrast.curve = curve->getCurve();
}

void LocalContrast::setDefaults(const ProcParams *defParams)
{
    radius->setDefault(defParams->localContrast.radius);
    amount->setDefault(defParams->localContrast.amount);
    darkness->setDefault(defParams->localContrast.darkness);
    lightness->setDefault(defParams->localContrast.lightness);
    contrast->setDefault(defParams->localContrast.contrast);
    curve->setResetCurve(FlatCurveType(defParams->localContrast.curve.at(0)), defParams->localContrast.curve);
}

void LocalContrast::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        if (a == radius) {
            listener->panelChanged(EvLocalContrastRadius, a->getTextValue());
        } else if (a == amount) {
            listener->panelChanged(EvLocalContrastAmount, a->getTextValue());
        } else if (a == darkness) {
            listener->panelChanged(EvLocalContrastDarkness, a->getTextValue());
        } else if (a == lightness) {
            listener->panelChanged(EvLocalContrastLightness, a->getTextValue());
        } else if (a == contrast) {
            listener->panelChanged(EvLocalContrastContrast, a->getTextValue());
        }
    }
}

void LocalContrast::adjusterAutoToggled(Adjuster* a, bool newval)
{
}

void LocalContrast::enabledChanged ()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(EvLocalContrastEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(EvLocalContrastEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvLocalContrastEnabled, M("GENERAL_DISABLED"));
        }
    }
}


void LocalContrast::modeChanged()
{
    removeIfThere(this, usm, false);
    removeIfThere(this, wavelets, false);

    if (mode->get_active_row_number() == 0) {
        pack_start(*usm);
    } else if (mode->get_active_row_number() == 1) {
        pack_start(*wavelets);
    }

    if (listener && getEnabled() ) {
        listener->panelChanged(EvLocalContrastMode, mode->get_active_text());
    }
}


void LocalContrast::curveChanged()
{
    if (listener) {
        listener->panelChanged(EvLocalContrastCurve, M("HISTORY_CUSTOMCURVE"));
    }
}


void LocalContrast::autoOpenCurve()
{
    curve->openIfNonlinear();
}
