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
    auto EVENT = DISPLAY;
    EvLocalContrastEnabled = m->newEvent(EVENT, "HISTORY_MSG_LOCALCONTRAST_ENABLED");
    EvLocalContrastContrast = m->newEvent(EVENT, "HISTORY_MSG_LOCALCONTRAST_CONTRAST");
    EvLocalContrastCurve = m->newEvent(EVENT, "HISTORY_MSG_LOCALCONTRAST_CURVE");

    Gtk::VBox *wavelets = Gtk::manage(new Gtk::VBox());

    contrast = Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_CONTRAST"), -100., 100., 0.1, 0.));

    const LocalContrastParams default_params;
    
    cg = Gtk::manage(new CurveEditorGroup(options.lastColorToningCurvesDir, M("TP_LOCALCONTRAST_CURVE"), 0.7));
    cg->setCurveListener(this);
    curve = static_cast<FlatCurveEditor *>(cg->addCurve(CT_Flat, "", nullptr, false, false));
    curve->setIdentityValue(0.);
    curve->setResetCurve(FlatCurveType(default_params.curve.at(0)), default_params.curve);
    cg->curveListComplete();
    cg->show();
    
    contrast->setAdjusterListener(this);

    wavelets->pack_start(*cg);
    wavelets->pack_start(*contrast);

    pack_start(*wavelets);
}


void LocalContrast::read(const ProcParams *pp)
{
    disableListener();

    setEnabled(pp->localContrast.enabled);
    contrast->setValue(pp->localContrast.contrast);
    curve->setCurve(pp->localContrast.curve);

    enableListener();
}


void LocalContrast::write(ProcParams *pp)
{
    pp->localContrast.enabled = getEnabled();
    pp->localContrast.contrast = contrast->getValue();
    pp->localContrast.curve = curve->getCurve();
}

void LocalContrast::setDefaults(const ProcParams *defParams)
{
    contrast->setDefault(defParams->localContrast.contrast);
    curve->setResetCurve(FlatCurveType(defParams->localContrast.curve.at(0)), defParams->localContrast.curve);
}

void LocalContrast::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        if (a == contrast) {
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
