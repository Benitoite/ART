/*
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
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
#include "perspective.h"
#include "rtimage.h"

using namespace rtengine;
using namespace rtengine::procparams;

PerspCorrection::PerspCorrection () : FoldableToolPanel(this, "perspective", M("TP_PERSPECTIVE_LABEL"))
{

    Gtk::Image* ipersHL =   Gtk::manage (new RTImage ("perspective-horizontal-left-small.png"));
    Gtk::Image* ipersHR =   Gtk::manage (new RTImage ("perspective-horizontal-right-small.png"));
    Gtk::Image* ipersVL =   Gtk::manage (new RTImage ("perspective-vertical-bottom-small.png"));
    Gtk::Image* ipersVR =   Gtk::manage (new RTImage ("perspective-vertical-top-small.png"));

    horiz = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_HORIZONTAL"), -100, 100, 0.1, 0, ipersHR, ipersHL));
    horiz->setAdjusterListener (this);

    vert = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_VERTICAL"), -100, 100, 0.1, 0, ipersVL, ipersVR));
    vert->setAdjusterListener (this);

    angle = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_ANGLE"), -20, 20, 0.01, 0));
    shear = Gtk::manage (new Adjuster (M("TP_PERSPECTIVE_SHEAR"), -50, 50, 0.1, 0));
    angle->setAdjusterListener(this);
    shear->setAdjusterListener(this);

    pack_start (*horiz);
    pack_start (*vert);
    pack_start(*shear);
    pack_start(*angle);

    horiz->setLogScale(2, 0);
    vert->setLogScale(2, 0);

    show_all();
}

void PerspCorrection::read(const ProcParams* pp)
{

    disableListener ();

    horiz->setValue (pp->perspective.horizontal);
    vert->setValue (pp->perspective.vertical);
    angle->setValue(pp->perspective.angle);
    shear->setValue(pp->perspective.shear);

    enableListener ();
}

void PerspCorrection::write(ProcParams* pp)
{

    pp->perspective.horizontal  = horiz->getValue ();
    pp->perspective.vertical = vert->getValue ();
    pp->perspective.angle = angle->getValue();
    pp->perspective.shear = shear->getValue();
}

void PerspCorrection::setDefaults(const ProcParams* defParams)
{

    horiz->setDefault (defParams->perspective.horizontal);
    vert->setDefault (defParams->perspective.vertical);
    angle->setDefault(defParams->perspective.angle);
    shear->setDefault(defParams->perspective.shear);
}

void PerspCorrection::adjusterChanged(Adjuster* a, double newval)
{
    if (listener) {
        listener->panelChanged (EvPerspCorr, Glib::ustring::compose ("%1=%2  %3=%4\n%5=%6  %7=%8", M("TP_PERSPECTIVE_HORIZONTAL"), horiz->getValue(), M("TP_PERSPECTIVE_VERTICAL"), vert->getValue(), M("TP_PERSPECTIVE_ANGLE"), angle->getValue(), M("TP_PERSPECTIVE_SHEAR"), shear->getValue()));
    }
}

void PerspCorrection::adjusterAutoToggled(Adjuster* a, bool newval)
{
}


void PerspCorrection::trimValues (rtengine::procparams::ProcParams* pp)
{

    horiz->trimValue(pp->perspective.horizontal);
    vert->trimValue(pp->perspective.vertical);
    angle->trimValue(pp->perspective.angle);
    shear->trimValue(pp->perspective.shear);
}

