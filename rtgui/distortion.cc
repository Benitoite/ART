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
#include "distortion.h"
#include <iomanip>
#include "rtimage.h"
#include "eventmapper.h"

using namespace rtengine;
using namespace rtengine::procparams;

Distortion::Distortion (): FoldableToolPanel(this, "distortion", M("TP_DISTORTION_LABEL"), false, true)
{
    rlistener = nullptr;

    EvToolEnabled.set_action(TRANSFORM);
    
    autoDistor = Gtk::manage (new Gtk::Button (M("GENERAL_AUTO")));
    autoDistor->set_image (*Gtk::manage (new RTImage ("distortion-auto-small.png")));
    autoDistor->get_style_context()->add_class("independent");
    autoDistor->set_alignment(0.5f, 0.5f);
    autoDistor->set_tooltip_text (M("TP_DISTORTION_AUTO_TIP"));
    idConn = autoDistor->signal_pressed().connect( sigc::mem_fun(*this, &Distortion::idPressed) );
    autoDistor->show();
    pack_start (*autoDistor);

    Gtk::Image* idistL =   Gtk::manage (new RTImage ("distortion-pincushion-small.png"));
    Gtk::Image* idistR =   Gtk::manage (new RTImage ("distortion-barrel-small.png"));

    distor = Gtk::manage (new Adjuster (M("TP_DISTORTION_AMOUNT"), -0.5, 0.5, 0.001, 0, idistL, idistR));
    distor->setAdjusterListener (this);

    distor->setLogScale(2, 0);
    
    distor->show();
    pack_start (*distor);
}

void Distortion::read(const ProcParams* pp)
{
    disableListener ();
    setEnabled(pp->distortion.enabled);
    distor->setValue (pp->distortion.amount);
    enableListener ();
}

void Distortion::write(ProcParams* pp)
{
    pp->distortion.enabled = getEnabled();
    pp->distortion.amount = distor->getValue ();
}

void Distortion::setDefaults(const ProcParams* defParams)
{
    distor->setDefault (defParams->distortion.amount);
}

void Distortion::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        listener->panelChanged (EvDISTAmount, Glib::ustring::format (std::setw(4), std::fixed, std::setprecision(3), a->getValue()));
    }
}

void Distortion::adjusterAutoToggled(Adjuster* a, bool newval)
{
}

void Distortion::idPressed ()
{
    setEnabled(true);
    if (rlistener) {
        double new_amount = rlistener->autoDistorRequested();
        distor->setValue(new_amount);
        adjusterChanged (distor, new_amount);
    }
}

void Distortion::trimValues (rtengine::procparams::ProcParams* pp)
{

    distor->trimValue(pp->distortion.amount);
}
