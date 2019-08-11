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
#include <cmath>
#include <sstream>
#include "sharpenmicro.h"
#include "guiutils.h"
#include "eventmapper.h"

using namespace rtengine;
using namespace rtengine::procparams;


SharpenMicro::SharpenMicro () : FoldableToolPanel(this, "sharpenmicro", M("TP_SHARPENMICRO_LABEL"), true, true)
{

    auto m = ProcEventMapper::getInstance();
    EvSharpenMicroContrast = m->newEvent(DIRPYREQUALIZER, "HISTORY_MSG_MICROCONTRAST_CONTRAST");

    contrast = Gtk::manage(new Adjuster (M("TP_SHARPENMICRO_CONTRAST"), 0, 200, 1, 20));
    contrast->setAdjusterListener (this);
    contrast->show();

    amount = Gtk::manage(new Adjuster (M("TP_SHARPENMICRO_AMOUNT"), 0, 100, 1, 20));
    amount->setAdjusterListener (this);
    amount->show();

    uniformity = Gtk::manage(new Adjuster (M("TP_SHARPENMICRO_UNIFORMITY"), 0, 10, 1, 5));
    uniformity->setAdjusterListener (this);
    uniformity->show();

    pack_start( *contrast, Gtk::PACK_SHRINK, 0);
    pack_start( *amount, Gtk::PACK_SHRINK, 0);
    pack_start( *uniformity, Gtk::PACK_SHRINK, 0);

    matrix = Gtk::manage (new Gtk::CheckButton (M("TP_SHARPENMICRO_MATRIX")));
    matrix->set_active (true);
    pack_start(*matrix, Gtk::PACK_SHRINK, 0);
    matrix->show ();

    matrixconn = matrix->signal_toggled().connect( sigc::mem_fun(*this, &SharpenMicro::matrix_toggled) );
}

void SharpenMicro::read(const ProcParams* pp)
{
    disableListener ();

    setEnabled(pp->sharpenMicro.enabled);

    matrixconn.block (true);
    matrix->set_active (pp->sharpenMicro.matrix);
    matrixconn.block (false);
    lastmatrix = pp->sharpenMicro.matrix;

    amount->setValue     (pp->sharpenMicro.amount);
    contrast->setValue   (pp->sharpenMicro.contrast);
    uniformity->setValue (pp->sharpenMicro.uniformity);

    enableListener ();
}

void SharpenMicro::write(ProcParams* pp)
{
    pp->sharpenMicro.enabled    = getEnabled();
    pp->sharpenMicro.matrix     = matrix->get_active ();
    pp->sharpenMicro.amount     = amount->getValue ();
    pp->sharpenMicro.contrast   = contrast->getValue ();
    pp->sharpenMicro.uniformity = uniformity->getValue ();
}

void SharpenMicro::enabledChanged ()
{

    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvSharpenMicroEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvSharpenMicroEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvSharpenMicroEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void SharpenMicro::matrix_toggled ()
{
    if (listener && getEnabled()) {
        if (matrix->get_active ()) {
            listener->panelChanged (EvSharpenMicroMatrix, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvSharpenMicroMatrix, M("GENERAL_DISABLED"));
        }
    }
}

void SharpenMicro::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        const Glib::ustring value = a->getTextValue();

        if (a == amount) {
            listener->panelChanged (EvSharpenMicroAmount,     value);
        } else if (a == contrast) {
            listener->panelChanged (EvSharpenMicroContrast,   value);
        } else if (a == uniformity) {
            listener->panelChanged (EvSharpenMicroUniformity, value);
        }
    }
}

void SharpenMicro::adjusterAutoToggled(Adjuster* a, bool newval)
{
}

void SharpenMicro::setDefaults(const ProcParams* defParams)
{
    amount->setDefault     (defParams->sharpenMicro.amount);
    contrast->setDefault   (defParams->sharpenMicro.contrast);
    uniformity->setDefault (defParams->sharpenMicro.uniformity);
}

void SharpenMicro::trimValues (ProcParams* pp)
{
    amount->trimValue     (pp->sharpenMicro.amount);
    contrast->trimValue   (pp->sharpenMicro.contrast);
    uniformity->trimValue (pp->sharpenMicro.uniformity);
}
