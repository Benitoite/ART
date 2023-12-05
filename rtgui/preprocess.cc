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
#include "preprocess.h"
#include "guiutils.h"
#include <sstream>
#include "../rtengine/refreshmap.h"

using namespace rtengine;
using namespace rtengine::procparams;

PreProcess::PreProcess () : FoldableToolPanel(this, "preprocess", M("TP_PREPROCESS_LABEL"), false, true, true)
{
    EvToolEnabled.set_action(DARKFRAME);
    EvToolReset.set_action(DARKFRAME);

    Gtk::HBox* hotdeadPixel = Gtk::manage( new Gtk::HBox () );
    hotdeadPixel->set_spacing(4);
    hotPixel = Gtk::manage(new Gtk::CheckButton((M("TP_PREPROCESS_HOTPIXFILT"))));
    deadPixel = Gtk::manage(new Gtk::CheckButton((M("TP_PREPROCESS_DEADPIXFILT"))));

    hotPixel->set_tooltip_markup (M("TP_PREPROCESS_HOTPIXFILT_TOOLTIP"));
    deadPixel->set_tooltip_markup (M("TP_PREPROCESS_DEADPIXFILT_TOOLTIP"));

    hotdeadPixel->pack_start( *hotPixel, Gtk::PACK_SHRINK);
    hotdeadPixel->pack_start( *deadPixel, Gtk::PACK_SHRINK, 0);
    pack_start(*hotdeadPixel, Gtk::PACK_SHRINK, 0);
    hdThreshold = Gtk::manage (new Adjuster (M("TP_RAW_HD"), 20, 200, 2, 100));
    hdThreshold->set_tooltip_markup (M("TP_RAW_HD_TOOLTIP"));
    hdThreshold->setAdjusterListener (this);

    if (hdThreshold->delay < options.adjusterMaxDelay) {
        hdThreshold->delay = options.adjusterMaxDelay;
    }

    hdThreshold->show();
    pack_start( *hdThreshold, Gtk::PACK_SHRINK, 4);

//  hotdeadPixel->show();
    hpixelconn = hotPixel->signal_toggled().connect ( sigc::mem_fun(*this, &PreProcess::hotPixelChanged), true);
    dpixelconn = deadPixel->signal_toggled().connect ( sigc::mem_fun(*this, &PreProcess::deadPixelChanged), true);
}

void PreProcess::read(const rtengine::procparams::ProcParams* pp)
{
    disableListener ();
    hpixelconn.block (true);
    dpixelconn.block (true);

    setEnabled(pp->raw.enable_hotdeadpix);

    lastHot = pp->raw.hotPixelFilter;
    lastDead = pp->raw.deadPixelFilter;
    hotPixel->set_active (pp->raw.hotPixelFilter);
    deadPixel->set_active (pp->raw.deadPixelFilter);
    hdThreshold->setValue (pp->raw.hotdeadpix_thresh);
    hpixelconn.block (false);
    dpixelconn.block (false);
    enableListener ();
}

void PreProcess::write( rtengine::procparams::ProcParams* pp)
{
    pp->raw.enable_hotdeadpix = getEnabled();
    pp->raw.hotPixelFilter = hotPixel->get_active();
    pp->raw.deadPixelFilter = deadPixel->get_active();
    pp->raw.hotdeadpix_thresh = hdThreshold->getIntValue();
}

void PreProcess::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        if (a == hdThreshold) {
            listener->panelChanged (EvPreProcessHotDeadThresh, a->getTextValue() );
        }
    }
}

void PreProcess::adjusterAutoToggled(Adjuster* a, bool newval)
{
}

void PreProcess::hotPixelChanged ()
{
    if (listener && getEnabled()) {
        listener->panelChanged (EvPreProcessHotPixel, hotPixel->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}

void PreProcess::deadPixelChanged ()
{
    if (listener && getEnabled()) {
        listener->panelChanged (EvPreProcessDeadPixel, deadPixel->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}


void PreProcess::setDefaults(const ProcParams *def)
{
    hdThreshold->setDefault(def->raw.hotdeadpix_thresh);
    initial_params = def->raw;
}


void PreProcess::toolReset(bool to_initial)
{
    ProcParams pp;
    if (to_initial) {
        pp.raw = initial_params;
    }
    pp.raw.enable_hotdeadpix = getEnabled();
    read(&pp);
}
