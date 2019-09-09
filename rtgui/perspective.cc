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
#include "eventmapper.h"

using namespace rtengine;
using namespace rtengine::procparams;

PerspCorrection::PerspCorrection() : FoldableToolPanel(this, "perspective", M("TP_PERSPECTIVE_LABEL"), false, true)
{
    EvToolEnabled.set_action(TRANSFORM);
    auto m = ProcEventMapper::getInstance();
    EvPerspCorrLens = m->newEvent(TRANSFORM, "HISTORY_MSG_PERSPECTIVE_LENS");
    
    lgl = nullptr;
    metadata = nullptr;

    Gtk::Image *ipersHL = Gtk::manage(new RTImage("perspective-horizontal-left-small.png"));
    Gtk::Image *ipersHR = Gtk::manage(new RTImage("perspective-horizontal-right-small.png"));
    Gtk::Image *ipersVL = Gtk::manage(new RTImage("perspective-vertical-bottom-small.png"));
    Gtk::Image *ipersVR = Gtk::manage(new RTImage("perspective-vertical-top-small.png"));
    Gtk::Image *ipersSL = Gtk::manage(new RTImage("perspective-shear-left-small.png"));
    Gtk::Image *ipersSR = Gtk::manage(new RTImage("perspective-shear-right-small.png"));
    Gtk::Image *irotateL = Gtk::manage(new RTImage("rotate-right-small.png"));
    Gtk::Image *irotateR = Gtk::manage(new RTImage("rotate-left-small.png"));
    Gtk::Image *iaspectL = Gtk::manage(new RTImage("perspective-aspect-vertical-small.png"));
    Gtk::Image *iaspectR = Gtk::manage(new RTImage("perspective-aspect-horizontal-small.png"));

    horiz = Gtk::manage(new Adjuster(M("TP_PERSPECTIVE_HORIZONTAL"), -200, 200, 0.1, 0, ipersHL, ipersHR));
    horiz->setAdjusterListener(this);

    vert = Gtk::manage(new Adjuster(M("TP_PERSPECTIVE_VERTICAL"), -200, 200, 0.1, 0, ipersVL, ipersVR));
    vert->setAdjusterListener(this);

    angle = Gtk::manage(new Adjuster(M("TP_PERSPECTIVE_ANGLE"), -20, 20, 0.01, 0, irotateL, irotateR));
    shear = Gtk::manage(new Adjuster(M("TP_PERSPECTIVE_SHEAR"), -50, 50, 0.1, 0, ipersSL, ipersSR));
    flength = Gtk::manage(new Adjuster(M("TP_PERSPECTIVE_FLENGTH"), 1, 1200, 0.1, 28));
    cropfactor = Gtk::manage(new Adjuster(M("TP_PERSPECTIVE_CROPFACTOR"), 0.1, 10, 0.01, 1));
    aspect = Gtk::manage(new Adjuster(M("TP_PERSPECTIVE_ASPECT"), 0.5, 2, 0.001, 1, iaspectL, iaspectR));
    angle->setAdjusterListener(this);
    shear->setAdjusterListener(this);
    flength->setAdjusterListener(this);
    cropfactor->setAdjusterListener(this);
    aspect->setAdjusterListener(this);

    pack_start(*horiz);
    pack_start(*vert);
    pack_start(*shear);
    pack_start(*angle);
    pack_start(*aspect);
    pack_start(*flength);
    pack_start(*cropfactor);

    horiz->setLogScale(2, 0);
    vert->setLogScale(2, 0);
    flength->setLogScale(100, 1);
    cropfactor->setLogScale(2, 1);
    aspect->setLogScale(100, 1, true);

    auto_horiz = Gtk::manage(new Gtk::Button());
    auto_horiz->add(*Gtk::manage(new RTImage("perspective-horizontal-left.png")));
    auto_horiz->get_style_context()->add_class("independent");
    auto_horiz->set_tooltip_markup(M("TP_PERSPECTIVE_AUTO_HORIZONTAL_TOOLTIP"));

    auto_vert = Gtk::manage(new Gtk::Button());
    auto_vert->add(*Gtk::manage(new RTImage("perspective-vertical-top.png")));
    auto_vert->get_style_context()->add_class("independent");
    auto_vert->set_tooltip_markup(M("TP_PERSPECTIVE_AUTO_VERTICAL_TOOLTIP"));

    auto_both = Gtk::manage(new Gtk::Button());
    auto_both->add(*Gtk::manage(new RTImage("perspective-horizontal-vertical.png")));
    auto_both->get_style_context()->add_class("independent");
    auto_both->set_tooltip_markup(M("TP_PERSPECTIVE_AUTO_BOTH_TOOLTIP"));

    Gtk::HBox *hb = Gtk::manage(new Gtk::HBox());
    hb->pack_start(*Gtk::manage(new Gtk::Label(M("TP_PERSPECTIVE_AUTO") + ": ")), Gtk::PACK_SHRINK, 4);
    hb->pack_start(*auto_horiz, Gtk::PACK_EXPAND_WIDGET, 2);
    hb->pack_start(*auto_vert, Gtk::PACK_EXPAND_WIDGET, 2);
    hb->pack_start(*auto_both, Gtk::PACK_EXPAND_WIDGET, 2);
    auto_horiz->show();
    auto_vert->show();
    auto_both->show();

    auto_horiz->signal_pressed().connect(sigc::bind(sigc::mem_fun(*this, &PerspCorrection::autoPressed), auto_horiz));
    auto_vert->signal_pressed().connect(sigc::bind(sigc::mem_fun(*this, &PerspCorrection::autoPressed), auto_vert));
    auto_both->signal_pressed().connect(sigc::bind(sigc::mem_fun(*this, &PerspCorrection::autoPressed), auto_both));

    pack_start(*hb, Gtk::PACK_EXPAND_WIDGET, 4);
    hb->show();
    
    show_all();
}

void PerspCorrection::read(const ProcParams* pp)
{
    disableListener ();

    horiz->setValue (pp->perspective.horizontal);
    vert->setValue (pp->perspective.vertical);
    angle->setValue(pp->perspective.angle);
    shear->setValue(pp->perspective.shear);
    aspect->setValue(pp->perspective.aspect);
    if (pp->perspective.flength > 0) {
        flength->setValue(pp->perspective.flength);
        cropfactor->setValue(pp->perspective.cropfactor);
    } else if (metadata) {
        do_set_metadata(metadata);
    } else {
        flength->setValue(28);
        cropfactor->setValue(1);
    }
    setEnabled(pp->perspective.enabled);

    enableListener ();
}

void PerspCorrection::write(ProcParams* pp)
{
    pp->perspective.enabled = getEnabled();
    pp->perspective.horizontal  = horiz->getValue ();
    pp->perspective.vertical = vert->getValue ();
    pp->perspective.angle = angle->getValue();
    pp->perspective.shear = shear->getValue();
    pp->perspective.flength = flength->getValue();
    pp->perspective.cropfactor = cropfactor->getValue();
    pp->perspective.aspect = aspect->getValue();
}

void PerspCorrection::setDefaults(const ProcParams* defParams)
{

    horiz->setDefault (defParams->perspective.horizontal);
    vert->setDefault (defParams->perspective.vertical);
    angle->setDefault(defParams->perspective.angle);
    shear->setDefault(defParams->perspective.shear);
    flength->setDefault(defParams->perspective.flength);
    cropfactor->setDefault(defParams->perspective.cropfactor);
    aspect->setDefault(defParams->perspective.aspect);
}

void PerspCorrection::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        if (a == flength || a == cropfactor || a == aspect) {
            listener->panelChanged(EvPerspCorrLens, Glib::ustring::compose("%1=%2  %3=%4\n%5=%6", M("TP_PERSPECTIVE_FLENGTH"), flength->getValue(), M("TP_PERSPECTIVE_CROPFACTOR"), cropfactor->getValue(), M("TP_PERSPECTIVE_ASPECT"), aspect->getValue()));
        } else {
            listener->panelChanged(EvPerspCorr, Glib::ustring::compose ("%1=%2  %3=%4\n%5=%6  %7=%8", M("TP_PERSPECTIVE_HORIZONTAL"), horiz->getValue(), M("TP_PERSPECTIVE_VERTICAL"), vert->getValue(), M("TP_PERSPECTIVE_ANGLE"), angle->getValue(), M("TP_PERSPECTIVE_SHEAR"), shear->getValue()));
        }
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
    flength->trimValue(pp->perspective.flength);
    cropfactor->trimValue(pp->perspective.cropfactor);
    aspect->trimValue(pp->perspective.aspect);
}


void PerspCorrection::autoPressed(Gtk::Button *which)
{
    if (!lgl || !getEnabled()) {
        return;
    }
    
    double a, h, v, s;
    bool hh = false;
    bool vv = false;

    if (which == auto_horiz) {
        hh = true;
    } else if (which == auto_vert) {
        vv = true;
    } else if (which == auto_both) {
        hh = true;
        vv = true;
    }

    lgl->autoPerspectiveRequested(hh, vv, a, h, v, s);

    disableListener();
    angle->setValue(a);
    horiz->setValue(h);
    vert->setValue(v);
    shear->setValue(s);
    enableListener();

    adjusterChanged(nullptr, 0);
}


void PerspCorrection::do_set_metadata(const rtengine::FramesMetaData *meta)
{
    metadata = meta;
    if (metadata) {
        double f = metadata->getFocalLen();
        double f35 = metadata->getFocalLen35mm();
        if (f > 0) {
            flength->setValue(f);
            flength->setDefault(f, true);
            if (f35 > 0) {
                double crop = f35 / f;
                cropfactor->setValue(crop);
                cropfactor->setDefault(crop, true);
            }
        }
    }
}


void PerspCorrection::setRawMeta(bool raw, const rtengine::FramesMetaData *meta)
{
    disableListener();
    do_set_metadata(meta);
    enableListener();
}
