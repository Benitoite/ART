/*
 *  This file is part of RawTherapee.
 */
#include "pcvignette.h"
#include "eventmapper.h"

using namespace rtengine;
using namespace rtengine::procparams;

PCVignette::PCVignette () : FoldableToolPanel(this, "pcvignette", M("TP_PCVIGNETTE_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    EvCenter = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_PCVIGNETTE_CENTER");

    strength = Gtk::manage (new Adjuster (M("TP_PCVIGNETTE_STRENGTH"), -6, 6, 0.01, 0));
    strength->set_tooltip_text (M("TP_PCVIGNETTE_STRENGTH_TOOLTIP"));
    strength->setAdjusterListener (this);

    feather = Gtk::manage (new Adjuster (M("TP_PCVIGNETTE_FEATHER"), 0, 100, 1, 50));
    feather->set_tooltip_text (M("TP_PCVIGNETTE_FEATHER_TOOLTIP"));
    feather->setAdjusterListener (this);

    roundness = Gtk::manage (new Adjuster (M("TP_PCVIGNETTE_ROUNDNESS"), 0, 100, 1, 50));
    roundness->set_tooltip_text (M("TP_PCVIGNETTE_ROUNDNESS_TOOLTIP"));
    roundness->setAdjusterListener (this);

    centerX = Gtk::manage(new Adjuster(M("TP_PCVIGNETTE_CENTER_X"), -100, 100, 1, 0));
    centerX->setAdjusterListener(this);
    centerY = Gtk::manage(new Adjuster(M("TP_PCVIGNETTE_CENTER_Y"), -100, 100, 1, 0));
    centerY->setAdjusterListener(this);

    pack_start(*strength);
    pack_start(*feather);
    pack_start(*roundness);
    pack_start(*centerX);
    pack_start(*centerY);

    show_all();
}

void PCVignette::read (const ProcParams* pp)
{
    disableListener ();

    setEnabled(pp->pcvignette.enabled);
    strength->setValue(pp->pcvignette.strength);
    feather->setValue(pp->pcvignette.feather);
    roundness->setValue(pp->pcvignette.roundness);
    centerX->setValue(pp->pcvignette.centerX);
    centerY->setValue(pp->pcvignette.centerY);

    enableListener ();
}

void PCVignette::write(ProcParams* pp)
{
    pp->pcvignette.strength = strength->getValue ();
    pp->pcvignette.feather = feather->getIntValue ();
    pp->pcvignette.roundness = roundness->getIntValue ();
    pp->pcvignette.enabled = getEnabled();
    pp->pcvignette.centerX = centerX->getIntValue();
    pp->pcvignette.centerY = centerY->getIntValue();
}


void PCVignette::setDefaults(const ProcParams* defParams)
{
    strength->setDefault (defParams->pcvignette.strength);
    feather->setDefault (defParams->pcvignette.feather);
    roundness->setDefault (defParams->pcvignette.roundness);
    centerX->setDefault(defParams->pcvignette.centerX);
    centerY->setDefault(defParams->pcvignette.centerY);
}

void PCVignette::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        if (a == strength) {
            listener->panelChanged (EvPCVignetteStrength, strength->getTextValue());
        } else if (a == feather) {
            listener->panelChanged (EvPCVignetteFeather, feather->getTextValue());
        } else if (a == roundness) {
            listener->panelChanged (EvPCVignetteRoundness, roundness->getTextValue());
        } else if (a == centerX || a == centerY) {
            listener->panelChanged(EvCenter, Glib::ustring::compose("X=%1 Y=%2", centerX->getTextValue(), centerY->getTextValue()));
        }
    }
}

void PCVignette::adjusterAutoToggled(Adjuster* a, bool newval)
{
}

void PCVignette::enabledChanged ()
{

    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvPCVignetteEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvPCVignetteEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvPCVignetteEnabled, M("GENERAL_DISABLED"));
        }
    }
}

void PCVignette::trimValues (rtengine::procparams::ProcParams* pp)
{
    strength->trimValue(pp->pcvignette.strength);
    feather->trimValue(pp->pcvignette.feather);
    roundness->trimValue(pp->pcvignette.roundness);
    centerX->trimValue(pp->pcvignette.centerX);
    centerY->trimValue(pp->pcvignette.centerY);
}

