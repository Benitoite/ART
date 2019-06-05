/*
 *  This file is part of RawTherapee.
 */
#include "pcvignette.h"

using namespace rtengine;
using namespace rtengine::procparams;

PCVignette::PCVignette () : FoldableToolPanel(this, "pcvignette", M("TP_PCVIGNETTE_LABEL"), false, true)
{
    strength = Gtk::manage (new Adjuster (M("TP_PCVIGNETTE_STRENGTH"), -6, 6, 0.01, 0));
    strength->set_tooltip_text (M("TP_PCVIGNETTE_STRENGTH_TOOLTIP"));
    strength->setAdjusterListener (this);

    feather = Gtk::manage (new Adjuster (M("TP_PCVIGNETTE_FEATHER"), 0, 100, 1, 50));
    feather->set_tooltip_text (M("TP_PCVIGNETTE_FEATHER_TOOLTIP"));
    feather->setAdjusterListener (this);

    roundness = Gtk::manage (new Adjuster (M("TP_PCVIGNETTE_ROUNDNESS"), 0, 100, 1, 50));
    roundness->set_tooltip_text (M("TP_PCVIGNETTE_ROUNDNESS_TOOLTIP"));
    roundness->setAdjusterListener (this);

    pack_start (*strength);
    pack_start (*feather);
    pack_start (*roundness);

    show_all();
}

void PCVignette::read (const ProcParams* pp)
{
    disableListener ();

    setEnabled(pp->pcvignette.enabled);
    strength->setValue (pp->pcvignette.strength);
    feather->setValue (pp->pcvignette.feather);
    roundness->setValue (pp->pcvignette.roundness);

    enableListener ();
}

void PCVignette::write(ProcParams* pp)
{
    pp->pcvignette.strength = strength->getValue ();
    pp->pcvignette.feather = feather->getIntValue ();
    pp->pcvignette.roundness = roundness->getIntValue ();
    pp->pcvignette.enabled = getEnabled();
}

void PCVignette::setDefaults(const ProcParams* defParams)
{
    strength->setDefault (defParams->pcvignette.strength);
    feather->setDefault (defParams->pcvignette.feather);
    roundness->setDefault (defParams->pcvignette.roundness);
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
}

