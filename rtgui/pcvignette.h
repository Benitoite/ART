/*
 *  This file is part of RawTherapee.
 */
#ifndef _PCVIGNETTE_H_
#define _PCVIGNETTE_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"

class PCVignette : public ToolParamBlock, public AdjusterListener, public FoldableToolPanel
{

protected:
    Adjuster* strength;
    Adjuster* feather;
    Adjuster* roundness;

public:

    PCVignette ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setDefaults    (const rtengine::procparams::ProcParams* defParams, const ParamsEdited* pedited = nullptr);
    void setBatchMode   (bool batchMode);

    void adjusterChanged (Adjuster* a, double newval);
    void adjusterAutoToggled(Adjuster* a, bool newval);
    void enabledChanged  ();
    void setAdjusterBehavior (bool strengthadd, bool featheradd, bool roundnessadd);
    void trimValues          (rtengine::procparams::ProcParams* pp);
};

#endif
