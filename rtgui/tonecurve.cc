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
#include "tonecurve.h"
#include "adjuster.h"
#include <sigc++/slot.h>
#include <iomanip>
#include "ppversion.h"
#include "edit.h"
#include "eventmapper.h"

using namespace rtengine;
using namespace rtengine::procparams;

namespace {

const std::vector<double> default_satcurve{
    FCT_MinMaxCPoints,
    0, 0.397626,
    0.843827, 0,
    0.5, 0.5,
    0.35, 0.35,
    1, 0.231454,
    0, 0.35
};

} // namespace


ToneCurve::ToneCurve():
    FoldableToolPanel(this, "tonecurve", M("TP_TONECURVE_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    EvHistMatching = m->newEvent(AUTOEXP, "HISTORY_MSG_HISTMATCHING");
    EvHistMatchingBatch = m->newEvent(M_VOID, "HISTORY_MSG_HISTMATCHING");
    EvSatCurve = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_TONECURVE_SATCURVE");
    EvToolEnabled.set_action(AUTOEXP);

    CurveListener::setMulti(true);

    std::vector<GradientMilestone> bottomMilestones;
    bottomMilestones.push_back( GradientMilestone(0., 0., 0., 0.) );
    bottomMilestones.push_back( GradientMilestone(1., 1., 1., 1.) );

//----------- Curve 1 ------------------------------
    pack_start (*Gtk::manage (new  Gtk::HSeparator()));

    histmatching = Gtk::manage(new Gtk::ToggleButton(M("TP_EXPOSURE_HISTMATCHING")));
    histmatching->set_tooltip_markup(M("TP_EXPOSURE_HISTMATCHING_TOOLTIP"));
    histmatchconn = histmatching->signal_toggled().connect(sigc::mem_fun(*this, &ToneCurve::histmatchingToggled));
    pack_start(*histmatching, true, true, 2);

    toneCurveMode = Gtk::manage (new MyComboBoxText ());
    toneCurveMode->append (M("TP_EXPOSURE_TCMODE_STANDARD"));
    toneCurveMode->append (M("TP_EXPOSURE_TCMODE_WEIGHTEDSTD"));
    toneCurveMode->append (M("TP_EXPOSURE_TCMODE_FILMLIKE"));
    toneCurveMode->append (M("TP_EXPOSURE_TCMODE_SATANDVALBLENDING"));
    toneCurveMode->append (M("TP_EXPOSURE_TCMODE_LUMINANCE"));
    toneCurveMode->append (M("TP_EXPOSURE_TCMODE_PERCEPTUAL"));
    toneCurveMode->set_active (0);
    toneCurveMode->set_tooltip_text(M("TP_EXPOSURE_TCMODE_LABEL1"));

    curveEditorG = new CurveEditorGroup (options.lastToneCurvesDir, M("TP_EXPOSURE_CURVEEDITOR1"));
    curveEditorG->setCurveListener (this);

    shape = static_cast<DiagonalCurveEditor*>(curveEditorG->addCurve(CT_Diagonal, "", toneCurveMode));
    shape->setEditID(EUID_ToneCurve1, BT_IMAGEFLOAT);
    shape->setBottomBarBgGradient(bottomMilestones);
    shape->setLeftBarBgGradient(bottomMilestones);

    // This will add the reset button at the end of the curveType buttons
    curveEditorG->curveListComplete();

    pack_start( *curveEditorG, Gtk::PACK_SHRINK, 2);

    tcmodeconn = toneCurveMode->signal_changed().connect( sigc::mem_fun(*this, &ToneCurve::curveMode1Changed), true );

//----------- Curve 2 ------------------------------

    toneCurveMode2 = Gtk::manage (new MyComboBoxText ());
    toneCurveMode2->append (M("TP_EXPOSURE_TCMODE_STANDARD"));
    toneCurveMode2->append (M("TP_EXPOSURE_TCMODE_WEIGHTEDSTD"));
    toneCurveMode2->append (M("TP_EXPOSURE_TCMODE_FILMLIKE"));
    toneCurveMode2->append (M("TP_EXPOSURE_TCMODE_SATANDVALBLENDING"));
    toneCurveMode2->append (M("TP_EXPOSURE_TCMODE_LUMINANCE"));
    toneCurveMode2->append (M("TP_EXPOSURE_TCMODE_PERCEPTUAL"));
    toneCurveMode2->set_active (0);
    toneCurveMode2->set_tooltip_text(M("TP_EXPOSURE_TCMODE_LABEL2"));

    curveEditorG2 = new CurveEditorGroup (options.lastToneCurvesDir, M("TP_EXPOSURE_CURVEEDITOR2"));
    curveEditorG2->setCurveListener (this);

    shape2 = static_cast<DiagonalCurveEditor*>(curveEditorG2->addCurve(CT_Diagonal, "", toneCurveMode2));
    shape2->setEditID(EUID_ToneCurve2, BT_IMAGEFLOAT);
    shape2->setBottomBarBgGradient(bottomMilestones);
    shape2->setLeftBarBgGradient(bottomMilestones);

    // This will add the reset button at the end of the curveType buttons
    curveEditorG2->curveListComplete();
    curveEditorG2->setTooltip(M("TP_EXPOSURE_CURVEEDITOR2_TOOLTIP"));

    pack_start( *curveEditorG2, Gtk::PACK_SHRINK, 2);

    tcmode2conn = toneCurveMode2->signal_changed().connect( sigc::mem_fun(*this, &ToneCurve::curveMode2Changed), true );


    satcurveG = Gtk::manage(new CurveEditorGroup(options.lastColorToningCurvesDir, M("TP_TONECURVE_SATCURVE"), 0.7));
    satcurveG->setCurveListener(this);
    satcurve = static_cast<FlatCurveEditor *>(satcurveG->addCurve(CT_Flat, "", nullptr, false, false));
    satcurve->setResetCurve(FlatCurveType(default_satcurve[0]), default_satcurve);
    satcurve->setBottomBarBgGradient(bottomMilestones);
    satcurveG->curveListComplete();
    satcurveG->show();

    pack_start(*satcurveG, Gtk::PACK_SHRINK, 2);
}


ToneCurve::~ToneCurve ()
{
    idle_register.destroy();

    delete curveEditorG;
    delete curveEditorG2;
}


void ToneCurve::read(const ProcParams* pp)
{
    disableListener ();

    tcmodeconn.block(true);
    tcmode2conn.block(true);

    setEnabled(pp->toneCurve.enabled);

    shape->setCurve (pp->toneCurve.curve);
    shape2->setCurve (pp->toneCurve.curve2);

    toneCurveMode->set_active(rtengine::toUnderlying(pp->toneCurve.curveMode));
    toneCurveMode2->set_active(rtengine::toUnderlying(pp->toneCurve.curveMode2));

    histmatching->set_active(pp->toneCurve.histmatching);
    fromHistMatching = pp->toneCurve.fromHistMatching;

    satcurve->setCurve(default_satcurve);
    satcurve->setCurve(pp->toneCurve.saturation);

    tcmode2conn.block(false);
    tcmodeconn.block(false);

    enableListener ();
}


void ToneCurve::autoOpenCurve  ()
{
    shape->openIfNonlinear();
    shape2->openIfNonlinear();
    satcurve->openIfNonlinear();
}


void ToneCurve::setEditProvider(EditDataProvider *provider)
{
    shape->setEditProvider(provider);
    shape2->setEditProvider(provider);
    //satcurve->setEditProvider(provider);
}


void ToneCurve::write(ProcParams* pp)
{
    pp->toneCurve.enabled = getEnabled();
    pp->toneCurve.curve = shape->getCurve ();
    pp->toneCurve.curve2 = shape2->getCurve ();

    int tcMode = toneCurveMode->get_active_row_number();

    if (tcMode == 0) {
        pp->toneCurve.curveMode = ToneCurveParams::TcMode::STD;
    } else if (tcMode == 1) {
        pp->toneCurve.curveMode = ToneCurveParams::TcMode::WEIGHTEDSTD;
    } else if (tcMode == 2) {
        pp->toneCurve.curveMode = ToneCurveParams::TcMode::FILMLIKE;
    } else if (tcMode == 3) {
        pp->toneCurve.curveMode = ToneCurveParams::TcMode::SATANDVALBLENDING;
    } else if (tcMode == 4) {
        pp->toneCurve.curveMode = ToneCurveParams::TcMode::LUMINANCE;
    } else if (tcMode == 5) {
        pp->toneCurve.curveMode = ToneCurveParams::TcMode::PERCEPTUAL;
    }

    tcMode = toneCurveMode2->get_active_row_number();

    if (tcMode == 0) {
        pp->toneCurve.curveMode2 = ToneCurveParams::TcMode::STD;
    } else if (tcMode == 1) {
        pp->toneCurve.curveMode2 = ToneCurveParams::TcMode::WEIGHTEDSTD;
    } else if (tcMode == 2) {
        pp->toneCurve.curveMode2 = ToneCurveParams::TcMode::FILMLIKE;
    } else if (tcMode == 3) {
        pp->toneCurve.curveMode2 = ToneCurveParams::TcMode::SATANDVALBLENDING;
    } else if (tcMode == 4) {
        pp->toneCurve.curveMode2 = ToneCurveParams::TcMode::LUMINANCE;
    } else if (tcMode == 5) {
        pp->toneCurve.curveMode2 = ToneCurveParams::TcMode::PERCEPTUAL;
    }

    pp->toneCurve.histmatching = histmatching->get_active();
    pp->toneCurve.fromHistMatching = fromHistMatching;

    pp->toneCurve.saturation = satcurve->getCurve();
}


void ToneCurve::setRaw(bool raw)
{
    disableListener();
    histmatching->set_sensitive(raw);
    enableListener();
}


void ToneCurve::setDefaults(const ProcParams* defParams)
{
}


void ToneCurve::curveChanged(CurveEditor *ce)
{
    if (listener && getEnabled()) {
        if (ce == shape) {
            setHistmatching(false);
            listener->panelChanged(EvToneCurve1, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == shape2) {
            setHistmatching(false);
            listener->panelChanged(EvToneCurve2, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == satcurve) {
            listener->panelChanged(EvSatCurve, M("HISTORY_CUSTOMCURVE"));
        }
    }
}

void ToneCurve::curveMode1Changed ()
{
    //if (listener)  listener->panelChanged (EvToneCurveMode, toneCurveMode->get_active_text());
    if (listener && getEnabled()) {
        //setHistmatching(false);
        Glib::signal_idle().connect (sigc::mem_fun(*this, &ToneCurve::curveMode1Changed_));
    }
}

bool ToneCurve::curveMode1Changed_ ()
{
    if (listener && getEnabled()) {
        listener->panelChanged(EvToneCurveMode1, toneCurveMode->get_active_text());
    }

    return false;
}

void ToneCurve::curveMode2Changed ()
{
    //if (listener)  listener->panelChanged (EvToneCurveMode, toneCurveMode->get_active_text());
    if (listener && getEnabled()) {
        //setHistmatching(false);
        Glib::signal_idle().connect (sigc::mem_fun(*this, &ToneCurve::curveMode2Changed_));
    }
}

bool ToneCurve::curveMode2Changed_ ()
{
    if (listener && getEnabled()) {
        listener->panelChanged(EvToneCurveMode2, toneCurveMode2->get_active_text());
    }

    return false;
}

float ToneCurve::blendPipetteValues(CurveEditor *ce, float chan1, float chan2, float chan3)
{
    // assuming that all the channels are used...
    if (ce == shape) {
        if (toneCurveMode->get_active_row_number() == 4) {
            return chan1 * 0.2126729f + chan2 * 0.7151521f + chan3 * 0.0721750f;
        }
    } else if (ce == shape2) {
        if (toneCurveMode2->get_active_row_number() == 4) {
            return chan1 * 0.2126729f + chan2 * 0.7151521f + chan3 * 0.0721750f;
        }
    }

    return CurveListener::blendPipetteValues(ce, chan1, chan2, chan3);
}


void ToneCurve::enableAll(bool yes)
{
    curveEditorG->set_sensitive(yes);
    toneCurveMode->set_sensitive(yes);
    curveEditorG2->set_sensitive(yes);
    toneCurveMode2->set_sensitive(yes);
    histmatching->set_sensitive(yes);
    satcurveG->set_sensitive(yes);
}


void ToneCurve::trimValues (rtengine::procparams::ProcParams* pp)
{
}


void ToneCurve::updateCurveBackgroundHistogram(
    const LUTu& histToneCurve,
    const LUTu& histLCurve,
    const LUTu& histCCurve,
    const LUTu& histLCAM,
    const LUTu& histCCAM,
    const LUTu& histRed,
    const LUTu& histGreen,
    const LUTu& histBlue,
    const LUTu& histLuma,
    const LUTu& histLRETI
)
{
    shape->updateBackgroundHistogram(histToneCurve);
    satcurve->updateBackgroundHistogram(histToneCurve);
}


void ToneCurve::setHistmatching(bool enabled)
{
    fromHistMatching = enabled;
    if (histmatching->get_active()) {
        histmatchconn.block(true);
        histmatching->set_active(enabled);
        histmatchconn.block(false);
        histmatching->set_inconsistent(false);
    }
}


void ToneCurve::histmatchingToggled()
{
    if (listener) {
        if (histmatching->get_active()) {
            fromHistMatching = false;
            disableListener();
            setEnabled(true);
            enableListener();
            listener->panelChanged(EvHistMatching, M("GENERAL_ENABLED"));
            enableAll(false);
        } else if (getEnabled()) {
            listener->panelChanged(EvHistMatching, M("GENERAL_DISABLED"));
        }
    }
}


void ToneCurve::autoMatchedToneCurveChanged(rtengine::procparams::ToneCurveParams::TcMode curveMode, const std::vector<double>& curve)
{
    nextToneCurveMode = curveMode;
    nextToneCurve = curve;

    idle_register.add(
        [this]() -> bool
        {
            GThreadLock lock;
            disableListener();
            enableAll();

            // toneCurveMode->set_active(rtengine::toUnderlying(nextToneCurveMode));
            shape->setCurve(nextToneCurve);
            // shape2->setCurve({ DCT_Linear });
            shape->openIfNonlinear();

            enableListener();
            fromHistMatching = true;

            return false;
        });
}
