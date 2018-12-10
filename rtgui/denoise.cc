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
#include "denoise.h"
#include <iomanip>
#include <cmath>
#include "edit.h"
#include "guiutils.h"
#include "eventmapper.h"

using namespace rtengine;
using namespace rtengine::procparams;
extern Options options;

Denoise::Denoise():
    FoldableToolPanel(this, "dirpyrdenoise", M("TP_DIRPYRDENOISE_LABEL"), true, true)
{
    auto m = ProcEventMapper::getInstance();
    EvSmoothingMethod = m->newEvent(ALLNORAW, "HISTORY_MSG_DENOISE_SMOOTHING_METHOD");
    EvGuidedLumaRadius = m->newEvent(ALLNORAW, "HISTORY_MSG_DENOISE_GUIDED_LUMA_RADIUS");
    EvGuidedLumaStrength = m->newEvent(ALLNORAW, "HISTORY_MSG_DENOISE_GUIDED_LUMA_STRENGTH");
    EvGuidedChromaRadius = m->newEvent(ALLNORAW, "HISTORY_MSG_DENOISE_GUIDED_CHROMA_RADIUS");
    EvGuidedChromaStrength = m->newEvent(ALLNORAW, "HISTORY_MSG_DENOISE_GUIDED_CHROMA_STRENGTH");

    std::vector<GradientMilestone> milestones;
    CurveListener::setMulti(true);

    std::vector<double> defaultCurve;

    Gtk::Frame *lumaFrame = Gtk::manage(new Gtk::Frame(M("TP_DIRPYRDENOISE_LUMINANCE_FRAME")));
    lumaFrame->set_label_align(0.025, 0.5);

    Gtk::VBox *lumaVBox = Gtk::manage(new Gtk::VBox());
    lumaVBox->set_spacing(2);

    Gtk::HBox *hb = Gtk::manage(new Gtk::HBox());
    hb->pack_start(*Gtk::manage(new Gtk::Label(M("TP_DIRPYRDENOISE_MAIN_COLORSPACE") + ":")), Gtk::PACK_SHRINK, 1);
    colorSpace = Gtk::manage(new MyComboBoxText());
    colorSpace->append(M("TP_DIRPYRDENOISE_MAIN_COLORSPACE_LAB"));
    colorSpace->append(M("TP_DIRPYRDENOISE_MAIN_COLORSPACE_RGB"));
    colorSpace->set_active(0);
    hb->pack_start(*colorSpace, Gtk::PACK_EXPAND_WIDGET, 1);
    hb->set_tooltip_markup (M("TP_DIRPYRDENOISE_MAIN_COLORSPACE_TOOLTIP"));
    pack_start(*hb, Gtk::PACK_SHRINK, 1);

    hb = Gtk::manage(new Gtk::HBox());
    hb->pack_start(*Gtk::manage(new Gtk::Label(M("TP_DIRPYRDENOISE_MAIN_MODE") + ": ")), Gtk::PACK_SHRINK, 1);
    hb->set_tooltip_markup(M("TP_DIRPYRDENOISE_MAIN_MODE_TOOLTIP"));

    aggressive = Gtk::manage(new MyComboBoxText ());
    aggressive->append(M("TP_DIRPYRDENOISE_MAIN_MODE_CONSERVATIVE"));
    aggressive->append(M("TP_DIRPYRDENOISE_MAIN_MODE_AGGRESSIVE"));
    aggressive->set_active(0);
    hb->pack_start(*aggressive, Gtk::PACK_EXPAND_WIDGET, 1);
    pack_start(*hb, Gtk::PACK_SHRINK, 1);

    gamma = Gtk::manage(new Adjuster(M("TP_DIRPYRDENOISE_MAIN_GAMMA"), 1.0, 3.0, 0.01, 1.7));
    gamma->set_tooltip_text(M("TP_DIRPYRDENOISE_MAIN_GAMMA_TOOLTIP"));
    gamma->setAdjusterListener(this);
    pack_start(*gamma, Gtk::PACK_EXPAND_WIDGET, 1);

    hb = Gtk::manage(new Gtk::HBox());
    hb->pack_start(*Gtk::manage(new Gtk::Label(M("TP_DIRPYRDENOISE_LUMINANCE_CONTROL") + ":")), Gtk::PACK_SHRINK, 1);
    luminanceMethod = Gtk::manage(new MyComboBoxText ());
    luminanceMethod->append(M("GENERAL_SLIDER"));
    luminanceMethod->append(M("CURVEEDITOR_CURVE"));
    luminanceMethod->set_active(0);
    hb->pack_start(*luminanceMethod);
    lumaVBox->pack_start(*hb);
    
    luminance = Gtk::manage(new Adjuster(M("TP_DIRPYRDENOISE_LUMINANCE_SMOOTHING"), 0, 100, 0.01, 0));
    luminanceDetail = Gtk::manage(new Adjuster(M("TP_DIRPYRDENOISE_LUMINANCE_DETAIL"), 0, 100, 0.01, 50));
    luminanceEditorGroup = Gtk::manage(new CurveEditorGroup(options.lastDenoiseCurvesDir, M("TP_DIRPYRDENOISE_LUMINANCE_CURVE"), 0.7));
    luminanceEditorGroup->setCurveListener (this);
    defaultCurve = rtengine::DenoiseParams().luminanceCurve;
    luminanceCurve = static_cast<FlatCurveEditor*>(luminanceEditorGroup->addCurve(CT_Flat, "", nullptr, false, false));
    luminanceCurve->setIdentityValue(0.);
    luminanceCurve->setResetCurve(FlatCurveType(defaultCurve.at(0)), defaultCurve);

    milestones.push_back(GradientMilestone(0., 0., 0., 0.));
    milestones.push_back(GradientMilestone(1., 1., 1., 1.));
    luminanceCurve->setBottomBarBgGradient(milestones);

    milestones.push_back(GradientMilestone(0., 0., 0., 0.));
    milestones.push_back(GradientMilestone(1., 1., 1., 1.));
    luminanceEditorGroup->curveListComplete();
    luminanceEditorGroup->show();

    lumaVBox->pack_start(*luminance);
    lumaVBox->pack_start(*luminanceEditorGroup, Gtk::PACK_SHRINK, 4);
    lumaVBox->pack_start(*luminanceDetail);
    lumaFrame->add(*lumaVBox);
    pack_start(*lumaFrame);

    Gtk::Frame *chromaFrame = Gtk::manage(new Gtk::Frame(M("TP_DIRPYRDENOISE_CHROMINANCE_FRAME")));
    chromaFrame->set_label_align(0.025, 0.5);

    Gtk::VBox *chromaVBox = Gtk::manage(new Gtk::VBox());
    chromaVBox->set_spacing(2);

    hb = Gtk::manage(new Gtk::HBox());
    hb->pack_start(*Gtk::manage (new Gtk::Label (M("TP_DIRPYRDENOISE_CHROMINANCE_METHOD") + ":")), Gtk::PACK_SHRINK, 1);

    chrominanceMethod = Gtk::manage (new MyComboBoxText ());
    chrominanceMethod->append(M("TP_DIRPYRDENOISE_CHROMINANCE_MANUAL"));
    chrominanceMethod->append(M("TP_DIRPYRDENOISE_CHROMINANCE_AUTOGLOBAL"));
    chrominanceMethod->set_active(0);
    chrominanceMethod->set_tooltip_markup (M("TP_DIRPYRDENOISE_CHROMINANCE_METHOD_TOOLTIP"));
    hb->pack_start(*chrominanceMethod);
    chromaVBox->pack_start(*hb);

    chrominance = Gtk::manage(new Adjuster(M("TP_DIRPYRDENOISE_CHROMINANCE_MASTER"), 0, 100, 0.01, 15));
    chrominanceRedGreen = Gtk::manage(new Adjuster(M("TP_DIRPYRDENOISE_CHROMINANCE_REDGREEN"), -100, 100, 0.1, 0));
    chrominanceBlueYellow = Gtk::manage(new Adjuster(M("TP_DIRPYRDENOISE_CHROMINANCE_BLUEYELLOW"), -100, 100, 0.1, 0));

    chrominance->setLogScale(10, 0);
    chrominanceRedGreen->setLogScale(100, 0, true);
    chrominanceBlueYellow->setLogScale(100, 0, true);

    chromaVBox->pack_start(*chrominance);
    chromaVBox->pack_start(*chrominanceRedGreen);
    chromaVBox->pack_start(*chrominanceBlueYellow);

    luminance->setAdjusterListener(this);
    luminanceDetail->setAdjusterListener(this);
    chrominance->setAdjusterListener(this);
    chrominanceRedGreen->setAdjusterListener(this);
    chrominanceBlueYellow->setAdjusterListener(this);

    chrominanceEditorGroup = Gtk::manage(new CurveEditorGroup(options.lastDenoiseCurvesDir, M("TP_DIRPYRDENOISE_CHROMINANCE_CURVE"), 0.7));
    chrominanceEditorGroup->setCurveListener (this);
    defaultCurve = rtengine::DenoiseParams().chrominanceCurve;
    chrominanceCurve = static_cast<FlatCurveEditor*>(chrominanceEditorGroup->addCurve(CT_Flat, "", nullptr, false, false));
    chrominanceCurve->setIdentityValue(0.);
    chrominanceCurve->setResetCurve(FlatCurveType(defaultCurve.at(0)), defaultCurve);
    chrominanceCurve->setTooltip(M("TP_DIRPYRDENOISE_CHROMINANCE_CURVE_TOOLTIP"));
    chrominanceCurve->setBottomBarColorProvider(this, 2);
    chrominanceEditorGroup->curveListComplete();
    chromaVBox->pack_start(*chrominanceEditorGroup, Gtk::PACK_SHRINK, 4);
    chromaFrame->add(*chromaVBox);
    pack_start(*chromaFrame);

    luminance->hide();
    luminanceDetail->show();

    chrominance->show();
    chrominanceRedGreen->show();
    chrominanceBlueYellow->show();

    smoothingEnabled = Gtk::manage(new MyExpander(true, M("TP_DENOISE_SMOOTHING")));
    ToolParamBlock *smoothing = Gtk::manage(new ToolParamBlock());

    smoothingMethod = Gtk::manage(new MyComboBoxText());
    smoothingMethod->append(M("TP_DENOISE_SMOOTHING_MEDIAN"));
    smoothingMethod->append(M("TP_DENOISE_SMOOTHING_GUIDED"));
    hb = Gtk::manage(new Gtk::HBox());
    hb->pack_start(*Gtk::manage (new Gtk::Label(M("TP_DENOISE_SMOOTHING_METHOD") + ":")), Gtk::PACK_SHRINK, 1);
    hb->pack_start(*smoothingMethod);
    smoothing->pack_start(*hb);

    medianBox = Gtk::manage(new Gtk::VBox());
    medianBox->set_spacing(2);

    medianMethod = Gtk::manage(new MyComboBoxText());
    medianMethod->append(M("TP_DIRPYRDENOISE_MEDIAN_METHOD_LUMINANCE"));
    medianMethod->append(M("TP_DIRPYRDENOISE_MEDIAN_METHOD_CHROMINANCE"));
    medianMethod->append(M("TP_DIRPYRDENOISE_MEDIAN_METHOD_WEIGHTED"));
    medianMethod->append(M("TP_DIRPYRDENOISE_MEDIAN_METHOD_LAB"));
    medianMethod->append(M("TP_DIRPYRDENOISE_MEDIAN_METHOD_RGB"));
    medianMethod->set_active(0);
    medianMethod->set_tooltip_text(M("TP_DIRPYRDENOISE_MEDIAN_METHOD_TOOLTIP"));

    medianType = Gtk::manage(new MyComboBoxText());
    medianType->append(M("TP_DIRPYRDENOISE_TYPE_3X3SOFT"));
    medianType->append(M("TP_DIRPYRDENOISE_TYPE_3X3"));
    medianType->append(M("TP_DIRPYRDENOISE_TYPE_5X5SOFT"));
    medianType->append(M("TP_DIRPYRDENOISE_TYPE_5X5"));
    medianType->append(M("TP_DIRPYRDENOISE_TYPE_7X7"));
    medianType->append(M("TP_DIRPYRDENOISE_TYPE_9X9"));
    medianType->set_active(0);
    medianType->set_tooltip_text(M("TP_DIRPYRDENOISE_MEDIAN_TYPE_TOOLTIP"));

    hb = Gtk::manage(new Gtk::HBox());
    hb->pack_start(*Gtk::manage(new Gtk::Label(M("TP_DIRPYRDENOISE_MEDIAN_METHOD") + ":")), Gtk::PACK_SHRINK, 1);
    hb->pack_start(*medianMethod);
    medianBox->pack_start(*hb);

    hb = Gtk::manage(new Gtk::HBox());
    hb->pack_start(*Gtk::manage(new Gtk::Label(M("TP_DIRPYRDENOISE_MEDIAN_TYPE") + ":")), Gtk::PACK_SHRINK, 1);
    hb->pack_start(*medianType);
    medianBox->pack_start(*hb);

    medianIterations = Gtk::manage(new Adjuster(M("TP_DIRPYRDENOISE_MEDIAN_PASSES"), 1.0, 3.0, 1., 1.));
    medianIterations->set_tooltip_text(M("TP_DIRPYRDENOISE_MEDIAN_PASSES_TOOLTIP"));
    medianIterations->setAdjusterListener(this);
    medianIterations->show();
    medianBox->pack_start(*medianIterations);

    smoothing->pack_start(*medianBox);

    guidedBox = Gtk::manage(new Gtk::VBox());

    lumaFrame = Gtk::manage(new Gtk::Frame(M("TP_DIRPYRDENOISE_LUMINANCE_FRAME")));
    lumaFrame->set_label_align(0.025, 0.5);
    lumaVBox = Gtk::manage(new Gtk::VBox());
    lumaVBox->set_spacing(2);
    
    guidedLumaRadius = Gtk::manage(new Adjuster(M("TP_DENOISE_GUIDED_RADIUS"), 0, 100, 1, 0));
    guidedLumaRadius->setLogScale(100, 0);
    guidedLumaRadius->setAdjusterListener(this);
    lumaVBox->pack_start(*guidedLumaRadius);
    
    guidedLumaStrength = Gtk::manage(new Adjuster(M("TP_DENOISE_GUIDED_STRENGTH"), 0, 100, 1, 0));
    guidedLumaStrength->setAdjusterListener(this);
    lumaVBox->pack_start(*guidedLumaStrength);

    lumaFrame->add(*lumaVBox);
    guidedBox->pack_start(*lumaFrame);

    chromaFrame = Gtk::manage(new Gtk::Frame(M("TP_DIRPYRDENOISE_CHROMINANCE_FRAME")));
    chromaFrame->set_label_align(0.025, 0.5);
    chromaVBox = Gtk::manage(new Gtk::VBox());
    chromaVBox->set_spacing(2);
    
    guidedChromaRadius = Gtk::manage(new Adjuster(M("TP_DENOISE_GUIDED_RADIUS"), 0, 100, 1, 0));
    guidedChromaRadius->setLogScale(100, 0);
    guidedChromaRadius->setAdjusterListener(this);
    chromaVBox->pack_start(*guidedChromaRadius);

    guidedChromaStrength = Gtk::manage(new Adjuster(M("TP_DENOISE_GUIDED_STRENGTH"), 0, 100, 1, 0));
    guidedChromaStrength->setAdjusterListener(this);
    chromaVBox->pack_start(*guidedChromaStrength);

    chromaFrame->add(*chromaVBox);
    guidedBox->pack_start(*chromaFrame);

    smoothing->pack_start(*guidedBox);
    smoothingEnabled->add(*smoothing, false);
    smoothingEnabled->setLevel(2);
    pack_start(*smoothingEnabled);

    colorSpace->signal_changed().connect(sigc::mem_fun(*this, &Denoise::colorSpaceChanged));
    aggressive->signal_changed().connect(sigc::mem_fun(*this, &Denoise::aggressiveChanged));
    luminanceMethod->signal_changed().connect(sigc::mem_fun(*this, &Denoise::luminanceMethodChanged));
    chrominanceMethod->signal_changed().connect(sigc::mem_fun(*this, &Denoise::chrominanceMethodChanged));
    medianType->signal_changed().connect(sigc::mem_fun(*this, &Denoise::medianTypeChanged));
    medianMethod->signal_changed().connect(sigc::mem_fun(*this, &Denoise::medianMethodChanged));
    smoothingMethod->signal_changed().connect(sigc::mem_fun(*this, &Denoise::smoothingMethodChanged));
    smoothingEnabled->signal_enabled_toggled().connect(sigc::mem_fun(*this, &Denoise::smoothingEnabledToggled));
}


Denoise::~Denoise ()
{
    idle_register.destroy();
}

void Denoise::chromaChanged (double autchroma, double autred, double autblue)
{
    struct Data {
        Denoise *dn;
        double chroma;
        double red;
        double blue;
    };
    Data *d = new Data{this, autchroma, autred, autblue};
        
    const auto func = [](gpointer data) -> gboolean {
        Data *d = static_cast<Data *>(data);
        d->dn->chromaComputed(d->chroma, d->red, d->blue);
        delete d;
        return false;
    };

    idle_register.add(func, d);
}


bool Denoise::chromaComputed(double chroma, double red, double blue)
{
    disableListener();
    chrominance->setValue(chroma);
    chrominanceRedGreen->setValue(red);
    chrominanceBlueYellow->setValue(blue);
    enableListener();
    return false;
}


void Denoise::read(const ProcParams *pp, const ParamsEdited *pedited)
{
    disableListener();

    setEnabled(pp->denoise.enabled);
    
    colorSpace->set_active(int(pp->denoise.colorSpace));
    aggressive->set_active(pp->denoise.aggressive ? 1 : 0);
    gamma->setValue(pp->denoise.gamma);
    luminanceMethod->set_active(int(pp->denoise.luminanceMethod));
    luminanceMethodChanged();
    luminance->setValue(pp->denoise.luminance);
    luminanceCurve->setCurve(pp->denoise.luminanceCurve);
    luminanceDetail->setValue(pp->denoise.luminanceDetail);

    chrominanceMethod->set_active(int(pp->denoise.chrominanceMethod));
    chrominanceMethodChanged();
    chrominance->setValue(pp->denoise.chrominance);
    chrominanceRedGreen->setValue(pp->denoise.chrominanceRedGreen);
    chrominanceBlueYellow->setValue(pp->denoise.chrominanceBlueYellow);
    chrominanceCurve->setCurve(pp->denoise.chrominanceCurve);

    smoothingEnabled->setEnabled(pp->denoise.smoothingEnabled);
    smoothingMethod->set_active(int(pp->denoise.smoothingMethod));
    smoothingMethodChanged();
    
    medianType->set_active(int(pp->denoise.medianType));
    medianMethod->set_active(int(pp->denoise.medianMethod));
    medianIterations->setValue(pp->denoise.medianIterations);

    guidedLumaRadius->setValue(pp->denoise.guidedLumaRadius);
    guidedLumaStrength->setValue(pp->denoise.guidedLumaStrength);
    guidedChromaRadius->setValue(pp->denoise.guidedChromaRadius);
    guidedChromaStrength->setValue(pp->denoise.guidedChromaStrength);

    if (pedited) {
        if (!pedited->denoise.colorSpace) {
            colorSpace->set_active(2);
        }

        if (!pedited->denoise.aggressive) {
            aggressive->set_active(2);
        }

        if (!pedited->denoise.luminanceMethod) {
            luminanceMethod->set_active(2);
        }

        if (!pedited->denoise.chrominanceMethod) {
            chrominanceMethod->set_active(2);
        }

        if (!pedited->denoise.medianType) {
            medianType->set_active(6);
        }

        if (!pedited->denoise.medianMethod) {
            medianMethod->set_active(5);
        }

        if (!pedited->denoise.smoothingMethod) {
            smoothingMethod->set_active(2);
        }

        luminance->setEditedState(pedited->denoise.luminance ? Edited : UnEdited);
        luminanceDetail->setEditedState(pedited->denoise.luminanceDetail ? Edited : UnEdited);
        chrominance->setEditedState(pedited->denoise.chrominance ? Edited : UnEdited);
        chrominanceRedGreen->setEditedState(pedited->denoise.chrominanceRedGreen ? Edited : UnEdited);
        chrominanceBlueYellow->setEditedState(pedited->denoise.chrominanceBlueYellow ? Edited : UnEdited);

        gamma->setEditedState(pedited->denoise.gamma ? Edited : UnEdited);
        medianIterations->setEditedState(pedited->denoise.medianIterations ? Edited : UnEdited);
        set_inconsistent(multiImage && !pedited->denoise.enabled);
        smoothingEnabled->set_inconsistent(!pedited->denoise.smoothingEnabled);
        chrominanceCurve->setUnChanged(!pedited->denoise.chrominanceCurve);
        luminanceCurve->setUnChanged(!pedited->denoise.luminanceCurve);

        guidedLumaRadius->setEditedState(pedited->denoise.guidedLumaRadius ? Edited : UnEdited);
        guidedLumaStrength->setEditedState(pedited->denoise.guidedLumaStrength ? Edited : UnEdited);
        guidedChromaRadius->setEditedState(pedited->denoise.guidedChromaRadius ? Edited : UnEdited);
        guidedChromaStrength->setEditedState(pedited->denoise.guidedChromaStrength ? Edited : UnEdited);
    }
    enableListener ();
}


void Denoise::setEditProvider(EditDataProvider *provider)
{
    luminanceCurve->setEditProvider(provider);
    chrominanceCurve->setEditProvider(provider);
}


void Denoise::autoOpenCurve ()
{
    luminanceCurve->openIfNonlinear();
    chrominanceCurve->openIfNonlinear();
}


void Denoise::write(ProcParams *pp, ParamsEdited *pedited)
{
    pp->denoise.enabled = getEnabled();
    if (colorSpace->get_active_row_number() < 2) {
        pp->denoise.colorSpace = static_cast<DenoiseParams::ColorSpace>(colorSpace->get_active_row_number());
    }
    if (aggressive->get_active_row_number() < 2) {
        pp->denoise.aggressive = aggressive->get_active_row_number();
    }
    pp->denoise.gamma = gamma->getValue();
    if (luminanceMethod->get_active_row_number() < 3) {
        pp->denoise.luminanceMethod = static_cast<DenoiseParams::LuminanceMethod>(luminanceMethod->get_active_row_number());
    }
    pp->denoise.luminance = luminance->getValue();
    pp->denoise.luminanceDetail = luminanceDetail->getValue();
    pp->denoise.luminanceCurve = luminanceCurve->getCurve();
    if (chrominanceMethod->get_active_row_number() < 2) {
        pp->denoise.chrominanceMethod = static_cast<DenoiseParams::ChrominanceMethod>(chrominanceMethod->get_active_row_number());
    }
    pp->denoise.chrominance = chrominance->getValue();
    pp->denoise.chrominanceRedGreen = chrominanceRedGreen->getValue();
    pp->denoise.chrominanceBlueYellow = chrominanceBlueYellow->getValue();
    pp->denoise.smoothingEnabled = smoothingEnabled->getEnabled();
    if (smoothingMethod->get_active_row_number() < 2) {
        pp->denoise.smoothingMethod = static_cast<DenoiseParams::SmoothingMethod>(smoothingMethod->get_active_row_number());
    }
    if (medianType->get_active_row_number() < 6) {
        pp->denoise.medianType = static_cast<DenoiseParams::MedianType>(medianType->get_active_row_number());
    }
    if (medianMethod->get_active_row_number() < 5) {
        pp->denoise.medianMethod = static_cast<DenoiseParams::MedianMethod>(medianMethod->get_active_row_number());
    }
    pp->denoise.medianIterations = medianIterations->getValue();
    pp->denoise.guidedLumaRadius = guidedLumaRadius->getValue();
    pp->denoise.guidedLumaStrength = guidedLumaStrength->getValue();
    pp->denoise.guidedChromaRadius = guidedChromaRadius->getValue();
    pp->denoise.guidedChromaStrength = guidedChromaStrength->getValue();

    if (pedited) {
        pedited->denoise.enabled  = !get_inconsistent();
        pedited->denoise.colorSpace = colorSpace->get_active_row_number() != 2;
        pedited->denoise.aggressive = aggressive->get_active_row_number() != 2;
        pedited->denoise.gamma = gamma->getEditedState();
        pedited->denoise.luminanceMethod = luminanceMethod->get_active_row_number() != 2;
        pedited->denoise.luminance = luminance->getEditedState();
        pedited->denoise.luminanceCurve = !luminanceCurve->isUnChanged();
        pedited->denoise.luminanceDetail = luminanceDetail->getEditedState();
        pedited->denoise.chrominanceMethod = chrominanceMethod->get_active_row_number() != 2;
        pedited->denoise.chrominance = chrominance->getEditedState();
        pedited->denoise.chrominanceCurve = !chrominanceCurve->isUnChanged();
        pedited->denoise.chrominanceRedGreen = chrominanceRedGreen->getEditedState();
        pedited->denoise.chrominanceBlueYellow = chrominanceBlueYellow->getEditedState();
        pedited->denoise.smoothingEnabled = !smoothingEnabled->get_inconsistent();
        pedited->denoise.smoothingMethod = smoothingMethod->get_active_row_number() != 2;
        pedited->denoise.medianType = medianType->get_active_row_number() != 6;
        pedited->denoise.medianMethod = medianMethod->get_active_row_number() != 5;
        pedited->denoise.medianIterations = medianIterations->getEditedState();
        pedited->denoise.guidedLumaRadius = guidedLumaRadius->getEditedState();
        pedited->denoise.guidedLumaStrength = guidedLumaStrength->getEditedState();
        pedited->denoise.guidedChromaRadius = guidedChromaRadius->getEditedState();
        pedited->denoise.guidedChromaStrength = guidedChromaStrength->getEditedState();
    }
}


void Denoise::curveChanged(CurveEditor* ce)
{

    if (listener && getEnabled()) {
        if (ce == luminanceCurve) {
            listener->panelChanged(EvDPDNLCurve, M("HISTORY_CUSTOMCURVE"));
        }

        if (ce == chrominanceCurve) {
            listener->panelChanged(EvDPDNCCCurve, M("HISTORY_CUSTOMCURVE"));
        }
    }
}


void Denoise::colorSpaceChanged()
{
    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged(EvDPDNmet, colorSpace->get_active_text());
    }
}


void Denoise::luminanceMethodChanged()
{
    if (!batchMode) {
        if (luminanceMethod->get_active_row_number() == 1) {
            luminance->hide();
            luminanceEditorGroup->show();
        } else if (luminanceMethod->get_active_row_number() == 0) {
            luminance->show();
            luminanceEditorGroup->hide();
        }
    }

    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged(EvDPDNLmet, luminanceMethod->get_active_text());
    }
}


void Denoise::chrominanceMethodChanged()
{
    if (!batchMode) {
        if (chrominanceMethod->get_active_row_number() == 0) {
            chrominance->set_sensitive(true);
            chrominanceRedGreen->set_sensitive(true);
            chrominanceBlueYellow->set_sensitive(true);
        } else if (chrominanceMethod->get_active_row_number() == 1) {
            chrominance->set_sensitive(false);
            chrominanceRedGreen->set_sensitive(false);
            chrominanceBlueYellow->set_sensitive(false);
        }
    }

    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged(EvDPDNCmet, chrominanceMethod->get_active_text());
    }
}


void Denoise::aggressiveChanged()
{
    if (listener && (multiImage || getEnabled()) ) {
        listener->panelChanged(EvDPDNsmet, aggressive->get_active_text());
    }
}


void Denoise::medianTypeChanged()
{
    if (listener && (multiImage || getEnabled())  && smoothingEnabled->getEnabled()) {
        listener->panelChanged(EvDPDNmedmet, medianType->get_active_text());
    }
}


void Denoise::medianMethodChanged()
{
    if (listener && (multiImage || getEnabled())  && smoothingEnabled->getEnabled()) {
        listener->panelChanged(EvDPDNmetmed, medianMethod->get_active_text());
    }
}


void Denoise::setDefaults(const ProcParams *defParams, const ParamsEdited *pedited)
{
    luminance->setDefault(defParams->denoise.luminance);
    luminanceDetail->setDefault(defParams->denoise.luminanceDetail);
    chrominance->setDefault(defParams->denoise.chrominance);
    chrominanceRedGreen->setDefault(defParams->denoise.chrominanceRedGreen);
    chrominanceBlueYellow->setDefault(defParams->denoise.chrominanceBlueYellow);
    gamma->setDefault(defParams->denoise.gamma);
    medianIterations->setDefault(defParams->denoise.medianIterations);
    guidedLumaRadius->setDefault(defParams->denoise.guidedLumaRadius);
    guidedLumaStrength->setDefault(defParams->denoise.guidedLumaStrength);
    guidedChromaRadius->setDefault(defParams->denoise.guidedChromaRadius);
    guidedChromaStrength->setDefault(defParams->denoise.guidedChromaStrength);

    if (pedited) {
        luminance->setDefaultEditedState(pedited->denoise.luminance ? Edited : UnEdited);
        luminanceDetail->setDefaultEditedState(pedited->denoise.luminanceDetail ? Edited : UnEdited);
        chrominance->setDefaultEditedState(pedited->denoise.chrominance ? Edited : UnEdited);
        chrominanceRedGreen->setDefaultEditedState(pedited->denoise.chrominanceRedGreen ? Edited : UnEdited);
        chrominanceBlueYellow->setDefaultEditedState(pedited->denoise.chrominanceBlueYellow ? Edited : UnEdited);
        gamma->setDefaultEditedState(pedited->denoise.gamma ? Edited : UnEdited);
        medianIterations->setDefaultEditedState(pedited->denoise.medianIterations ? Edited : UnEdited);
        guidedLumaRadius->setDefaultEditedState(pedited->denoise.guidedLumaRadius ? Edited : UnEdited);
        guidedLumaStrength->setDefaultEditedState(pedited->denoise.guidedLumaStrength ? Edited : UnEdited);
        guidedChromaRadius->setDefaultEditedState(pedited->denoise.guidedChromaRadius ? Edited : UnEdited);
        guidedChromaStrength->setDefaultEditedState(pedited->denoise.guidedChromaStrength ? Edited : UnEdited);
    } else {
        luminance->setDefaultEditedState(Irrelevant);
        luminanceDetail->setDefaultEditedState(Irrelevant);
        chrominance->setDefaultEditedState(Irrelevant);
        chrominanceRedGreen->setDefaultEditedState(Irrelevant);
        chrominanceBlueYellow->setDefaultEditedState(Irrelevant);
        gamma->setDefaultEditedState(Irrelevant);
        medianIterations->setDefaultEditedState(Irrelevant);
        guidedLumaRadius->setDefaultEditedState(Irrelevant);
        guidedLumaStrength->setDefaultEditedState(Irrelevant);
        guidedChromaRadius->setDefaultEditedState(Irrelevant);
        guidedChromaStrength->setDefaultEditedState(Irrelevant);
    }
}


void Denoise::adjusterChanged(Adjuster* a, double newval)
{
    const Glib::ustring costr = Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), a->getValue());

    if (listener && getEnabled()) {
        if (a == luminanceDetail) {
            listener->panelChanged(EvDPDNLdetail, costr);
        } else if (a == luminance) {
            listener->panelChanged(EvDPDNLuma, costr);
        } else if (a == chrominance) {
            listener->panelChanged(EvDPDNChroma, costr);
        } else if (a == chrominanceRedGreen) {
            listener->panelChanged(EvDPDNredchro, costr);
        } else if (a == chrominanceBlueYellow) {
            listener->panelChanged(EvDPDNbluechro, costr);
        } else if (a == gamma) {
            listener->panelChanged(EvDPDNGamma, costr);
        } else if (a == medianIterations && smoothingEnabled->getEnabled()) {
            listener->panelChanged(EvDPDNpasses, costr);
        } else if (a == guidedLumaRadius && smoothingEnabled->getEnabled()) {
            listener->panelChanged(EvGuidedLumaRadius, costr);
        } else if (a == guidedLumaStrength && smoothingEnabled->getEnabled()) {
            listener->panelChanged(EvGuidedLumaStrength, costr);
        } else if (a == guidedChromaRadius && smoothingEnabled->getEnabled()) {
            listener->panelChanged(EvGuidedChromaRadius, costr);
        } else if (a == guidedChromaStrength && smoothingEnabled->getEnabled()) {
            listener->panelChanged(EvGuidedChromaStrength, costr);
        } 
    }
}

void Denoise::adjusterAutoToggled(Adjuster* a, bool newval)
{
}


void Denoise::enabledChanged()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(EvDPDNEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(EvDPDNEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvDPDNEnabled, M("GENERAL_DISABLED"));
        }
    }
}


void Denoise::smoothingEnabledToggled()
{
    if (listener) {
        if (smoothingEnabled->getEnabled()) {
            listener->panelChanged(EvDPDNmedian, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvDPDNmedian, M("GENERAL_DISABLED"));
        }
    }
}


void Denoise::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);
    luminance->showEditedCB ();
    luminanceDetail->showEditedCB ();
    chrominance->showEditedCB ();
    chrominanceRedGreen->showEditedCB ();
    chrominanceBlueYellow->showEditedCB ();
    gamma->showEditedCB ();
    medianIterations->showEditedCB ();
    luminanceEditorGroup->setBatchMode (batchMode);
    chrominanceEditorGroup->setBatchMode (batchMode);
    guidedLumaRadius->showEditedCB();
    guidedLumaStrength->showEditedCB();
    guidedChromaRadius->showEditedCB();
    guidedChromaStrength->showEditedCB();

    colorSpace->append (M("GENERAL_UNCHANGED"));
    aggressive->append (M("GENERAL_UNCHANGED"));
    luminanceMethod->append (M("GENERAL_UNCHANGED"));
    chrominanceMethod->append (M("GENERAL_UNCHANGED"));
    medianMethod->append (M("GENERAL_UNCHANGED"));
    medianType->append (M("GENERAL_UNCHANGED"));
    smoothingMethod->append(M("GENERAL_UNCHANGED"));
}


void Denoise::setAdjusterBehavior (bool lumaadd, bool lumdetadd, bool chromaadd, bool chromaredadd, bool chromablueadd, bool gammaadd, bool passesadd)
{
    luminance->setAddMode(lumaadd);
    luminanceDetail->setAddMode(lumdetadd);
    chrominance->setAddMode(chromaadd);
    chrominanceRedGreen->setAddMode(chromaredadd);
    chrominanceBlueYellow->setAddMode(chromablueadd);
    gamma->setAddMode(gammaadd);
    medianIterations->setAddMode(passesadd);
}


void Denoise::colorForValue (double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller *caller)
{

    float R = 0.f, G = 0.f, B = 0.f;

    if (elemType == ColorCaller::CCET_VERTICAL_BAR) {
        valY = 0.5;
    }

    if (callerId == 1) {         // ch - main curve

        Color::hsv2rgb01(float(valX), float(valY), 0.5f, R, G, B);
    } else if (callerId == 2) {  // cc - bottom bar

        float value = (1.f - 0.7f) * float(valX) + 0.7f;
        // whole hue range
        // Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
        Color::hsv2rgb01(float(valY), float(valX), value, R, G, B);
    } else if (callerId == 3) {  // lc - bottom bar

        float value = (1.f - 0.7f) * float(valX) + 0.7f;
        // whole hue range
        // Y axis / from 0.15 up to 0.75 (arbitrary values; was 0.45 before)
        Color::hsv2rgb01(float(valY), float(valX), value, R, G, B);
    }

    else if (callerId == 4) {    // LH - bottom bar
        Color::hsv2rgb01(float(valX), 0.5f, float(valY), R, G, B);
    } else if (callerId == 5) {  // HH - bottom bar
        float h = float((valY - 0.5) * 0.3 + valX);

        if (h > 1.0f) {
            h -= 1.0f;
        } else if (h < 0.0f) {
            h += 1.0f;
        }

        Color::hsv2rgb01(h, 0.5f, 0.5f, R, G, B);
    }

    caller->ccRed = double(R);
    caller->ccGreen = double(G);
    caller->ccBlue = double(B);
}


void Denoise::trimValues (rtengine::procparams::ProcParams* pp)
{
    luminance->trimValue(pp->denoise.luminance);
    luminanceDetail->trimValue(pp->denoise.luminanceDetail);
    chrominance->trimValue(pp->denoise.chrominance);
    chrominanceRedGreen->trimValue(pp->denoise.chrominanceRedGreen);
    chrominanceBlueYellow->trimValue(pp->denoise.chrominanceBlueYellow);
    gamma->trimValue(pp->denoise.gamma);
    medianIterations->trimValue(pp->denoise.medianIterations);
    guidedLumaRadius->trimValue(pp->denoise.guidedLumaRadius);
    guidedLumaStrength->trimValue(pp->denoise.guidedLumaStrength);
    guidedChromaRadius->trimValue(pp->denoise.guidedChromaRadius);
    guidedChromaStrength->trimValue(pp->denoise.guidedChromaStrength);
}


void Denoise::smoothingMethodChanged()
{
    if (!batchMode) {
        if (smoothingMethod->get_active_row_number() == 0) {
            medianBox->show();
            guidedBox->hide();
        } else if (smoothingMethod->get_active_row_number() == 1) {
            medianBox->hide();
            guidedBox->show();
        }
    }
    if (listener) {
        listener->panelChanged(EvSmoothingMethod, M("GENERAL_CHANGED"));
    }
}
