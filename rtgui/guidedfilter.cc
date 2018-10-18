/** -*- C++ -*-
 *  
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2018 Alberto Griggio <alberto.griggio@gmail.com>
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
#include "guidedfilter.h"
#include "eventmapper.h"
#include <iomanip>
#include <cmath>

using namespace rtengine;
using namespace rtengine::procparams;

GuidedFilter::GuidedFilter(): FoldableToolPanel(this, "guidedfilter", M("TP_GUIDED_FILTER_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    EvEnabled = m->newEvent(RGBCURVE, "HISTORY_MSG_GUIDED_FILTER_ENABLED");
    EvSmoothingRadius = m->newEvent(RGBCURVE, "HISTORY_MSG_GUIDED_FILTER_SMOOTHING_RADIUS");
    EvSmoothingEpsilon = m->newEvent(RGBCURVE, "HISTORY_MSG_GUIDED_FILTER_SMOOTHING_EPSILON");
    EvSmoothingIterations = m->newEvent(RGBCURVE, "HISTORY_MSG_GUIDED_FILTER_SMOOTHING_ITERATIONS");
    EvSmoothingLumaBlend = m->newEvent(RGBCURVE, "HISTORY_MSG_GUIDED_FILTER_SMOOTHING_LUMA_BLEND");
    EvSmoothingChromaBlend = m->newEvent(RGBCURVE, "HISTORY_MSG_GUIDED_FILTER_SMOOTHING_CHROMA_BLEND");
    EvDecompRadius = m->newEvent(RGBCURVE, "HISTORY_MSG_GUIDED_FILTER_DECOMP_RADIUS");
    EvDecompEpsilon = m->newEvent(RGBCURVE, "HISTORY_MSG_GUIDED_FILTER_DECOMP_EPSILON");
    EvDecompDetailBoost = m->newEvent(RGBCURVE, "HISTORY_MSG_GUIDED_FILTER_DECOMP_DETAIL_BOOST");
    EvDecompBaseCurve = m->newEvent(RGBCURVE, "HISTORY_MSG_GUIDED_FILTER_DECOMP_BASE_CURVE");

    MyExpander *smoothingExp = Gtk::manage(new MyExpander(false, M("TP_GUIDED_FILTER_SMOOTHING")));
    ToolParamBlock *smoothing = Gtk::manage(new ToolParamBlock());

    MyExpander *decompExp = Gtk::manage(new MyExpander(false, M("TP_GUIDED_FILTER_DECOMP")));
    ToolParamBlock *decomp = Gtk::manage(new ToolParamBlock());

    smoothingRadius = Gtk::manage(new Adjuster(M("TP_GUIDED_FILTER_RADIUS"), 0, 1000, 1, 0));
    smoothingRadius->setLogScale(100, 0);
    smoothingRadius->setAdjusterListener(this);
    smoothing->pack_start(*smoothingRadius);
    
    smoothingEpsilon = Gtk::manage(new Adjuster(M("TP_GUIDED_FILTER_EPSILON"), 1e-6, 0.1, 1e-6, 0.02));
    smoothingEpsilon->setLogScale(1e4, 1e-6);
    smoothingEpsilon->setAdjusterListener(this);
    smoothing->pack_start(*smoothingEpsilon);
    
    smoothingIterations = Gtk::manage(new Adjuster(M("TP_GUIDED_FILTER_SMOOTHING_ITERATIONS"), 1, 4, 1, 1));
    smoothingIterations->setAdjusterListener(this);
    smoothing->pack_start(*smoothingIterations);
    
    smoothingLumaBlend = Gtk::manage(new Adjuster(M("TP_GUIDED_FILTER_SMOOTHING_LUMA_BLEND"), 0, 100, 1, 0));
    smoothingLumaBlend->setAdjusterListener(this);
    smoothing->pack_start(*smoothingLumaBlend);
    
    smoothingChromaBlend = Gtk::manage(new Adjuster(M("TP_GUIDED_FILTER_SMOOTHING_CHROMA_BLEND"), 0, 100, 1, 0));
    smoothingChromaBlend->setAdjusterListener(this);
    smoothing->pack_start(*smoothingChromaBlend);

    decompRadius = Gtk::manage(new Adjuster(M("TP_GUIDED_FILTER_RADIUS"), 0, 1000, 1, 0));
    decompRadius->setLogScale(100, 0);
    decompRadius->setAdjusterListener(this);
    decomp->pack_start(*decompRadius);

    decompEpsilon = Gtk::manage(new Adjuster(M("TP_GUIDED_FILTER_EPSILON"), 1e-6, 0.1, 1e-6, 0.02));
    decompEpsilon->setLogScale(1e4, 1e-6);
    decompEpsilon->setAdjusterListener(this);
    decomp->pack_start(*decompEpsilon);
    
    decompDetailBoost = Gtk::manage(new Adjuster(M("TP_GUIDED_FILTER_DECOMP_DETAIL_BOOST"), -10.0, 10.0, 0.01, 0));
    decompDetailBoost->setLogScale(10, 0, true);
    decompDetailBoost->setAdjusterListener(this);
    decomp->pack_start(*decompDetailBoost);
    
    curveEditorG = new CurveEditorGroup(options.lastToneCurvesDir, M("TP_GUIDED_FILTER_DECOMP_BASE_CURVE"));
    curveEditorG->setCurveListener(this);

    decompBaseCurve = static_cast<DiagonalCurveEditor *>(curveEditorG->addCurve(CT_Diagonal, ""));
    std::vector<GradientMilestone> bottomMilestones;
    bottomMilestones.push_back(GradientMilestone(0., 0., 0., 0.));
    bottomMilestones.push_back(GradientMilestone(1., 1., 1., 1.));    
    decompBaseCurve->setBottomBarBgGradient(bottomMilestones);
    decompBaseCurve->setLeftBarBgGradient(bottomMilestones);
    curveEditorG->curveListComplete();

    decomp->pack_start(*curveEditorG, Gtk::PACK_SHRINK, 2);

    smoothingExp->add(*smoothing, false);
    smoothingExp->setLevel(2);
    decompExp->add(*decomp, false);
    decompExp->setLevel(2);
    pack_start(*smoothingExp);
    pack_start(*decompExp);

    smoothingRadius->delay = options.adjusterMaxDelay;
    smoothingEpsilon->delay = options.adjusterMaxDelay;
    smoothingIterations->delay = options.adjusterMaxDelay;
    smoothingLumaBlend->delay = options.adjusterMaxDelay;
    smoothingChromaBlend->delay = options.adjusterMaxDelay;
    decompRadius->delay = options.adjusterMaxDelay;
    decompEpsilon->delay = options.adjusterMaxDelay;
    decompDetailBoost->delay = options.adjusterMaxDelay;

    show_all_children();
}


void GuidedFilter::read(const ProcParams *pp, const ParamsEdited *pedited)
{
    disableListener();

    if (pedited) {
        set_inconsistent(multiImage && !pedited->guidedfilter.enabled);
        smoothingRadius->setEditedState(pedited->guidedfilter.smoothingRadius ? Edited : UnEdited);
        smoothingEpsilon->setEditedState(pedited->guidedfilter.smoothingEpsilon ? Edited : UnEdited);
        smoothingIterations->setEditedState(pedited->guidedfilter.smoothingIterations ? Edited : UnEdited);
        smoothingLumaBlend->setEditedState(pedited->guidedfilter.smoothingLumaBlend ? Edited : UnEdited);
        smoothingChromaBlend->setEditedState(pedited->guidedfilter.smoothingChromaBlend ? Edited : UnEdited);        
        decompRadius->setEditedState(pedited->guidedfilter.decompRadius ? Edited : UnEdited);
        decompEpsilon->setEditedState(pedited->guidedfilter.decompEpsilon ? Edited : UnEdited);
        decompDetailBoost->setEditedState(pedited->guidedfilter.decompDetailBoost ? Edited : UnEdited);    
        decompBaseCurve->setUnChanged(!pedited->guidedfilter.decompBaseCurve);
    }

    setEnabled(pp->guidedfilter.enabled);
    smoothingRadius->setValue(pp->guidedfilter.smoothingRadius);
    smoothingEpsilon->setValue(pp->guidedfilter.smoothingEpsilon);
    smoothingIterations->setValue(pp->guidedfilter.smoothingIterations);
    smoothingLumaBlend->setValue(pp->guidedfilter.smoothingLumaBlend);
    smoothingChromaBlend->setValue(pp->guidedfilter.smoothingChromaBlend);
    decompRadius->setValue(pp->guidedfilter.decompRadius);
    decompEpsilon->setValue(pp->guidedfilter.decompEpsilon);
    decompDetailBoost->setValue(pp->guidedfilter.decompDetailBoost);    
    decompBaseCurve->setCurve(pp->guidedfilter.decompBaseCurve);

    enableListener();
}


void GuidedFilter::write(ProcParams *pp, ParamsEdited *pedited)
{
    pp->guidedfilter.enabled = getEnabled();
    pp->guidedfilter.smoothingRadius = smoothingRadius->getValue();
    pp->guidedfilter.smoothingEpsilon = smoothingEpsilon->getValue();
    pp->guidedfilter.smoothingIterations = smoothingIterations->getValue();
    pp->guidedfilter.smoothingLumaBlend = smoothingLumaBlend->getValue();
    pp->guidedfilter.smoothingChromaBlend = smoothingChromaBlend->getValue();
    pp->guidedfilter.decompRadius = decompRadius->getValue();
    pp->guidedfilter.decompEpsilon = decompEpsilon->getValue();
    pp->guidedfilter.decompDetailBoost = decompDetailBoost->getValue();
    pp->guidedfilter.decompBaseCurve = decompBaseCurve->getCurve();

    if (pedited) {
        pedited->guidedfilter.enabled = !get_inconsistent();
        pedited->guidedfilter.smoothingRadius = smoothingRadius->getEditedState();
        pedited->guidedfilter.smoothingEpsilon = smoothingEpsilon->getEditedState();
        pedited->guidedfilter.smoothingIterations = smoothingIterations->getEditedState();
        pedited->guidedfilter.smoothingLumaBlend = smoothingLumaBlend->getEditedState();
        pedited->guidedfilter.smoothingChromaBlend = smoothingChromaBlend->getEditedState();
        pedited->guidedfilter.decompRadius = decompRadius->getEditedState();
        pedited->guidedfilter.decompEpsilon = decompEpsilon->getEditedState();
        pedited->guidedfilter.decompDetailBoost = decompDetailBoost->getEditedState();
        pedited->guidedfilter.decompBaseCurve = !decompBaseCurve->isUnChanged();
    }
}

void GuidedFilter::setDefaults(const ProcParams *defParams, const ParamsEdited *pedited)
{
    smoothingRadius->setDefault(defParams->guidedfilter.smoothingRadius);
    smoothingEpsilon->setDefault(defParams->guidedfilter.smoothingEpsilon);
    smoothingIterations->setDefault(defParams->guidedfilter.smoothingIterations);
    smoothingLumaBlend->setDefault(defParams->guidedfilter.smoothingLumaBlend);
    smoothingChromaBlend->setDefault(defParams->guidedfilter.smoothingChromaBlend);
    decompRadius->setDefault(defParams->guidedfilter.decompRadius);
    decompEpsilon->setDefault(defParams->guidedfilter.decompEpsilon);
    decompDetailBoost->setDefault(defParams->guidedfilter.decompDetailBoost);    
    if (pedited) {
        smoothingRadius->setDefaultEditedState(pedited->guidedfilter.smoothingRadius ? Edited : UnEdited);
        smoothingEpsilon->setDefaultEditedState(pedited->guidedfilter.smoothingEpsilon ? Edited : UnEdited);
        smoothingIterations->setDefaultEditedState(pedited->guidedfilter.smoothingIterations ? Edited : UnEdited);
        smoothingLumaBlend->setDefaultEditedState(pedited->guidedfilter.smoothingLumaBlend ? Edited : UnEdited);
        smoothingChromaBlend->setDefaultEditedState(pedited->guidedfilter.smoothingChromaBlend ? Edited : UnEdited);        
        decompRadius->setDefaultEditedState(pedited->guidedfilter.decompRadius ? Edited : UnEdited);
        decompEpsilon->setDefaultEditedState(pedited->guidedfilter.decompEpsilon ? Edited : UnEdited);
        decompDetailBoost->setDefaultEditedState(pedited->guidedfilter.decompDetailBoost ? Edited : UnEdited);    
    } else {
        smoothingRadius->setDefaultEditedState(Irrelevant);
        smoothingEpsilon->setDefaultEditedState(Irrelevant);
        smoothingIterations->setDefaultEditedState(Irrelevant);
        smoothingLumaBlend->setDefaultEditedState(Irrelevant);
        smoothingChromaBlend->setDefaultEditedState(Irrelevant);        
        decompRadius->setDefaultEditedState(Irrelevant);
        decompEpsilon->setDefaultEditedState(Irrelevant);
        decompDetailBoost->setDefaultEditedState(Irrelevant);    
    }
}


void GuidedFilter::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        if (a == smoothingRadius) {
            listener->panelChanged(EvSmoothingRadius, a->getTextValue());
        } else if (a == smoothingEpsilon) {
            listener->panelChanged(EvSmoothingEpsilon, a->getTextValue());
        } else if (a == smoothingIterations) {
            listener->panelChanged(EvSmoothingIterations, a->getTextValue());
        } else if (a == smoothingLumaBlend) {
            listener->panelChanged(EvSmoothingLumaBlend, a->getTextValue());
        } else if (a == smoothingChromaBlend) {
            listener->panelChanged(EvSmoothingChromaBlend, a->getTextValue());
        } else if (a == decompRadius) {
            listener->panelChanged(EvDecompRadius, a->getTextValue());
        } else if (a == decompEpsilon) {
            listener->panelChanged(EvDecompEpsilon, a->getTextValue());
        } else if (a == decompDetailBoost) {
            listener->panelChanged(EvDecompDetailBoost, a->getTextValue());
        }
    }
}


void GuidedFilter::enabledChanged ()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(EvEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(EvEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvEnabled, M("GENERAL_DISABLED"));
        }
    }
}


void GuidedFilter::curveChanged()
{
    if (listener && getEnabled()) {
        listener->panelChanged(EvDecompBaseCurve, M("HISTORY_CUSTOMCURVE"));
    }
}


void GuidedFilter::autoOpenCurve()
{
    decompBaseCurve->openIfNonlinear();
}


void GuidedFilter::setEditProvider(EditDataProvider *provider)
{
    decompBaseCurve->setEditProvider(provider);
}


void GuidedFilter::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);

    smoothingRadius->showEditedCB();
    smoothingEpsilon->showEditedCB();
    smoothingIterations->showEditedCB();
    smoothingLumaBlend->showEditedCB();
    smoothingChromaBlend->showEditedCB();
    decompRadius->showEditedCB();
    decompEpsilon->showEditedCB();
    decompDetailBoost->showEditedCB();
}
