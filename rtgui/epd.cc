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
#include "epd.h"
#include "eventmapper.h"
#include <iomanip>
#include <cmath>

using namespace rtengine;
using namespace rtengine::procparams;

//-----------------------------------------------------------------------------
// EPDMasksContentProvider
//-----------------------------------------------------------------------------

class EPDMasksContentProvider: public LabMasksContentProvider {
public:
    EPDMasksContentProvider(EdgePreservingDecompositionUI *parent):
        parent_(parent)
    {
    }

    Gtk::Widget *getWidget() override
    {
        return parent_->box;
    }

    void getEvents(rtengine::ProcEvent &mask_list, rtengine::ProcEvent &h_mask, rtengine::ProcEvent &c_mask, rtengine::ProcEvent &l_mask, rtengine::ProcEvent &blur, rtengine::ProcEvent &show, rtengine::ProcEvent &area_mask) override
    {
        mask_list = parent_->EvList;
        h_mask = parent_->EvHueMask;
        c_mask = parent_->EvChromaticityMask;
        l_mask = parent_->EvLightnessMask;
        blur = parent_->EvMaskBlur;
        show = parent_->EvShowMask;
        area_mask = parent_->EvAreaMask;
    }

    ToolPanelListener *listener() override
    {
        if (parent_->getEnabled()) {
            return parent_->listener;
        }
        return nullptr;
    }

    void selectionChanging(int idx) override
    {
        parent_->regionGet(idx);
    }

    void selectionChanged(int idx) override
    {
        parent_->regionShow(idx);
    }

    bool addPressed() override
    {
        parent_->data.push_back(EPDParams::Region());
        return true;
    }

    bool removePressed(int idx) override
    {
        parent_->data.erase(parent_->data.begin() + idx);
        return true;
    }
    
    bool copyPressed(int idx) override
    {
        parent_->data.push_back(parent_->data[idx]);
        return true;
    }

    bool resetPressed() override
    {
        parent_->data = { EPDParams::Region() };
        parent_->labMasks->setMasks({ LabCorrectionMask() }, -1);
        return true;
    }
    
    bool moveUpPressed(int idx) override
    {
        auto r = parent_->data[idx];
        parent_->data.erase(parent_->data.begin() + idx);
        --idx;
        parent_->data.insert(parent_->data.begin() + idx, r);
        return true;
    }
    
    bool moveDownPressed(int idx) override
    {
        auto r = parent_->data[idx];
        parent_->data.erase(parent_->data.begin() + idx);
        ++idx;
        parent_->data.insert(parent_->data.begin() + idx, r);
        return true;
    }

    int getColumnCount() override
    {
        return 1;
    }
    
    Glib::ustring getColumnHeader(int col) override
    {
        return M("TP_EPD_LIST_TITLE");
    }
    
    Glib::ustring getColumnContent(int col, int row) override
    {
        auto &r = parent_->data[row];

        return Glib::ustring::compose(
            "%1 %2 %3 %4 %5", r.strength, r.gamma, r.edgeStopping, r.scale, r.reweightingIterates); 
    }

    void getEditIDs(EditUniqueID &hcurve, EditUniqueID &ccurve, EditUniqueID &lcurve) override
    {
        hcurve = EUID_LabMasks_H4;
        ccurve = EUID_LabMasks_C4;
        lcurve = EUID_LabMasks_L4;
    }

private:
    EdgePreservingDecompositionUI *parent_;
};


//-----------------------------------------------------------------------------
// EPD
//-----------------------------------------------------------------------------

EdgePreservingDecompositionUI::EdgePreservingDecompositionUI () : FoldableToolPanel(this, "epd", M("TP_EPD_LABEL"), true, true)
{
    auto m = ProcEventMapper::getInstance();
    EvList = m->newEvent(DIRPYREQUALIZER, "HISTORY_MSG_EPD_LIST");
    EvHueMask = m->newEvent(DIRPYREQUALIZER, "HISTORY_MSG_EPD_HUEMASK");
    EvChromaticityMask = m->newEvent(DIRPYREQUALIZER, "HISTORY_MSG_EPD_CHROMATICITYMASK");
    EvLightnessMask = m->newEvent(DIRPYREQUALIZER, "HISTORY_MSG_EPD_LIGHTNESSMASK");
    EvMaskBlur = m->newEvent(DIRPYREQUALIZER, "HISTORY_MSG_EPD_MASKBLUR");
    EvShowMask = m->newEvent(DIRPYREQUALIZER, "HISTORY_MSG_EPD_SHOWMASK");
    EvAreaMask = m->newEvent(DIRPYREQUALIZER, "HISTORY_MSG_EPD_AREAMASK");

    strength = Gtk::manage(new Adjuster (M("TP_EPD_STRENGTH"), -1.0, 2.0, 0.01, 0.5));
    // gamma = Gtk::manage(new Adjuster (M("TP_EPD_GAMMA"), 0.8, 1.5, 0.01, 1.));
    edgeStopping = Gtk::manage(new Adjuster (M("TP_EPD_EDGESTOPPING"), 0.1, 4.0, 0.01, 0.5));
    scale = Gtk::manage(new Adjuster (M("TP_EPD_SCALE"), 0.1, 10.0, 0.01, 0.1));
    // reweightingIterates = Gtk::manage(new Adjuster (M("TP_EPD_REWEIGHTINGITERATES"), 0, 9, 1, 0));

    box = Gtk::manage(new Gtk::VBox());

    strength->setAdjusterListener(this);
    // gamma->setAdjusterListener(this);
    edgeStopping->setAdjusterListener(this);
    scale->setAdjusterListener(this);
    // reweightingIterates->setAdjusterListener(this);

    strength->show();
    // gamma->show();
    edgeStopping->show();
    scale->show();
    // reweightingIterates->show();

    box->pack_start(*strength);
    // box->pack_start(*gamma);
    box->pack_start(*edgeStopping);
    box->pack_start(*scale);
    // box->pack_start(*reweightingIterates);

    labMasksContentProvider.reset(new EPDMasksContentProvider(this));
    labMasks = Gtk::manage(new LabMasksPanel(labMasksContentProvider.get()));
    pack_start(*labMasks, Gtk::PACK_EXPAND_WIDGET, 4);   

    show_all_children();
}

void EdgePreservingDecompositionUI::read(const ProcParams *pp)
{
    disableListener();

    setEnabled(pp->epd.enabled);
    data = pp->epd.regions;
    auto m = pp->epd.labmasks;
    if (data.empty()) {
        data.emplace_back(rtengine::EPDParams::Region());
        m.emplace_back(rtengine::LabCorrectionMask());
    }
    labMasks->updateAreaMaskDefaults(pp);
    labMasks->setMasks(m, pp->epd.showMask);

    enableListener();
}

void EdgePreservingDecompositionUI::write(ProcParams *pp)
{
    pp->epd.enabled = getEnabled();

    regionGet(labMasks->getSelected());
    pp->epd.regions = data;

    labMasks->getMasks(pp->epd.labmasks, pp->epd.showMask);
    assert(pp->epd.regions.size() == pp->epd.labmasks.size());

    labMasks->updateSelected();
}

void EdgePreservingDecompositionUI::setDefaults(const ProcParams *defParams)
{
    strength->setDefault(defParams->epd.regions[0].strength);
    // gamma->setDefault(defParams->epd.regions[0].gamma);
    edgeStopping->setDefault(defParams->epd.regions[0].edgeStopping);
    scale->setDefault(defParams->epd.regions[0].scale);
    // reweightingIterates->setDefault(defParams->epd.regions[0].reweightingIterates);
}

void EdgePreservingDecompositionUI::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        labMasks->setEdited(true);

        if(a == strength) {
            listener->panelChanged(EvEPDStrength, Glib::ustring::format(std::setw(2), std::fixed, std::setprecision(2), a->getValue()));
        // } else if(a == gamma) {
        //     listener->panelChanged(EvEPDgamma, Glib::ustring::format(std::setw(2), std::fixed, std::setprecision(2), a->getValue()));
        } else if(a == edgeStopping) {
            listener->panelChanged(EvEPDEdgeStopping, Glib::ustring::format(std::setw(2), std::fixed, std::setprecision(2), a->getValue()));
        } else if(a == scale) {
            listener->panelChanged(EvEPDScale, Glib::ustring::format(std::setw(2), std::fixed, std::setprecision(2), a->getValue()));
        // } else if(a == reweightingIterates) {
        //     listener->panelChanged(EvEPDReweightingIterates, Glib::ustring::format((int)a->getValue()));
        }
    }
}

void EdgePreservingDecompositionUI::adjusterAutoToggled(Adjuster* a, bool newval)
{
}

void EdgePreservingDecompositionUI::enabledChanged ()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvEPDEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvEPDEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvEPDEnabled, M("GENERAL_DISABLED"));
        }
    }
}


void EdgePreservingDecompositionUI::setEditProvider(EditDataProvider *provider)
{
    labMasks->setEditProvider(provider);
}


void EdgePreservingDecompositionUI::procParamsChanged(
    const rtengine::procparams::ProcParams* params,
    const rtengine::ProcEvent& ev,
    const Glib::ustring& descr,
    const ParamsEdited* paramsEdited)
{
    labMasks->updateAreaMaskDefaults(params);
}


void EdgePreservingDecompositionUI::updateGeometry(int fw, int fh)
{
    labMasks->updateGeometry(fw, fh);
}


void EdgePreservingDecompositionUI::regionGet(int idx)
{
    if (idx < 0 || size_t(idx) >= data.size()) {
        return;
    }
    
    auto &r = data[idx];
    r.strength = strength->getValue();
    // r.gamma = gamma->getValue();
    r.edgeStopping = edgeStopping->getValue();
    r.scale = scale->getValue();
    // r.reweightingIterates = reweightingIterates->getValue();
}


void EdgePreservingDecompositionUI::regionShow(int idx)
{
    const bool disable = listener;
    if (disable) {
        disableListener();
    }

    auto &r = data[idx];
    strength->setValue(r.strength);
    // gamma->setValue(r.gamma);
    edgeStopping->setValue(r.edgeStopping);
    scale->setValue(r.scale);
    // reweightingIterates->setValue(r.reweightingIterates);
    
    if (disable) {
        enableListener();
    }
}
