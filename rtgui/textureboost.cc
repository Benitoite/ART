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
#include "textureboost.h"
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
    EPDMasksContentProvider(TextureBoost *parent):
        parent_(parent)
    {
    }

    Gtk::Widget *getWidget() override
    {
        return parent_->box;
    }

    void getEvents(rtengine::ProcEvent &mask_list, rtengine::ProcEvent &parametric_mask, rtengine::ProcEvent &h_mask, rtengine::ProcEvent &c_mask, rtengine::ProcEvent &l_mask, rtengine::ProcEvent &blur, rtengine::ProcEvent &show, rtengine::ProcEvent &area_mask, rtengine::ProcEvent &deltaE_mask, rtengine::ProcEvent &contrastThreshold_mask, rtengine::ProcEvent &drawn_mask) override
    {
        mask_list = parent_->EvList;
        parametric_mask = parent_->EvParametricMask;
        h_mask = parent_->EvHueMask;
        c_mask = parent_->EvChromaticityMask;
        l_mask = parent_->EvLightnessMask;
        blur = parent_->EvMaskBlur;
        show = parent_->EvShowMask;
        area_mask = parent_->EvAreaMask;
        deltaE_mask = parent_->EvDeltaEMask;
        contrastThreshold_mask = parent_->EvContrastThresholdMask;
        drawn_mask = parent_->EvDrawnMask;
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
        parent_->data.push_back(TextureBoostParams::Region());
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
        parent_->data = { TextureBoostParams::Region() };
        parent_->labMasks->setMasks({ Mask() }, -1);
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
            "%1 %2 %3", r.strength, r.edgeStopping, r.scale); 
    }

    void getEditIDs(EditUniqueID &hcurve, EditUniqueID &ccurve, EditUniqueID &lcurve, EditUniqueID &deltaE) override
    {
        hcurve = EUID_LabMasks_H4;
        ccurve = EUID_LabMasks_C4;
        lcurve = EUID_LabMasks_L4;
        deltaE = EUID_LabMasks_DE4;
    }

private:
    TextureBoost *parent_;
};


//-----------------------------------------------------------------------------
// EPD
//-----------------------------------------------------------------------------

TextureBoost::TextureBoost () : FoldableToolPanel(this, "epd", M("TP_EPD_LABEL"), true, true)
{
    auto m = ProcEventMapper::getInstance();
    EvList = m->newEvent(DISPLAY, "HISTORY_MSG_EPD_LIST");
    EvParametricMask = m->newEvent(DISPLAY, "HISTORY_MSG_EPD_PARAMETRICMASK");
    EvHueMask = m->newEvent(DISPLAY, "HISTORY_MSG_EPD_HUEMASK");
    EvChromaticityMask = m->newEvent(DISPLAY, "HISTORY_MSG_EPD_CHROMATICITYMASK");
    EvLightnessMask = m->newEvent(DISPLAY, "HISTORY_MSG_EPD_LIGHTNESSMASK");
    EvMaskBlur = m->newEvent(DISPLAY, "HISTORY_MSG_EPD_MASKBLUR");
    EvShowMask = m->newEvent(DISPLAY, "HISTORY_MSG_EPD_SHOWMASK");
    EvAreaMask = m->newEvent(DISPLAY, "HISTORY_MSG_EPD_AREAMASK");
    EvDeltaEMask = m->newEvent(DISPLAY, "HISTORY_MSG_EPD_DELTAEMASK");
    EvContrastThresholdMask = m->newEvent(DISPLAY, "HISTORY_MSG_EPD_CONTRASTTHRESHOLDMASK");
    EvDrawnMask = m->newEvent(DISPLAY, "HISTORY_MSG_EPD_DRAWNMASK");

    strength = Gtk::manage(new Adjuster (M("TP_EPD_STRENGTH"), -1.0, 2.0, 0.01, 0.5));
    edgeStopping = Gtk::manage(new Adjuster (M("TP_EPD_EDGESTOPPING"), 0.1, 4.0, 0.01, 0.5));
    scale = Gtk::manage(new Adjuster (M("TP_EPD_SCALE"), 0.1, 10.0, 0.01, 0.1));

    box = Gtk::manage(new Gtk::VBox());

    strength->setAdjusterListener(this);
    edgeStopping->setAdjusterListener(this);
    scale->setAdjusterListener(this);

    strength->show();
    edgeStopping->show();
    scale->show();

    box->pack_start(*strength);
    box->pack_start(*edgeStopping);
    box->pack_start(*scale);

    labMasksContentProvider.reset(new EPDMasksContentProvider(this));
    labMasks = Gtk::manage(new LabMasksPanel(labMasksContentProvider.get()));
    pack_start(*labMasks, Gtk::PACK_EXPAND_WIDGET, 4);   

    show_all_children();
}

void TextureBoost::read(const ProcParams *pp)
{
    disableListener();

    setEnabled(pp->textureBoost.enabled);
    data = pp->textureBoost.regions;
    auto m = pp->textureBoost.labmasks;
    if (data.empty()) {
        data.emplace_back(rtengine::TextureBoostParams::Region());
        m.emplace_back(rtengine::Mask());
    }
    labMasks->setMasks(m, pp->textureBoost.showMask);

    enableListener();
}

void TextureBoost::write(ProcParams *pp)
{
    pp->textureBoost.enabled = getEnabled();

    regionGet(labMasks->getSelected());
    pp->textureBoost.regions = data;

    labMasks->getMasks(pp->textureBoost.labmasks, pp->textureBoost.showMask);
    assert(pp->textureBoost.regions.size() == pp->textureBoost.labmasks.size());

    labMasks->updateSelected();
}

void TextureBoost::setDefaults(const ProcParams *defParams)
{
    strength->setDefault(defParams->textureBoost.regions[0].strength);
    edgeStopping->setDefault(defParams->textureBoost.regions[0].edgeStopping);
    scale->setDefault(defParams->textureBoost.regions[0].scale);
}

void TextureBoost::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        labMasks->setEdited(true);

        if(a == strength) {
            listener->panelChanged(EvEPDStrength, Glib::ustring::format(std::setw(2), std::fixed, std::setprecision(2), a->getValue()));
        } else if(a == edgeStopping) {
            listener->panelChanged(EvEPDEdgeStopping, Glib::ustring::format(std::setw(2), std::fixed, std::setprecision(2), a->getValue()));
        } else if(a == scale) {
            listener->panelChanged(EvEPDScale, Glib::ustring::format(std::setw(2), std::fixed, std::setprecision(2), a->getValue()));
        }
    }
}

void TextureBoost::adjusterAutoToggled(Adjuster* a, bool newval)
{
}

void TextureBoost::enabledChanged ()
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


void TextureBoost::setEditProvider(EditDataProvider *provider)
{
    labMasks->setEditProvider(provider);
}


void TextureBoost::procParamsChanged(
    const rtengine::procparams::ProcParams* params,
    const rtengine::ProcEvent& ev,
    const Glib::ustring& descr,
    const ParamsEdited* paramsEdited)
{
}


void TextureBoost::updateGeometry(int fw, int fh)
{
    labMasks->updateGeometry(fw, fh);
}


void TextureBoost::regionGet(int idx)
{
    if (idx < 0 || size_t(idx) >= data.size()) {
        return;
    }
    
    auto &r = data[idx];
    r.strength = strength->getValue();
    r.edgeStopping = edgeStopping->getValue();
    r.scale = scale->getValue();
}


void TextureBoost::regionShow(int idx)
{
    const bool disable = listener;
    if (disable) {
        disableListener();
    }

    auto &r = data[idx];
    strength->setValue(r.strength);
    edgeStopping->setValue(r.edgeStopping);
    scale->setValue(r.scale);
    
    if (disable) {
        enableListener();
    }
}


void TextureBoost::setAreaDrawListener(AreaDrawListener *l)
{
    labMasks->setAreaDrawListener(l);
}


void TextureBoost::setDeltaEColorProvider(DeltaEColorProvider *p)
{
    labMasks->setDeltaEColorProvider(p);
}
