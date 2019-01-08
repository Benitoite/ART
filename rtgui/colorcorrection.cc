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
#include "colorcorrection.h"
#include "eventmapper.h"
#include <iomanip>
#include <cmath>

using namespace rtengine;
using namespace rtengine::procparams;


//-----------------------------------------------------------------------------
// ColorCorrectionMasksContentProvider
//-----------------------------------------------------------------------------

class ColorCorrectionMasksContentProvider: public LabMasksContentProvider {
public:
    ColorCorrectionMasksContentProvider(ColorCorrection *parent):
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
        parent_->data.push_back(rtengine::ColorCorrectionParams::LabCorrectionRegion());
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

    bool resetPressed() override
    {
        parent_->data = { ColorCorrectionParams::LabCorrectionRegion() };
        parent_->labMasks->setMasks({ LabCorrectionMask() }, -1);
        return true;
    }

    int getColumnCount() override
    {
        return 1;
    }
    
    Glib::ustring getColumnHeader(int col) override
    {
        return M("TP_COLORCORRECTION_LIST_TITLE");
    }
    
    Glib::ustring getColumnContent(int col, int row) override
    {
        auto &r = parent_->data[row];
        const char *ch = "";
        switch (r.channel) {
        case rtengine::ColorCorrectionParams::LabCorrectionRegion::CHAN_R:
            ch = " [R]"; break;
        case rtengine::ColorCorrectionParams::LabCorrectionRegion::CHAN_G:
            ch = " [G]"; break;
        case rtengine::ColorCorrectionParams::LabCorrectionRegion::CHAN_B:
            ch = " [B]"; break;
        default:
            ch = "";
        }

        const auto round_ab = [](float v) -> float
            {
                return int(v * 1000) / 1000.f;
            };
        return Glib::ustring::compose("a=%1 b=%2 S=%3\ns=%4 o=%5 p=%6%7", round_ab(r.a), round_ab(r.b), r.saturation, r.slope, r.offset, r.power, ch);
    }

    void getEditIDs(EditUniqueID &hcurve, EditUniqueID &ccurve, EditUniqueID &lcurve) override
    {
        hcurve = EUID_LabMasks_H1;
        ccurve = EUID_LabMasks_C1;
        lcurve = EUID_LabMasks_L1;
    }

private:
    ColorCorrection *parent_;
};

//-----------------------------------------------------------------------------
// ColorCorrection
//-----------------------------------------------------------------------------

ColorCorrection::ColorCorrection(): FoldableToolPanel(this, "colorcorrection", M("TP_COLORCORRECTION_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    EvEnabled = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORCORRECTION_ENABLED");
    EvAB = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORCORRECTION_AB");
    EvSaturation = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORCORRECTION_SATURATION");
    EvLightness = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORCORRECTION_LIGHTNESS");
    EvSlope = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORCORRECTION_SLOPE");
    EvOffset = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORCORRECTION_OFFSET");
    EvPower = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORCORRECTION_POWER");
    EvChannel = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORCORRECTION_CHANNEL");

    EvList = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORCORRECTION_LIST");
    EvHueMask = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORCORRECTION_HUEMASK");
    EvChromaticityMask = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORCORRECTION_CHROMATICITYMASK");
    EvLightnessMask = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORCORRECTION_LIGHTNESSMASK");
    EvMaskBlur = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORCORRECTION_MASKBLUR");
    EvShowMask = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORCORRECTION_SHOWMASK");
    EvAreaMask = m->newEvent(LUMINANCECURVE, "HISTORY_MSG_COLORCORRECTION_AREAMASK");

    box = Gtk::manage(new Gtk::VBox());

    gridAB = Gtk::manage(new LabGrid(EvAB, M("TP_COLORCORRECTION_ABVALUES"), false));
    box->pack_start(*gridAB);

    saturation = Gtk::manage(new Adjuster(M("TP_COLORCORRECTION_SATURATION"), -100, 100, 1, 0));
    saturation->setAdjusterListener(this);
    box->pack_start(*saturation);

    slope = Gtk::manage(new Adjuster(M("TP_COLORCORRECTION_SLOPE"), 0.1, 4.0, 0.001, 1));
    slope->setLogScale(4, 0.1);
    slope->setAdjusterListener(this);
    box->pack_start(*slope);
    offset = Gtk::manage(new Adjuster(M("TP_COLORCORRECTION_OFFSET"), -0.1, 0.1, 0.001, 0));
    offset->setAdjusterListener(this);
    box->pack_start(*offset);
    power = Gtk::manage(new Adjuster(M("TP_COLORCORRECTION_POWER"), 0.1, 4.0, 0.001, 1));
    power->setAdjusterListener(this);
    power->setLogScale(4, 0.1);
    box->pack_start(*power);

    Gtk::HBox *hb = Gtk::manage(new Gtk::HBox());
    channel = Gtk::manage(new MyComboBoxText());
    channel->append(M("TP_COLORCORRECTION_CHANNEL_ALL"));
    channel->append(M("TP_COLORCORRECTION_CHANNEL_R"));
    channel->append(M("TP_COLORCORRECTION_CHANNEL_G"));
    channel->append(M("TP_COLORCORRECTION_CHANNEL_B"));
    channel->set_active(0);
    channel->signal_changed().connect(sigc::mem_fun(*this, &ColorCorrection::channelChanged));
    
    hb->pack_start(*Gtk::manage(new Gtk::Label(M("TP_COLORCORRECTION_CHANNEL") + ": ")), Gtk::PACK_SHRINK);
    hb->pack_start(*channel);
    box->pack_start(*hb);

    saturation->delay = options.adjusterMaxDelay;
    slope->delay = options.adjusterMaxDelay;
    offset->delay = options.adjusterMaxDelay;
    power->delay = options.adjusterMaxDelay;

    labMasksContentProvider.reset(new ColorCorrectionMasksContentProvider(this));
    labMasks = Gtk::manage(new LabMasksPanel(labMasksContentProvider.get()));
    pack_start(*labMasks, Gtk::PACK_EXPAND_WIDGET, 4);

    show_all_children();
}


void ColorCorrection::read(const ProcParams *pp, const ParamsEdited *pedited)
{
    disableListener();

    setEnabled(pp->colorcorrection.enabled);
    data = pp->colorcorrection.regions;
    auto m = pp->colorcorrection.labmasks;
    if (data.empty()) {
        data.emplace_back(rtengine::ColorCorrectionParams::LabCorrectionRegion());
        m.emplace_back(rtengine::LabCorrectionMask());
    }
    labMasks->updateAreaMaskDefaults(pp);
    labMasks->setMasks(m, pp->colorcorrection.showMask);

    if (pedited) {
        set_inconsistent(multiImage && !pedited->colorcorrection.enabled);
        if (pedited->colorcorrection.regions) {
            labMasks->setEdited(true);
        } else {
            labMasks->setEdited(false);
        }
    }

    enableListener();
}


void ColorCorrection::write(ProcParams *pp, ParamsEdited *pedited)
{
    pp->colorcorrection.enabled = getEnabled();

    regionGet(labMasks->getSelected());
    pp->colorcorrection.regions = data;

    labMasks->getMasks(pp->colorcorrection.labmasks, pp->colorcorrection.showMask);
    assert(pp->colorcorrection.regions.size() == pp->colorcorrection.labmasks.size());

    labMasks->updateSelected();
        
    if (pedited) {
        pedited->colorcorrection.enabled = !get_inconsistent();
        pedited->colorcorrection.regions = labMasks->getEdited();
    }
}


void ColorCorrection::setDefaults(const ProcParams *defParams, const ParamsEdited *pedited)
{
    saturation->setDefault(defParams->colorcorrection.regions[0].saturation);
    slope->setDefault(defParams->colorcorrection.regions[0].slope);
    offset->setDefault(defParams->colorcorrection.regions[0].offset);
    power->setDefault(defParams->colorcorrection.regions[0].power);
    if (pedited) {
        saturation->setDefaultEditedState(pedited->colorcorrection.regions ? Edited : UnEdited);
        slope->setDefaultEditedState(pedited->colorcorrection.regions ? Edited : UnEdited);
        offset->setDefaultEditedState(pedited->colorcorrection.regions ? Edited : UnEdited);
        power->setDefaultEditedState(pedited->colorcorrection.regions ? Edited : UnEdited);
    } else {
        saturation->setDefaultEditedState(Irrelevant);
        slope->setDefaultEditedState(Irrelevant);
        offset->setDefaultEditedState(Irrelevant);
        power->setDefaultEditedState(Irrelevant);
    }
}


void ColorCorrection::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        labMasks->setEdited(true);
        
        if (a == saturation) {
            gridAB->setEdited(true);
            listener->panelChanged(EvSaturation, a->getTextValue());
        } else if (a == slope) {
            gridAB->setEdited(true);
            listener->panelChanged(EvSlope, a->getTextValue());
        } else if (a == offset) {
            gridAB->setEdited(true);
            listener->panelChanged(EvOffset, a->getTextValue());
        } else if (a == power) {
            gridAB->setEdited(true);
            listener->panelChanged(EvPower, a->getTextValue());
        }
    }
}


void ColorCorrection::enabledChanged ()
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


void ColorCorrection::setBatchMode(bool batchMode)
{
    ToolPanel::setBatchMode(batchMode);
    labMasks->setBatchMode();
}


void ColorCorrection::setEditProvider(EditDataProvider *provider)
{
    labMasks->setEditProvider(provider);
}


void ColorCorrection::procParamsChanged(
    const rtengine::procparams::ProcParams* params,
    const rtengine::ProcEvent& ev,
    const Glib::ustring& descr,
    const ParamsEdited* paramsEdited)
{
    labMasks->updateAreaMaskDefaults(params);
}


void ColorCorrection::updateGeometry(int fw, int fh)
{
    labMasks->updateGeometry(fw, fh);
}


void ColorCorrection::regionGet(int idx)
{
    if (idx < 0 || size_t(idx) >= data.size()) {
        return;
    }
    
    auto &r = data[idx];
    double la, lb;
    gridAB->getParams(la, lb, r.a, r.b);
    r.saturation = saturation->getValue();
    r.slope = slope->getValue();
    r.offset = offset->getValue();
    r.power = power->getValue();
    r.channel = channel->get_active_row_number() - 1;
}


void ColorCorrection::regionShow(int idx)
{
    const bool disable = listener;
    if (disable) {
        disableListener();
    }
    auto &r = data[idx];
    gridAB->setParams(0, 0, r.a, r.b, false);
    saturation->setValue(r.saturation);
    slope->setValue(r.slope);
    offset->setValue(r.offset);
    power->setValue(r.power);
    channel->set_active(r.channel+1);
    if (disable) {
        enableListener();
    }
}


void ColorCorrection::channelChanged()
{
    if (listener && (multiImage || getEnabled()) ) {
        labMasks->setEdited(true);        
        listener->panelChanged(EvChannel, channel->get_active_text());
    }
}


void ColorCorrection::setListener(ToolPanelListener *tpl)
{
     ToolPanel::setListener(tpl);
     gridAB->setListener(tpl);
}