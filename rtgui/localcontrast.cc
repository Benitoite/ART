/** -*- C++ -*-
 *  
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2017 Alberto Griggio <alberto.griggio@gmail.com>
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
#include "localcontrast.h"
#include "eventmapper.h"
#include <iomanip>
#include <cmath>

using namespace rtengine;
using namespace rtengine::procparams;

//-----------------------------------------------------------------------------
// LocalContrastMasksContentProvider
//-----------------------------------------------------------------------------

class LocalContrastMasksContentProvider: public LabMasksContentProvider {
public:
    LocalContrastMasksContentProvider(LocalContrast *parent):
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
        parent_->regionData.push_back(LocalContrastParams::Region());
        return true;
    }

    bool removePressed(int idx) override
    {
        parent_->regionData.erase(parent_->regionData.begin() + idx);
        return true;
    }
    
    bool copyPressed(int idx) override
    {
        parent_->regionData.push_back(parent_->regionData[idx]);
        return true;
    }

    bool resetPressed() override
    {
        parent_->regionData = { LocalContrastParams::Region() };
        parent_->labMasks->setMasks({ Mask() }, -1);
        return true;
    }
    
    bool moveUpPressed(int idx) override
    {
        auto r = parent_->regionData[idx];
        parent_->regionData.erase(parent_->regionData.begin() + idx);
        --idx;
        parent_->regionData.insert(parent_->regionData.begin() + idx, r);
        return true;
    }
    
    bool moveDownPressed(int idx) override
    {
        auto r = parent_->regionData[idx];
        parent_->regionData.erase(parent_->regionData.begin() + idx);
        ++idx;
        parent_->regionData.insert(parent_->regionData.begin() + idx, r);
        return true;
    }

    int getColumnCount() override
    {
        return 1;
    }
    
    Glib::ustring getColumnHeader(int col) override
    {
        return M("TP_LOCALCONTRAST_LIST_TITLE");
    }
    
    Glib::ustring getColumnContent(int col, int row) override
    {
        auto &r = parent_->regionData[row];

        std::string c = "---";
        if (r.curve.size() > 1) {
            std::ostringstream buf;
            const char *sep = "";
            for (size_t i = 1; i+3 < r.curve.size();  i += 4) {
                buf << sep << "("
                    << std::fixed << std::setprecision(2) << r.curve[i]
                    << ", "
                    << std::fixed << std::setprecision(2) << r.curve[i+1]
                    << ")";
                sep = " ";
            }
            c = buf.str();
        }
        return Glib::ustring::compose("%1\n[%2]", c, r.contrast);
    }

    void getEditIDs(EditUniqueID &hcurve, EditUniqueID &ccurve, EditUniqueID &lcurve, EditUniqueID &deltaE) override
    {
        hcurve = EUID_LabMasks_H2;
        ccurve = EUID_LabMasks_C2;
        lcurve = EUID_LabMasks_L2;
        deltaE = EUID_LabMasks_DE2;
    }

private:
    LocalContrast *parent_;
};


//-----------------------------------------------------------------------------
// LocalContrast
//-----------------------------------------------------------------------------

LocalContrast::LocalContrast(): FoldableToolPanel(this, "localcontrast", M("TP_LOCALCONTRAST_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    auto EVENT = DISPLAY;
    EvLocalContrastEnabled = m->newEvent(EVENT, "HISTORY_MSG_LOCALCONTRAST_ENABLED");
    EvLocalContrastContrast = m->newEvent(EVENT, "HISTORY_MSG_LOCALCONTRAST_CONTRAST");
    EvLocalContrastCurve = m->newEvent(EVENT, "HISTORY_MSG_LOCALCONTRAST_CURVE");

    EvList = m->newEvent(EVENT, "HISTORY_MSG_LOCALCONTRAST_LIST");
    EvParametricMask = m->newEvent(EVENT, "HISTORY_MSG_LOCALCONTRAST_PARAMETRICMASK");
    EvHueMask = m->newEvent(EVENT, "HISTORY_MSG_LOCALCONTRAST_HUEMASK");
    EvChromaticityMask = m->newEvent(EVENT, "HISTORY_MSG_LOCALCONTRAST_CHROMATICITYMASK");
    EvLightnessMask = m->newEvent(EVENT, "HISTORY_MSG_LOCALCONTRAST_LIGHTNESSMASK");
    EvMaskBlur = m->newEvent(EVENT, "HISTORY_MSG_LOCALCONTRAST_MASKBLUR");
    EvShowMask = m->newEvent(EVENT, "HISTORY_MSG_LOCALCONTRAST_SHOWMASK");
    EvAreaMask = m->newEvent(EVENT, "HISTORY_MSG_LOCALCONTRAST_AREAMASK");
    EvDeltaEMask = m->newEvent(EVENT, "HISTORY_MSG_LOCALCONTRAST_DELTAEMASK");
    EvContrastThresholdMask = m->newEvent(EVENT, "HISTORY_MSG_LOCALCONTRAST_CONTRASTTHRESHOLDMASK");
    EvDrawnMask = m->newEvent(EVENT, "HISTORY_MSG_LOCALCONTRAST_DRAWNMASK");
    
    box = Gtk::manage(new Gtk::VBox());

    contrast = Gtk::manage(new Adjuster(M("TP_LOCALCONTRAST_CONTRAST"), -100., 100., 0.1, 0.));

    const LocalContrastParams default_params;
    
    cg = Gtk::manage(new CurveEditorGroup(options.lastColorToningCurvesDir, M("TP_LOCALCONTRAST_CURVE"), 0.7));
    cg->setCurveListener(this);
    curve = static_cast<FlatCurveEditor *>(cg->addCurve(CT_Flat, "", nullptr, false, false));
    curve->setIdentityValue(0.);
    curve->setResetCurve(FlatCurveType(default_params.regions[0].curve.at(0)), default_params.regions[0].curve);
    cg->curveListComplete();
    cg->show();
    
    contrast->setAdjusterListener(this);

    box->pack_start(*cg);
    box->pack_start(*contrast);

    labMasksContentProvider.reset(new LocalContrastMasksContentProvider(this));
    labMasks = Gtk::manage(new LabMasksPanel(labMasksContentProvider.get()));
    pack_start(*labMasks, Gtk::PACK_EXPAND_WIDGET, 4);   
}


void LocalContrast::read(const ProcParams *pp)
{
    disableListener();

    setEnabled(pp->localContrast.enabled);
    regionData = pp->localContrast.regions;

    auto m = pp->localContrast.labmasks;
    if (regionData.empty()) {
        regionData.emplace_back(rtengine::LocalContrastParams::Region());
        m.emplace_back(rtengine::Mask());
    }
    labMasks->setMasks(m, pp->localContrast.showMask);

    enableListener();
}


void LocalContrast::write(ProcParams *pp)
{
    pp->localContrast.enabled = getEnabled();
    regionGet(labMasks->getSelected());
    pp->localContrast.regions = regionData;
    labMasks->getMasks(pp->localContrast.labmasks, pp->localContrast.showMask);
    assert(pp->localContrast.regions.size() == pp->localContrast.labmasks.size());
    labMasks->updateSelected();
}

void LocalContrast::setDefaults(const ProcParams *defParams)
{
    // if (defParams->localContrast.regions.size() == 1) {
    //     contrast->setDefault(defParams->localContrast.regions[0].contrast);
    //     curve->setResetCurve(FlatCurveType(defParams->localContrast.regions[0].curve.at(0)), defParams->localContrast.regions[0].curve);
    // }
}

void LocalContrast::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        if (a == contrast) {
            listener->panelChanged(EvLocalContrastContrast, a->getTextValue());
        }
    }
}

void LocalContrast::adjusterAutoToggled(Adjuster* a, bool newval)
{
}

void LocalContrast::enabledChanged ()
{
    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged(EvLocalContrastEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged(EvLocalContrastEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvLocalContrastEnabled, M("GENERAL_DISABLED"));
        }
    }
}


void LocalContrast::curveChanged()
{
    if (listener) {
        listener->panelChanged(EvLocalContrastCurve, M("HISTORY_CUSTOMCURVE"));
    }
}


void LocalContrast::autoOpenCurve()
{
    curve->openIfNonlinear();
}


void LocalContrast::setEditProvider(EditDataProvider *provider)
{
    labMasks->setEditProvider(provider);
}


void LocalContrast::procParamsChanged(
    const rtengine::procparams::ProcParams* params,
    const rtengine::ProcEvent& ev,
    const Glib::ustring& descr,
    const ParamsEdited* paramsEdited)
{
}


void LocalContrast::updateGeometry(int fw, int fh)
{
    labMasks->updateGeometry(fw, fh);
}


void LocalContrast::regionGet(int idx)
{
    if (idx < 0 || size_t(idx) >= regionData.size()) {
        return;
    }
    
    auto &r = regionData[idx];
    r.contrast = contrast->getValue();
    r.curve = curve->getCurve();
}


void LocalContrast::regionShow(int idx)
{
    const bool disable = listener;
    if (disable) {
        disableListener();
    }

    auto &r = regionData[idx];
    contrast->setValue(r.contrast);
    curve->setCurve(r.curve);
    
    if (disable) {
        enableListener();
    }
}


void LocalContrast::setAreaDrawListener(AreaDrawListener *l)
{
    labMasks->setAreaDrawListener(l);
}


void LocalContrast::setDeltaEColorProvider(DeltaEColorProvider *p)
{
    labMasks->setDeltaEColorProvider(p);
}
