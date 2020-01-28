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

    void getEvents(rtengine::ProcEvent &mask_list, rtengine::ProcEvent &h_mask, rtengine::ProcEvent &c_mask, rtengine::ProcEvent &l_mask, rtengine::ProcEvent &blur, rtengine::ProcEvent &show, rtengine::ProcEvent &area_mask, rtengine::ProcEvent &deltaE_mask, rtengine::ProcEvent &contrastThreshold_mask) override
    {
        mask_list = parent_->EvList;
        h_mask = parent_->EvHueMask;
        c_mask = parent_->EvChromaticityMask;
        l_mask = parent_->EvLightnessMask;
        blur = parent_->EvMaskBlur;
        show = parent_->EvShowMask;
        area_mask = parent_->EvAreaMask;
        deltaE_mask = parent_->EvDeltaEMask;
        contrastThreshold_mask = parent_->EvContrastThresholdMask;
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
        parent_->data.push_back(rtengine::ColorCorrectionParams::Region());
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
        parent_->data = { ColorCorrectionParams::Region() };
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
        switch (r.mode) {
        case rtengine::procparams::ColorCorrectionParams::Mode::RGB: {
            const auto lbl =
                [](const std::array<double, 3> &v) -> Glib::ustring
                {
                    return Glib::ustring::compose("{%1,%2,%3}", v[0], v[1], v[2]);
                };
            return Glib::ustring::compose("RGB s=%1 o=%2\np=%3 P=%4", lbl(r.slope), lbl(r.offset), lbl(r.power), lbl(r.pivot));
        }   break;
        case rtengine::procparams::ColorCorrectionParams::Mode::HSL: {
            const auto lbl =
                [](const std::array<double, 3> &v) -> Glib::ustring
                {
                    return Glib::ustring::compose("{s=%1,o=%2,p=%3}", v[0], v[1], v[2]);
                };
            return Glib::ustring::compose("HSL H=%1\nS=%2 L=%3", lbl(r.hue), lbl(r.sat), lbl(r.factor));
        }   break;
        default: {
            const auto round_ab = [](float v) -> float
                                  {
                                      return int(v * 1000) / 1000.f;
                                  };
            return Glib::ustring::compose("a=%1 b=%2 S=%3\ns=%4 o=%5 p=%6 P=%7", round_ab(r.a), round_ab(r.b), r.saturation, r.slope[0], r.offset[0], r.power[0], r.pivot[0]);            
        }   break;
        }
    }

    void getEditIDs(EditUniqueID &hcurve, EditUniqueID &ccurve, EditUniqueID &lcurve, EditUniqueID &deltaE) override
    {
        hcurve = EUID_LabMasks_H1;
        ccurve = EUID_LabMasks_C1;
        lcurve = EUID_LabMasks_L1;
        deltaE = EUID_LabMasks_DE1;
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
    auto EVENT = LUMINANCECURVE | M_LUMACURVE;
    EvEnabled = m->newEvent(EVENT, "HISTORY_MSG_COLORCORRECTION_ENABLED");
    EvAB = m->newEvent(EVENT, "HISTORY_MSG_COLORCORRECTION_AB");
    EvSaturation = m->newEvent(EVENT, "HISTORY_MSG_COLORCORRECTION_SATURATION");
    EvLightness = m->newEvent(EVENT, "HISTORY_MSG_COLORCORRECTION_LIGHTNESS");
    EvSlope = m->newEvent(EVENT, "HISTORY_MSG_COLORCORRECTION_SLOPE");
    EvOffset = m->newEvent(EVENT, "HISTORY_MSG_COLORCORRECTION_OFFSET");
    EvPower = m->newEvent(EVENT, "HISTORY_MSG_COLORCORRECTION_POWER");
    EvPivot = m->newEvent(EVENT, "HISTORY_MSG_COLORCORRECTION_PIVOT");
    EvMode = m->newEvent(EVENT, "HISTORY_MSG_COLORCORRECTION_MODE");

    EvList = m->newEvent(EVENT, "HISTORY_MSG_COLORCORRECTION_LIST");
    EvHueMask = m->newEvent(EVENT, "HISTORY_MSG_COLORCORRECTION_HUEMASK");
    EvChromaticityMask = m->newEvent(EVENT, "HISTORY_MSG_COLORCORRECTION_CHROMATICITYMASK");
    EvLightnessMask = m->newEvent(EVENT, "HISTORY_MSG_COLORCORRECTION_LIGHTNESSMASK");
    EvMaskBlur = m->newEvent(EVENT, "HISTORY_MSG_COLORCORRECTION_MASKBLUR");
    EvShowMask = m->newEvent(EVENT, "HISTORY_MSG_COLORCORRECTION_SHOWMASK");
    EvAreaMask = m->newEvent(EVENT, "HISTORY_MSG_COLORCORRECTION_AREAMASK");
    EvDeltaEMask = m->newEvent(EVENT, "HISTORY_MSG_COLORCORRECTION_DELTAEMASK");
    EvContrastThresholdMask = m->newEvent(EVENT, "HISTORY_MSG_COLORCORRECTION_CONTRASTTHRESHOLDMASK");

    box = Gtk::manage(new Gtk::VBox());

    Gtk::HBox *hb = Gtk::manage(new Gtk::HBox());
    mode = Gtk::manage(new MyComboBoxText());
    mode->append(M("TP_COLORCORRECTION_MODE_COMBINED"));
    mode->append(M("TP_COLORCORRECTION_MODE_RGBCHANNELS"));
    mode->append(M("TP_COLORCORRECTION_MODE_HSL"));
    mode->set_active(0);
    mode->signal_changed().connect(sigc::mem_fun(*this, &ColorCorrection::modeChanged));
    
    hb->pack_start(*Gtk::manage(new Gtk::Label(M("TP_COLORCORRECTION_MODE") + ": ")), Gtk::PACK_SHRINK);
    hb->pack_start(*mode);
    box->pack_start(*hb);
    
    box_combined = Gtk::manage(new Gtk::VBox());
    box_rgb = Gtk::manage(new Gtk::VBox());

    gridAB = Gtk::manage(new LabGrid(EvAB, M("TP_COLORCORRECTION_ABVALUES"), false));
    box_combined->pack_start(*gridAB);

    saturation = Gtk::manage(new Adjuster(M("TP_COLORCORRECTION_SATURATION"), -100, 100, 1, 0));
    saturation->setAdjusterListener(this);
    box_combined->pack_start(*saturation);

    slope = Gtk::manage(new Adjuster(M("TP_COLORCORRECTION_SLOPE"), 0.01, 10.0, 0.001, 1));
    slope->setLogScale(10, 1, true);
    slope->setAdjusterListener(this);
    box_combined->pack_start(*slope);
    offset = Gtk::manage(new Adjuster(M("TP_COLORCORRECTION_OFFSET"), -0.1, 0.1, 0.001, 0));
    offset->setAdjusterListener(this);
    box_combined->pack_start(*offset);
    power = Gtk::manage(new Adjuster(M("TP_COLORCORRECTION_POWER"), 0.1, 4.0, 0.001, 1));
    power->setAdjusterListener(this);
    power->setLogScale(10, 1, true);
    box_combined->pack_start(*power);
    pivot = Gtk::manage(new Adjuster(M("TP_COLORCORRECTION_PIVOT"), 0.001, 2.0, 0.001, 1));
    pivot->setAdjusterListener(this);
    pivot->setLogScale(100, 0.18, true);
    box_combined->pack_start(*pivot);

    for (int c = 0; c < 3; ++c) {
        const char *chan = (c == 0 ? "R" : (c == 1 ? "G" : "B"));
        const char *img = (c == 0 ? "red" : (c == 1 ? "green" : "blue"));
        Gtk::Frame *f = Gtk::manage(new Gtk::Frame(""));
        Gtk::HBox *lbl = Gtk::manage(new Gtk::HBox());
        lbl->pack_start(*Gtk::manage(new RTImage(std::string("circle-") + img + "-small.png")), Gtk::PACK_SHRINK, 2);
        lbl->pack_start(*Gtk::manage(new Gtk::Label(M(std::string("TP_COLORCORRECTION_CHANNEL_") + chan))));
        f->set_label_align(0.025, 0.5);
        f->set_label_widget(*lbl);
        Gtk::VBox *vb = Gtk::manage(new Gtk::VBox());
        vb->set_spacing(2);
        vb->set_border_width(2);
    
        slope_rgb[c] = Gtk::manage(new Adjuster(M("TP_COLORCORRECTION_SLOPE"), 0.01, 10.0, 0.001, 1));
        slope_rgb[c]->setLogScale(10, 1, true);
        slope_rgb[c]->setAdjusterListener(this);
        vb->pack_start(*slope_rgb[c]);
        offset_rgb[c] = Gtk::manage(new Adjuster(M("TP_COLORCORRECTION_OFFSET"), -0.1, 0.1, 0.001, 0));
        offset_rgb[c]->setAdjusterListener(this);
        vb->pack_start(*offset_rgb[c]);
        power_rgb[c] = Gtk::manage(new Adjuster(M("TP_COLORCORRECTION_POWER"), 0.1, 4.0, 0.001, 1));
        power_rgb[c]->setAdjusterListener(this);
        power_rgb[c]->setLogScale(10, 1, true);
        vb->pack_start(*power_rgb[c]);
        pivot_rgb[c] = Gtk::manage(new Adjuster(M("TP_COLORCORRECTION_PIVOT"), 0.001, 2.0, 0.001, 1));
        pivot_rgb[c]->setAdjusterListener(this);
        pivot_rgb[c]->setLogScale(100, 0.18, true);
        vb->pack_start(*pivot_rgb[c]);

        f->add(*vb);
        box_rgb->pack_start(*f);
    }

    sync_rgb_sliders = Gtk::manage(new Gtk::CheckButton(M("TP_COLORCORRECTION_SYNC_SLIDERS")));
    box_rgb->pack_start(*sync_rgb_sliders, Gtk::PACK_SHRINK, 4);

    box_hsl = Gtk::manage(new Gtk::VBox());
    for (int c = 0; c < 3; ++c) {
        const char *chan = (c == 0 ? "SLOPE" : (c == 1 ? "OFFSET" : "POWER"));
        Gtk::Frame *f = Gtk::manage(new Gtk::Frame(M(std::string("TP_COLORCORRECTION_") + chan)));
        Gtk::VBox *vb = Gtk::manage(new Gtk::VBox());
        vb->set_spacing(2);
        vb->set_border_width(2);
    
        huesat[c] = Gtk::manage(new ThresholdAdjuster(M("TP_COLORCORRECTION_HUESAT"), 0., 100., 0., M("TP_COLORCORRECTION_S"), 1., 0., 360., 0., M("TP_COLORCORRECTION_H"), 1., nullptr, false));
        huesat[c]->setAdjusterListener(this);
        huesat[c]->setBgColorProvider(this, c+1);
        huesat[c]->setUpdatePolicy(RTUP_DYNAMIC);
        vb->pack_start(*huesat[c]);

        lfactor[c] = Gtk::manage(new Adjuster(M("TP_COLORCORRECTION_L"), -50.0, 50.0, c == 1 ? 0.01 : 0.1, 0));
        if (c == 1) {
            lfactor[c]->setLogScale(100, 0, true);
        }
        lfactor[c]->setAdjusterListener(this);
        vb->pack_start(*lfactor[c]);

        f->add(*vb);
        box_hsl->pack_start(*f);
    }

        
    saturation->delay = options.adjusterMaxDelay;
    slope->delay = options.adjusterMaxDelay;
    offset->delay = options.adjusterMaxDelay;
    power->delay = options.adjusterMaxDelay;
    pivot->delay = options.adjusterMaxDelay;
    for (int c = 0; c < 3; ++c) {
        slope_rgb[c]->delay = options.adjusterMaxDelay;
        offset_rgb[c]->delay = options.adjusterMaxDelay;
        power_rgb[c]->delay = options.adjusterMaxDelay;
        pivot_rgb[c]->delay = options.adjusterMaxDelay;
        lfactor[c]->delay = options.adjusterMaxDelay;
    }

    labMasksContentProvider.reset(new ColorCorrectionMasksContentProvider(this));
    labMasks = Gtk::manage(new LabMasksPanel(labMasksContentProvider.get()));
    pack_start(*labMasks, Gtk::PACK_EXPAND_WIDGET, 4);

    box->pack_start(*box_combined);
    box->pack_start(*box_rgb);
    box->pack_start(*box_hsl);

    show_all_children();
}


void ColorCorrection::read(const ProcParams *pp)
{
    disableListener();

    setEnabled(pp->colorcorrection.enabled);
    data = pp->colorcorrection.regions;
    auto m = pp->colorcorrection.labmasks;
    if (data.empty()) {
        data.emplace_back(rtengine::ColorCorrectionParams::Region());
        m.emplace_back(rtengine::LabCorrectionMask());
    }
    labMasks->setMasks(m, pp->colorcorrection.showMask);

    modeChanged();

    enableListener();
}


void ColorCorrection::write(ProcParams *pp)
{
    pp->colorcorrection.enabled = getEnabled();

    regionGet(labMasks->getSelected());
    pp->colorcorrection.regions = data;

    labMasks->getMasks(pp->colorcorrection.labmasks, pp->colorcorrection.showMask);
    assert(pp->colorcorrection.regions.size() == pp->colorcorrection.labmasks.size());

    labMasks->updateSelected();
}


void ColorCorrection::setDefaults(const ProcParams *defParams)
{
    saturation->setDefault(defParams->colorcorrection.regions[0].saturation);
    slope->setDefault(defParams->colorcorrection.regions[0].slope[0]);
    offset->setDefault(defParams->colorcorrection.regions[0].offset[0]);
    power->setDefault(defParams->colorcorrection.regions[0].power[0]);
    pivot->setDefault(defParams->colorcorrection.regions[0].pivot[0]);
    for (int c = 0; c < 3; ++c) {
        slope_rgb[c]->setDefault(defParams->colorcorrection.regions[0].slope[c]);
        offset_rgb[c]->setDefault(defParams->colorcorrection.regions[0].offset[c]);
        power_rgb[c]->setDefault(defParams->colorcorrection.regions[0].power[c]);
        pivot_rgb[c]->setDefault(defParams->colorcorrection.regions[0].pivot[c]);
    }
}


void ColorCorrection::adjusterChanged(Adjuster* a, double newval)
{
    rtengine::ProcEvent evt;
    Glib::ustring msg = a->getTextValue();
    if (a == saturation) {
        evt = EvSaturation;
    } else if (a == slope) {
        evt = EvSlope;
    } else if (a == offset) {
        evt = EvOffset;
    } else if (a == power) {
        evt = EvPower;
    } else if (a == pivot) {
        evt = EvPivot;
    } else {
        Adjuster **targets = nullptr;
        for (int c = 0; c < 3; ++c) {
            if (a == slope_rgb[c]) {
                evt = EvSlope;
                targets = slope_rgb;
                break;
            } else if (a == offset_rgb[c]) {
                evt = EvOffset;
                targets = offset_rgb;
                break;
            } else if (a == power_rgb[c]) {
                evt = EvPower;
                targets = power_rgb;
                break;
            } else if (a == pivot_rgb[c]) {
                evt = EvPivot;
                targets = pivot_rgb;
                break;
            } else if (a == lfactor[c]) {
                evt = c == 0 ? EvSlope : (c == 1 ? EvOffset : EvPower);
                break;
            }
        }
        if (targets && sync_rgb_sliders->get_active()) {
            for (int c = 0; c < 3; ++c) {
                targets[c]->setAdjusterListener(nullptr);
                targets[c]->setValue(newval);
                targets[c]->setAdjusterListener(this);
            }
        }
        if (evt != 0) {
            if (targets) {
                msg = Glib::ustring::compose("R=%1 G=%2 B=%3",
                                             targets[0]->getTextValue(),
                                             targets[1]->getTextValue(),
                                             targets[2]->getTextValue());
            } else {
                double h, s, l;
                if (evt == EvSlope) {
                    l = lfactor[0]->getValue();
                    huesat[0]->getValue(s, h);
                } else if (evt == EvOffset) {
                    l = lfactor[1]->getValue();
                    huesat[1]->getValue(s, h);
                } else {
                    l = lfactor[2]->getValue();
                    huesat[2]->getValue(s, h);
                }
                msg = Glib::ustring::compose("H=%1 S=%2 L=%3", h, s, l);
            }
        }
    }
        
    if (listener && getEnabled() && evt != 0) {
        listener->panelChanged(evt, msg);
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
    switch (mode->get_active_row_number()) {
    case 1:
        r.mode = rtengine::procparams::ColorCorrectionParams::Mode::RGB;
        break;
    case 2:
        r.mode = rtengine::procparams::ColorCorrectionParams::Mode::HSL;
        break;
    default:
        r.mode = rtengine::procparams::ColorCorrectionParams::Mode::YUV;
    }
    if (r.mode != rtengine::procparams::ColorCorrectionParams::Mode::RGB) {
        double la, lb;
        gridAB->getParams(la, lb, r.a, r.b);
        r.saturation = saturation->getValue();
        for (int c = 0; c < 3; ++c) {
            r.slope[c] = slope->getValue();
            r.offset[c] = offset->getValue();
            r.power[c] = power->getValue();
            r.pivot[c] = pivot->getValue();
        }
    } else {
        r.saturation = 0.f;
        r.a = r.b = 0.f;
        for (int c = 0; c < 3; ++c) {
            r.slope[c] = slope_rgb[c]->getValue();
            r.offset[c] = offset_rgb[c]->getValue();
            r.power[c] = power_rgb[c]->getValue();
            r.pivot[c] = pivot_rgb[c]->getValue();
        }
    }
    for (int c = 0; c < 3; ++c) {
        double b, t;
        huesat[c]->getValue(b, t);
        r.hue[c] = t;
        r.sat[c] = b;
        r.factor[c] = lfactor[c]->getValue();
    }
}


void ColorCorrection::regionShow(int idx)
{
    const bool disable = listener;
    if (disable) {
        disableListener();
    }
    auto &r = data[idx];
    switch (r.mode) {
    case rtengine::procparams::ColorCorrectionParams::Mode::HSL:
        mode->set_active(2);
        break;
    case rtengine::procparams::ColorCorrectionParams::Mode::RGB:
        mode->set_active(1);
        break;
    default:
        mode->set_active(0);
    }
    if (r.mode != rtengine::procparams::ColorCorrectionParams::Mode::RGB) {
        gridAB->setParams(0, 0, r.a, r.b, false);
        saturation->setValue(r.saturation);
        slope->setValue(r.slope[0]);
        offset->setValue(r.offset[0]);
        power->setValue(r.power[0]);
        pivot->setValue(r.pivot[0]);
    } else {
        for (int c = 0; c < 3; ++c) {
            slope_rgb[c]->setValue(r.slope[c]);
            offset_rgb[c]->setValue(r.offset[c]);
            power_rgb[c]->setValue(r.power[c]);
            pivot_rgb[c]->setValue(r.pivot[c]);
        }
    }
    for (int c = 0; c < 3; ++c) {
        huesat[c]->setValue(r.sat[c], r.hue[c]);
        lfactor[c]->setValue(r.factor[c]);
    }
    modeChanged();
    if (disable) {
        enableListener();
    }
}


void ColorCorrection::modeChanged()
{
    removeIfThere(box, box_combined);
    removeIfThere(box, box_rgb);
    removeIfThere(box, box_hsl);
    if (mode->get_active_row_number() == 0) {
        box->pack_start(*box_combined);
    } else if (mode->get_active_row_number() == 1) {
        box->pack_start(*box_rgb);
    } else {
        box->pack_start(*box_hsl);
    }
    if (listener && getEnabled()) {
        labMasks->setEdited(true);        
        listener->panelChanged(EvMode, mode->get_active_text());
    }
}


void ColorCorrection::setListener(ToolPanelListener *tpl)
{
     ToolPanel::setListener(tpl);
     gridAB->setListener(tpl);
}


void ColorCorrection::setAreaDrawListener(AreaDrawListener *l)
{
    labMasks->setAreaDrawListener(l);
}


void ColorCorrection::setDeltaEColorProvider(DeltaEColorProvider *p)
{
    labMasks->setDeltaEColorProvider(p);
}


void ColorCorrection::adjusterChanged(ThresholdAdjuster *a, double newBottom, double newTop)
{
    if (listener && getEnabled()) {
        rtengine::ProcEvent evt;
        float lval = 0.f;
        for (int c = 0; c < 3; ++c) {
            if (a == huesat[c]) {
                evt = (c == 0 ? EvSlope : (c == 1 ? EvOffset : EvPower));
                lval = lfactor[c]->getValue();
                break;
            }
        }
        listener->panelChanged(evt, Glib::ustring::compose("H=%1 S=%2 L=%3", newTop, newBottom, lval));
    }
}


void ColorCorrection::colorForValue(double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller *caller)
{

    float R = 0.f, G = 0.f, B = 0.f;

    const auto to_hue =
        [](float x) -> float
        {
            x -= 0.05f;
            if (x < 0.f) {
                x += 1.f;
            }
            return x;
        };

    if (callerId >= 1 && callerId <= 3) {
        if (valY <= 0.5) {
            // the hue range
            Color::hsv2rgb01(to_hue(valX), 1.0f, 0.65f, R, G, B);
        } else {
            // the strength applied to the current hue
            double strength, hue;
            huesat[callerId-1]->getValue(strength, hue);
            Color::hsv2rgb01(to_hue(hue / 360.f), 1.f, 1.f, R, G, B);
            const double gray = 0.46;
            R = (gray * (1.0 - valX)) + R * valX;
            G = (gray * (1.0 - valX)) + G * valX;
            B = (gray * (1.0 - valX)) + B * valX;
        }
    }

    caller->ccRed = double(R);
    caller->ccGreen = double(G);
    caller->ccBlue = double(B);
}
