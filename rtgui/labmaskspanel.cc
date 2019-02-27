/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright 2018 Alberto Griggio <alberto.griggio@gmail.com>
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

#include "labmaskspanel.h"

namespace {

constexpr int ID_HUE_MASK = 5;

inline bool hasMask(const std::vector<double> &dflt, const std::vector<double> &mask)
{
    return !(mask.empty() || mask[0] == FCT_Linear || mask == dflt);
}

} // namespace


LabMasksPanel::LabMasksPanel(LabMasksContentProvider *cp):
    Gtk::VBox(),
    cp_(cp),
    masks_(),
    selected_(0),
    area_shape_index_(0),
    listEdited(false)
{
    Gtk::Widget *child = cp_->getWidget();
    cp_->getEvents(EvMaskList, EvHMask, EvCMask, EvLMask, EvMaskBlur, EvShowMask, EvAreaMask);
    
    CurveListener::setMulti(true);
    
    const auto add_button =
        [](Gtk::Button *btn, Gtk::Box *box, int h=20) -> void
        {
            setExpandAlignProperties(btn, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_START);
            btn->set_relief(Gtk::RELIEF_NONE);
            btn->get_style_context()->add_class(GTK_STYLE_CLASS_FLAT);
            btn->set_can_focus(false);
            if (h > 0) {
                btn->set_size_request(-1, h);
            }
            box->pack_start(*btn, false, false);
        };
    
    int n = cp_->getColumnCount();
    list = Gtk::manage(new Gtk::ListViewText(n+2));
    list->set_size_request(-1, 150);
    list->set_can_focus(false);
    list->set_column_title(0, "#");
    for (int i = 0; i < n; ++i) {
        list->set_column_title(i+1, cp_->getColumnHeader(i));
    }
    list->set_column_title(n+1, M("TP_LABMASKS_MASK"));
    list->set_activate_on_single_click(true);
    selectionConn = list->get_selection()->signal_changed().connect(sigc::mem_fun(this, &LabMasksPanel::onSelectionChanged));
    Gtk::HBox *hb = Gtk::manage(new Gtk::HBox());
    Gtk::ScrolledWindow *scroll = Gtk::manage(new Gtk::ScrolledWindow());
    scroll->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_NEVER);
    scroll->add(*list);
    hb->pack_start(*scroll, Gtk::PACK_EXPAND_WIDGET);
    Gtk::VBox *vb = Gtk::manage(new Gtk::VBox());
    reset = Gtk::manage(new Gtk::Button());
    reset->add(*Gtk::manage(new RTImage("undo-small.png")));
    reset->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onResetPressed));
    add_button(reset, vb);
    add = Gtk::manage(new Gtk::Button());
    add->add(*Gtk::manage(new RTImage("add-small.png")));
    add->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onAddPressed));
    add_button(add, vb);
    remove = Gtk::manage(new Gtk::Button());
    remove->add(*Gtk::manage(new RTImage("remove-small.png")));
    remove->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onRemovePressed));
    add_button(remove, vb);
    up = Gtk::manage(new Gtk::Button());
    up->add(*Gtk::manage(new RTImage("arrow-up-small.png")));
    up->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onUpPressed));
    add_button(up, vb);
    down = Gtk::manage(new Gtk::Button());
    down->add(*Gtk::manage(new RTImage("arrow-down-small.png")));
    down->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onDownPressed));
    add_button(down, vb);
    copy = Gtk::manage(new Gtk::Button());
    copy->add(*Gtk::manage(new RTImage("copy-small.png")));
    copy->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onCopyPressed));
    add_button(copy, vb);
    hb->pack_start(*vb, Gtk::PACK_SHRINK);
    pack_start(*hb, true, true);

    pack_start(*Gtk::manage(new Gtk::HSeparator()));
    pack_start(*child);
    pack_start(*Gtk::manage(new Gtk::HSeparator()));

    maskEditorGroup = Gtk::manage(new CurveEditorGroup(options.lastColorToningCurvesDir, M("TP_LABMASKS_MASK"), 0.7));
    maskEditorGroup->setCurveListener(this);

    rtengine::LabCorrectionMask default_params;

    EditUniqueID eh, ec, el;
    cp_->getEditIDs(eh, ec, el);

    hueMask = static_cast<FlatCurveEditor *>(maskEditorGroup->addCurve(CT_Flat, M("TP_LABMASKS_HUE"), nullptr, false, true));
    hueMask->setIdentityValue(0.);
    hueMask->setResetCurve(FlatCurveType(default_params.hueMask[0]), default_params.hueMask);
    hueMask->setCurveColorProvider(this, ID_HUE_MASK);
    hueMask->setBottomBarColorProvider(this, ID_HUE_MASK);
    hueMask->setEditID(eh, BT_SINGLEPLANE_FLOAT);

    chromaticityMask = static_cast<FlatCurveEditor *>(maskEditorGroup->addCurve(CT_Flat, M("TP_LABMASKS_CHROMATICITY"), nullptr, false, false));
    chromaticityMask->setIdentityValue(0.);
    chromaticityMask->setResetCurve(FlatCurveType(default_params.chromaticityMask[0]), default_params.chromaticityMask);
    chromaticityMask->setBottomBarColorProvider(this, ID_HUE_MASK+1);
    chromaticityMask->setEditID(ec, BT_SINGLEPLANE_FLOAT);

    lightnessMask = static_cast<FlatCurveEditor *>(maskEditorGroup->addCurve(CT_Flat, M("TP_LABMASKS_LIGHTNESS"), nullptr, false, false));
    lightnessMask->setIdentityValue(0.);
    lightnessMask->setResetCurve(FlatCurveType(default_params.lightnessMask[0]), default_params.lightnessMask);
    std::vector<GradientMilestone> milestones = {
        GradientMilestone(0., 0., 0., 0.),
        GradientMilestone(1., 1., 1., 1.)
    };
    lightnessMask->setBottomBarBgGradient(milestones);
    lightnessMask->setEditID(el, BT_SINGLEPLANE_FLOAT);

    maskEditorGroup->curveListComplete();
    maskEditorGroup->show();
    pack_start(*maskEditorGroup, Gtk::PACK_SHRINK, 2);

    maskBlur = Gtk::manage(new Adjuster(M("TP_LABMASKS_BLUR"), -10, 100, 0.1, 0));
    maskBlur->setLogScale(10, 0);
    maskBlur->setAdjusterListener(this);
    pack_start(*maskBlur);

    areaMask = Gtk::manage(new MyExpander(true, M("TP_LABMASKS_AREA")));
    areaMask->signal_enabled_toggled().connect(sigc::mem_fun(*this, &LabMasksPanel::onAreaMaskEnableToggled));
    ToolParamBlock *area = Gtk::manage(new ToolParamBlock());
    hb = Gtk::manage(new Gtk::HBox());
    areaMaskButtonsHb = hb;

    areaMaskInverted = Gtk::manage(new Gtk::CheckButton(M("TP_LABMASKS_AREA_INVERTED")));
    areaMaskInverted->signal_toggled().connect(sigc::mem_fun(*this, &LabMasksPanel::onAreaMaskInvertedChanged));
    hb->pack_start(*areaMaskInverted, Gtk::PACK_EXPAND_WIDGET);

    areaMaskCopy = Gtk::manage(new Gtk::Button());
    areaMaskCopy->add(*Gtk::manage(new RTImage("copy.png")));
    areaMaskCopy->set_tooltip_text(M("TP_LABMASKS_AREA_MASK_COPY_TOOLTIP"));
    areaMaskCopy->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onAreaMaskCopyPressed));
    add_button(areaMaskCopy, hb, 24);
    
    areaMaskPaste = Gtk::manage(new Gtk::Button());
    areaMaskPaste->add(*Gtk::manage(new RTImage("paste.png")));
    areaMaskPaste->set_tooltip_text(M("TP_LABMASKS_AREA_MASK_PASTE_TOOLTIP"));
    areaMaskPaste->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onAreaMaskPastePressed));
    add_button(areaMaskPaste, hb, 24);

    areaMaskToggle = new Gtk::ToggleButton();
    areaMaskToggle->get_style_context()->add_class("independent");
    areaMaskToggle->add(*Gtk::manage(new RTImage("crosshair-adjust.png")));
    areaMaskToggle->set_tooltip_text(M("EDIT_OBJECT_TOOLTIP"));
    areaMaskToggle->signal_toggled().connect(sigc::mem_fun(*this, &LabMasksPanel::onAreaMaskToggleChanged));
    hb->pack_start(*areaMaskToggle, Gtk::PACK_SHRINK, 0);
    area->pack_start(*hb);

    const auto add_adjuster =
        [&](Adjuster *a, Gtk::Box *b) -> void
        {
            areaMaskAdjusters.push_back(a);
            a->setAdjusterListener(this);
            if (b) {
                b->pack_start(*a);
            }
        };
    
    areaMaskFeather = Gtk::manage(new Adjuster(M("TP_LABMASKS_AREA_FEATHER"), 0, 100, 0.1, 0));
    add_adjuster(areaMaskFeather, area);

    CurveEditorGroup *cg = Gtk::manage(new CurveEditorGroup(options.lastToneCurvesDir, M("TP_LABMASKS_AREA_CONTRAST")));
    cg->setCurveListener(this);

    areaMaskContrast = static_cast<DiagonalCurveEditor *>(cg->addCurve(CT_Diagonal, ""));
    cg->curveListComplete();
    if (area) {
        area->pack_start(*cg, Gtk::PACK_SHRINK, 2);
    }
    
    // areaMaskContrast = Gtk::manage(new Adjuster(M("TP_EXPOSURE_CONTRAST"), 0, 100, 1, 0));
    // add_adjuster(areaMaskContrast, area);

    Gtk::Frame *areaFrame = Gtk::manage(new Gtk::Frame(M("TP_LABMASKS_AREA_SHAPES")));
    areaFrame->set_label_align(0.025, 0.5);
    areaFrame->set_border_width(4);

    areaMaskShapes = Gtk::manage(new Gtk::ListViewText(2));
    areaMaskShapes->set_size_request(-1, 80);
    areaMaskShapes->set_can_focus(false);
    areaMaskShapes->set_column_title(0, "#");
    areaMaskShapes->set_column_title(1, M("TP_LABMASKS_AREA_SHAPE"));
    areaMaskShapes->set_activate_on_single_click(true);
    shapeSelectionConn = areaMaskShapes->get_selection()->signal_changed().connect(sigc::mem_fun(this, &LabMasksPanel::onAreaShapeSelectionChanged));
    hb = Gtk::manage(new Gtk::HBox());
    scroll = Gtk::manage(new Gtk::ScrolledWindow());
    scroll->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_NEVER);
    scroll->add(*areaMaskShapes);
    hb->pack_start(*scroll, Gtk::PACK_EXPAND_WIDGET);
    vb = Gtk::manage(new Gtk::VBox());
    areaMaskReset = Gtk::manage(new Gtk::Button());
    areaMaskReset->add(*Gtk::manage(new RTImage("undo-small.png")));
    areaMaskReset->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onAreaShapeResetPressed));
    add_button(areaMaskReset, vb);
    areaMaskAdd = Gtk::manage(new Gtk::Button());
    areaMaskAdd->add(*Gtk::manage(new RTImage("add-small.png")));
    areaMaskAdd->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onAreaShapeAddPressed));
    add_button(areaMaskAdd, vb);
    areaMaskRemove = Gtk::manage(new Gtk::Button());
    areaMaskRemove->add(*Gtk::manage(new RTImage("remove-small.png")));
    areaMaskRemove->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onAreaShapeRemovePressed));
    add_button(areaMaskRemove, vb);
    hb->pack_start(*vb, Gtk::PACK_SHRINK);
    
    vb = Gtk::manage(new Gtk::VBox());
    vb->set_spacing(2);
    vb->pack_start(*hb);
    
    areaMaskX = Gtk::manage(new Adjuster(M("TP_LABMASKS_AREA_X"), -100, 100, 0.1, 0));
    add_adjuster(areaMaskX, vb);
    areaMaskY = Gtk::manage(new Adjuster(M("TP_LABMASKS_AREA_Y"), -100, 100, 0.1, 0));
    add_adjuster(areaMaskY, vb);
    areaMaskWidth = Gtk::manage(new Adjuster(M("TP_LABMASKS_AREA_WIDTH"), 1, 200, 0.1, 100));
    areaMaskWidth->setLogScale(10, 1);
    add_adjuster(areaMaskWidth, vb);
    areaMaskHeight = Gtk::manage(new Adjuster(M("TP_LABMASKS_AREA_HEIGHT"), 1, 200, 0.1, 100));
    areaMaskHeight->setLogScale(10, 1);
    add_adjuster(areaMaskHeight, vb);
    areaMaskAngle = Gtk::manage(new Adjuster(M("TP_LABMASKS_AREA_ANGLE"), 0, 180, 0.1, 0));
    add_adjuster(areaMaskAngle, vb);
    areaMaskRoundness = Gtk::manage(new Adjuster(M("TP_LABMASKS_AREA_ROUNDNESS"), 0, 100, 0.1, 0));
    add_adjuster(areaMaskRoundness, vb);

    areaFrame->add(*vb);
    vb = Gtk::manage(new Gtk::VBox());
    vb->pack_start(*area);
    vb->pack_start(*areaFrame);
    areaMask->add(*vb, false);
    areaMask->setLevel(2);
    pack_start(*areaMask);

    showMask = Gtk::manage(new Gtk::CheckButton(M("TP_LABMASKS_SHOW")));
    showMask->signal_toggled().connect(sigc::mem_fun(*this, &LabMasksPanel::onShowMaskChanged));
    pack_start(*showMask, Gtk::PACK_SHRINK, 4);

    maskBlur->delay = options.adjusterMaxDelay;
}


LabMasksPanel::~LabMasksPanel()
{
    delete areaMaskToggle;
}


ToolPanelListener *LabMasksPanel::getListener()
{
    if (listenerDisabled.empty() || !listenerDisabled.back()) {
        return cp_->listener();
    }
    return nullptr;
}


void LabMasksPanel::disableListener()
{
    listenerDisabled.push_back(true);
}


void LabMasksPanel::enableListener()
{
    listenerDisabled.pop_back();
}


void LabMasksPanel::onSelectionChanged()
{
    auto s = list->get_selected();
    if (!s.empty()) {
        int idx = s[0];
        // update the selected values
        maskGet(selected_);
        cp_->selectionChanging(selected_);
        selected_ = idx;
        area_shape_index_ = 0;
        maskShow(selected_);
        if (showMask->get_active()) {
            onShowMaskChanged();
        }
    }
}


void LabMasksPanel::maskGet(int idx)
{
    if (idx < 0 || size_t(idx) >= masks_.size()) {
        return;
    }
    
    auto &r = masks_[idx];
    r.hueMask = hueMask->getCurve();
    r.chromaticityMask = chromaticityMask->getCurve();
    r.lightnessMask = lightnessMask->getCurve();
    r.maskBlur = maskBlur->getValue();
    r.areaEnabled = areaMask->getEnabled();
    r.areaMask.inverted = areaMaskInverted->get_active();
    r.areaMask.feather = areaMaskFeather->getValue();
    r.areaMask.contrast = areaMaskContrast->getCurve();
    if (area_shape_index_ < r.areaMask.shapes.size()) {
        auto &a = r.areaMask.shapes[area_shape_index_];
        a.x = areaMaskX->getValue();
        a.y = areaMaskY->getValue();
        a.width = areaMaskWidth->getValue();
        a.height = areaMaskHeight->getValue();
        a.angle = areaMaskAngle->getValue();
        a.roundness = areaMaskRoundness->getValue();
    }
}


void LabMasksPanel::onAddPressed()
{
    if (!cp_->addPressed()) {
        return;
    }

    listEdited = true;
    selected_ = masks_.size();
    masks_.push_back(rtengine::procparams::LabCorrectionMask());
    masks_.back().areaMask.shapes = {defaultAreaShape};
    populateList();
    area_shape_index_ = 0;
    maskShow(selected_);

    auto l = getListener();
    if (l) {
        l->panelChanged(EvMaskList, M("HISTORY_CHANGED"));
    }
}


void LabMasksPanel::onRemovePressed()
{
    if (list->size() <= 1 || !cp_->removePressed(selected_)) {
        return;
    }
    
    listEdited = true;
    masks_.erase(masks_.begin() + selected_);
    selected_ = rtengine::LIM(selected_-1, 0u, unsigned(masks_.size()-1));
    populateList();
    area_shape_index_ = 0;
    maskShow(selected_);

    auto l = getListener();
    if (l) {
        l->panelChanged(EvMaskList, M("HISTORY_CHANGED"));
    }
}


void LabMasksPanel::onUpPressed()
{
    if (selected_ > 0 && cp_->moveUpPressed(selected_)) {
        listEdited = true;
        
        auto r = masks_[selected_];
        masks_.erase(masks_.begin() + selected_);
        --selected_;
        masks_.insert(masks_.begin() + selected_, r);
        populateList();

        auto l = getListener();
        if (l) {
            l->panelChanged(EvMaskList, M("HISTORY_CHANGED"));
        }
    }
}


void LabMasksPanel::onDownPressed()
{
    if (selected_ < masks_.size()-1 && cp_->moveDownPressed(selected_)) {
        listEdited = true;
        
        auto r = masks_[selected_];
        masks_.erase(masks_.begin() + selected_);
        ++selected_;
        masks_.insert(masks_.begin() + selected_, r);
        populateList();

        auto l = getListener();
        if (l) {
            l->panelChanged(EvMaskList, M("HISTORY_CHANGED"));
        }
    }
}


void LabMasksPanel::onCopyPressed()
{
    if (selected_ < masks_.size() && cp_->copyPressed(selected_)) {
        listEdited = true;
        
        auto r = masks_[selected_];
        masks_.push_back(r);
        selected_ = masks_.size()-1;
        populateList();

        auto l = getListener();
        if (l) {
            l->panelChanged(EvMaskList, M("HISTORY_CHANGED"));
        }
    }
}


void LabMasksPanel::onResetPressed()
{
    if (cp_->resetPressed()) {
        listEdited = true;
        auto l = getListener();
        if (l) {
            l->panelChanged(EvMaskList, M("HISTORY_CHANGED"));
        }
    }
}


void LabMasksPanel::onShowMaskChanged()
{
    auto l = getListener();
    if (l) {
        l->panelChanged(EvShowMask, showMask->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}


void LabMasksPanel::populateList()
{
    ConnectionBlocker b(selectionConn);
    list->clear_items();
    rtengine::procparams::LabCorrectionMask dflt;

    int n = cp_->getColumnCount();
    for (size_t i = 0; i < masks_.size(); ++i) {
        auto &r = masks_[i];
        auto j = list->append(std::to_string(i+1));
        for (int c = 0; c < n; ++c) {
            list->set_text(j, c+1, cp_->getColumnContent(c, j));
        }
        Glib::ustring am("");
        if (!r.areaMask.isTrivial()) {
            am = Glib::ustring::compose("\n%1 shape%2", r.areaMask.shapes.size(), r.areaMask.shapes.size() > 1 ? "s" : "");
        }
        list->set_text(
            j, n+1, Glib::ustring::compose(
                "%1%2%3%4%5",
                hasMask(dflt.hueMask, r.hueMask) ? "H" : "",
                hasMask(dflt.chromaticityMask, r.chromaticityMask) ? "C" : "",
                hasMask(dflt.lightnessMask, r.lightnessMask) ? "L" : "",
                r.maskBlur ? Glib::ustring::compose(" b=%1", r.maskBlur) : "",
                am));
    }
}


void LabMasksPanel::maskShow(int idx, bool list_only, bool unsub)
{
    disableListener();
    rtengine::procparams::LabCorrectionMask dflt;
    auto &r = masks_[idx];
    if (!list_only) {
        cp_->selectionChanged(idx);
        hueMask->setCurve(r.hueMask);
        chromaticityMask->setCurve(r.chromaticityMask);
        lightnessMask->setCurve(r.lightnessMask);
        maskBlur->setValue(r.maskBlur);

        if (unsub && isCurrentSubscriber()) {
            unsubscribe();
        }
        areaMaskToggle->set_active(false);
        areaMask->setEnabled(r.areaEnabled);
        areaMaskInverted->set_active(r.areaMask.inverted);
        areaMaskFeather->setValue(r.areaMask.feather);
        areaMaskContrast->setCurve(r.areaMask.contrast);
        if (area_shape_index_ < r.areaMask.shapes.size()) {
            auto &a = r.areaMask.shapes[area_shape_index_];
            areaMaskX->setValue(a.x);
            areaMaskY->setValue(a.y);
            areaMaskWidth->setValue(a.width);
            areaMaskHeight->setValue(a.height);
            areaMaskAngle->setValue(a.angle);
            areaMaskRoundness->setValue(a.roundness);
        }
        populateShapeList(idx, area_shape_index_);
        updateAreaMask(false);
    }

    int n = cp_->getColumnCount();
    for (int c = 0; c < n; ++c) {
        list->set_text(idx, c+1, cp_->getColumnContent(c, idx));
    }
    Glib::ustring am("");
    if (r.areaEnabled && !r.areaMask.isTrivial()) {
        am = Glib::ustring::compose("\n%1 shape%2", r.areaMask.shapes.size(),
                                    r.areaMask.shapes.size() > 1 ? "s" : "");
    }
    list->set_text(
        idx, n+1, Glib::ustring::compose(
            "%1%2%3%4%5",
            hasMask(dflt.hueMask, r.hueMask) ? "H" : "",
            hasMask(dflt.chromaticityMask, r.chromaticityMask) ? "C" : "",
            hasMask(dflt.lightnessMask, r.lightnessMask) ? "L" : "",
            r.maskBlur ? Glib::ustring::compose(" b=%1", r.maskBlur) : "", am));
    Gtk::TreePath pth;
    pth.push_back(idx);
    list->get_selection()->select(pth);
    enableListener();
}


void LabMasksPanel::setEditProvider(EditDataProvider *provider)
{
    hueMask->setEditProvider(provider);
    chromaticityMask->setEditProvider(provider);
    lightnessMask->setEditProvider(provider);
    AreaMask::setEditProvider(provider);
}


void LabMasksPanel::onAreaMaskToggleChanged()
{
    if (areaMaskToggle->get_active()) {
        subscribe();
    } else {
        unsubscribe();
    }
}


void LabMasksPanel::onAreaMaskInvertedChanged()
{
    auto l = getListener();
    if (l) {
        l->panelChanged(EvAreaMask, M("GENERAL_CHANGED"));
    }
}


void LabMasksPanel::updateAreaMask(bool from_mask)
{
    disableListener();
    if (from_mask) {
        areaMaskX->setValue(center_x_);
        areaMaskY->setValue(center_y_);
        areaMaskWidth->setValue(width_);
        areaMaskHeight->setValue(height_);
        areaMaskAngle->setValue(angle_);
    } else {
        center_x_ = areaMaskX->getValue();
        center_y_ = areaMaskY->getValue();
        width_ = areaMaskWidth->getValue();
        height_ = areaMaskHeight->getValue();
        angle_ = areaMaskAngle->getValue();
        updateGeometry();
    }
    enableListener();
}


bool LabMasksPanel::button1Released()
{
    if (last_object_ != -1) {
        updateAreaMask(true);
        onAreaShapeSelectionChanged();
        populateShapeList(selected_, area_shape_index_);
        auto l = getListener();
        if (l) {
            l->panelChanged(EvAreaMask, M("GENERAL_CHANGED"));
        }
    }
    return AreaMask::button1Released();
}


void LabMasksPanel::switchOffEditMode()
{
    areaMaskToggle->set_active(false);
    AreaMask::switchOffEditMode();
}


void LabMasksPanel::onAreaMaskEnableToggled()
{
    auto l = getListener();
    if (l) {
        l->panelChanged(EvAreaMask, areaMask->getEnabled() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}


void LabMasksPanel::updateAreaMaskDefaults(const rtengine::procparams::ProcParams *params)
{
    EditDataProvider *provider = getEditProvider();

    if (!provider) {
        return;
    }

    int imW = 0, imH = 0;
    provider->getImageSize(imW, imH);
    if (!imW || !imH) {
        return;
    }

    defaultAreaShape = rtengine::procparams::AreaMask::Shape();
    if (!params->crop.enabled) {
        return;
    }

    defaultAreaShape.width = double(params->crop.w)/double(imW) * 100.0;
    defaultAreaShape.height = double(params->crop.h)/double(imH) * 100.0;
    defaultAreaShape.x = (double(params->crop.x + params->crop.w * 0.5) / (imW * 0.5) - 1) * 100.0;
    defaultAreaShape.y = (double(params->crop.y + params->crop.h * 0.5) / (imH * 0.5) - 1) * 100.0;

    if (masks_.size() == 1 && !masks_[0].areaEnabled && masks_[0].areaMask.isTrivial()) {
        masks_[0].areaMask.shapes = {defaultAreaShape};
        areaShapeSelect(0, false);
        populateShapeList(selected_, area_shape_index_);
    }
}


void LabMasksPanel::curveChanged(CurveEditor* ce)
{
    auto l = getListener();
    if (l) {
        if (ce == areaMaskContrast) {
            l->panelChanged(EvAreaMask, M("GENERAL_CHANGED"));
        } else if (ce == hueMask) {
            l->panelChanged(EvHMask, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == chromaticityMask) {
            l->panelChanged(EvCMask, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == lightnessMask) {
            l->panelChanged(EvLMask, M("HISTORY_CUSTOMCURVE"));
        }
    }
}


void LabMasksPanel::adjusterChanged(Adjuster *a, double newval)
{
    auto l = getListener();
    if (!l) {
        return;
    }

    if (a == maskBlur) {
        l->panelChanged(EvMaskBlur, a->getTextValue());
    } else if (std::find(areaMaskAdjusters.begin(), areaMaskAdjusters.end(), a) != areaMaskAdjusters.end()) {
        onAreaShapeSelectionChanged();
        populateShapeList(selected_, area_shape_index_);
        l->panelChanged(EvAreaMask, M("GENERAL_CHANGED"));
    }
    maskShow(selected_, true);
}


void LabMasksPanel::adjusterAutoToggled(Adjuster *a, bool newval)
{
}


void LabMasksPanel::setMasks(const std::vector<rtengine::procparams::LabCorrectionMask> &masks, int show_mask_idx)
{
    disableListener();
    ConnectionBlocker b(selectionConn);
    
    masks_ = masks;
    selected_ = 0;
    if (show_mask_idx >= 0) {
        selected_ = show_mask_idx;
        showMask->set_active(true);
    } else {
        showMask->set_active(false);
    }
    populateList();
    area_shape_index_ = 0;
    maskShow(selected_);
    enableListener();
}
        

void LabMasksPanel::getMasks(std::vector<rtengine::procparams::LabCorrectionMask> &masks, int &show_mask_idx)
{
    maskGet(selected_);
    masks = masks_;
    if (showMask->get_active()) {
        show_mask_idx = selected_;
    } else {
        show_mask_idx = -1;
    }
}


int LabMasksPanel::getSelected()
{
    return selected_;
}


void LabMasksPanel::setBatchMode()
{
    removeIfThere(areaMaskButtonsHb, areaMaskToggle, false);
}


void LabMasksPanel::colorForValue(double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller *caller)
{
    float R = 0.f, G = 0.f, B = 0.f;

    if (callerId == ID_HUE_MASK) {
        float x = valX - 1.f/6.f;
        if (x < 0.f) {
            x += 1.f;
        }
        x = rtengine::log2lin(x, 3.f);
        rtengine::Color::hsv2rgb01(x, 0.5f, 0.5f, R, G, B);        
    } else if (callerId == ID_HUE_MASK+1) {
        rtengine::Color::hsv2rgb01(float(valY), float(valX), 0.5f, R, G, B);
    }

    caller->ccRed = double(R);
    caller->ccGreen = double(G);
    caller->ccBlue = double(B);
}


float LabMasksPanel::blendPipetteValues(CurveEditor *ce, float chan1, float chan2, float chan3)
{
    if (ce == chromaticityMask && chan1 > 0.f) {
        return rtengine::lin2log(chan1, 10.f);
    } else if (ce == hueMask && chan1 > 0.f) {
        float x = chan1 + 1.f/6.f;
        if (x > 1.f) {
            x -= 1.f;
        }
        return rtengine::lin2log(x, 3.f);
    }
    return CurveListener::blendPipetteValues(ce, chan1, chan2, chan3);
}


void LabMasksPanel::setEdited(bool yes)
{
    listEdited = yes;
    hueMask->setUnChanged(!yes);
    chromaticityMask->setUnChanged(!yes);
    lightnessMask->setUnChanged(!yes);
    maskBlur->setEditedState(yes ? Edited : UnEdited);
    showMask->set_inconsistent(!yes);
    areaMask->set_inconsistent(!yes);
    areaMaskInverted->set_inconsistent(!yes);
    for (auto a : areaMaskAdjusters) {
        a->setEditedState(yes ? Edited : UnEdited);
    }
    areaMaskContrast->setUnChanged(!yes);
}


bool LabMasksPanel::getEdited()
{
    for (auto a : areaMaskAdjusters) {
        if (a->getEditedState() == Edited) {
            return true;
        }
    }
    return listEdited
        || !hueMask->isUnChanged()
        || !chromaticityMask->isUnChanged()
        || !lightnessMask->isUnChanged()
        || maskBlur->getEditedState() == Edited
        || !showMask->get_inconsistent()
        || !areaMask->get_inconsistent()
        || !areaMaskInverted->get_inconsistent()
        || !areaMaskContrast->isUnChanged();
}


void LabMasksPanel::updateSelected()
{
    maskShow(selected_, true, false);
}


void LabMasksPanel::onAreaShapeSelectionChanged()
{
    if (selected_ < masks_.size() && area_shape_index_ < masks_[selected_].areaMask.shapes.size()) {
        disableListener();
        
        auto &s = masks_[selected_].areaMask.shapes[area_shape_index_];
        updateAreaMask(false);
        s.x = center_x_;
        s.y = center_y_;
        s.width = width_;
        s.height = height_;
        s.angle = angle_;
        s.roundness = areaMaskRoundness->getValue();

        auto sel = areaMaskShapes->get_selected();
        unsigned int newidx = sel.empty() ? area_shape_index_ : sel[0];

        if (newidx != area_shape_index_) {
            areaShapeSelect(newidx, false);
        }
        enableListener();
    }
}


void LabMasksPanel::onAreaShapeResetPressed()
{
    if (selected_ < masks_.size()) {
        listEdited = true;
        disableListener();
        masks_[selected_].areaMask.shapes = {defaultAreaShape};
        area_shape_index_ = 0;
        center_x_ = defaultAreaShape.x;
        center_y_ = defaultAreaShape.y;
        width_ = defaultAreaShape.width;
        height_ = defaultAreaShape.height;
        angle_ = defaultAreaShape.angle;
        updateGeometry();
        updateAreaMask(true);
        populateShapeList(selected_, area_shape_index_);
        maskShow(selected_, true);        
        enableListener();
        auto l = getListener();
        if (l && areaMask->getEnabled()) {
            l->panelChanged(EvAreaMask, M("GENERAL_CHANGED"));
        }
    }
}


void LabMasksPanel::onAreaShapeAddPressed()
{
    if (selected_ < masks_.size()) {
        listEdited = true;
        auto &am = masks_[selected_].areaMask;
        auto s = am.shapes.empty() ? defaultAreaShape : am.shapes.back();
        am.shapes.push_back(s);
        populateShapeList(selected_, -1);
        areaShapeSelect(am.shapes.size()-1, true);
        maskShow(selected_, true);        
        auto l = getListener();
        if (l && areaMask->getEnabled()) {
            l->panelChanged(EvAreaMask, M("GENERAL_CHANGED"));
        }
    }
}


void LabMasksPanel::onAreaShapeRemovePressed()
{
    if (selected_ < masks_.size() && area_shape_index_ < masks_[selected_].areaMask.shapes.size() && masks_[selected_].areaMask.shapes.size() > 1) {
        listEdited = true;
        masks_[selected_].areaMask.shapes.erase(masks_[selected_].areaMask.shapes.begin() + area_shape_index_);
        populateShapeList(selected_, -1);
        areaShapeSelect(area_shape_index_ > 0 ? area_shape_index_ - 1 : 0, true);
        maskShow(selected_, true);        
        auto l = getListener();
        if (l && areaMask->getEnabled()) {
            l->panelChanged(EvAreaMask, M("GENERAL_CHANGED"));
        }
    }
}


void LabMasksPanel::populateShapeList(int idx, int sel)
{
    ConnectionBlocker b(shapeSelectionConn);
    areaMaskShapes->clear_items();
    auto &r = masks_[idx];
    const auto rd =
        [](double v) -> double
        {
            return int(v * 10) / 10.0;
        };
    
    for (size_t i = 0; i < r.areaMask.shapes.size(); ++i) {
        auto &a = r.areaMask.shapes[i];
        auto j = areaMaskShapes->append(std::to_string(i+1));
        areaMaskShapes->set_text(j, 1, Glib::ustring::compose("%1 %2 %3 %4 %5 %6", rd(a.x), rd(a.y), rd(a.width), rd(a.height), rd(a.angle), rd(a.roundness)));
    }
    if (sel >= 0) {
        Gtk::TreePath pth;
        pth.push_back(sel);
        areaMaskShapes->get_selection()->select(pth);
    }
}


void LabMasksPanel::onAreaMaskCopyPressed()
{
    if (selected_ < masks_.size()) {
        clipboard.setAreaMask(masks_[selected_].areaMask);
    }
}


void LabMasksPanel::onAreaMaskPastePressed()
{
    if (selected_ < masks_.size() && clipboard.hasAreaMask()) {
        listEdited = true;
        disableListener();
        masks_[selected_].areaMask = clipboard.getAreaMask();
        area_shape_index_ = 0;
        if (area_shape_index_ < masks_[selected_].areaMask.shapes.size()) {
            auto &a = masks_[selected_].areaMask;
            auto &s = a.shapes[0];
            center_x_ = s.x;
            center_y_ = s.y;
            width_ = s.width;
            height_ = s.height;
            angle_ = s.angle;
            updateGeometry();
            updateAreaMask(true);
            areaMaskFeather->setValue(a.feather);
            areaMaskContrast->setCurve(a.contrast);
            areaMaskRoundness->setValue(s.roundness);
            areaMaskInverted->set_active(a.inverted);
        }
        populateShapeList(selected_, area_shape_index_);
        maskShow(selected_, true);        
        enableListener();
        auto l = getListener();
        if (l && areaMask->getEnabled()) {
            l->panelChanged(EvAreaMask, M("GENERAL_CHANGED"));
        }
    }
}


void LabMasksPanel::areaShapeSelect(int sel, bool update_list)
{
    area_shape_index_ = sel;
    auto &ns = masks_[selected_].areaMask.shapes[sel];
    center_x_ = ns.x;
    center_y_ = ns.y;
    width_ = ns.width;
    height_ = ns.height;
    angle_ = ns.angle;
    updateGeometry();
    updateAreaMask(true);
    areaMaskRoundness->setValue(ns.roundness);
    if (areaMaskToggle->get_active()) {
        areaMaskToggle->set_active(false);
        areaMaskToggle->set_active(true);
    }

    if (update_list) {
        ConnectionBlocker b(shapeSelectionConn);
        Gtk::TreePath pth;
        pth.push_back(sel);
        areaMaskShapes->get_selection()->select(pth);
    }
}
