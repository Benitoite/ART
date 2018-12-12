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

LabMasksPanel::LabMasksPanel(Gtk::Widget *child, rtengine::ProcEvent h_mask, rtengine::ProcEvent c_mask, rtengine::ProcEvent l_mask, rtengine::ProcEvent blur, rtengine::ProcEvent area_mask):
    Gtk::VBox(),
    masks_(),
    selected_(0),
    EvHMask(h_mask),
    EvCMask(c_mask),
    EvLMask(l_mask),
    EvMaskBlur(blur),
    EvAreaMask(area_mask)
{
    const auto add_button =
        [](Gtk::Button *btn, Gtk::Box *box) -> void
        {
            setExpandAlignProperties(btn, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_START);
            btn->set_relief(Gtk::RELIEF_NONE);
            btn->get_style_context()->add_class(GTK_STYLE_CLASS_FLAT);
            btn->set_can_focus(false);
            btn->set_size_request(-1, 20);
            box->pack_start(*btn, false, false);
        };
    
    int n = getColumnCount();
    list = Gtk::manage(new Gtk::ListViewText(n+2));
    list->set_size_request(-1, 150);
    list->set_can_focus(false);
    list->set_column_title(0, "#");
    for (int i = 0, i < n; ++i) {
        list->set_column_title(i+1, getColumnHeader(i));
    }
    list->set_column_title(n+1, M("TP_LABMASKS_MASK_COLUMN"));
    list->set_activate_on_single_click(true);
    list->get_selection()->signal_changed().connect(sigc::mem_fun(this, &LabMasksPanel::onSelectionChanged));
    Gtk::HBox *hb = Gtk::manage(new Gtk::HBox());
    Gtk::ScrolledWindow *scroll = Gtk::manage(new Gtk::ScrolledWindow());
    scroll->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_NEVER);
    scroll->add(*list);
    hb->pack_start(*scroll, Gtk::PACK_EXPAND_WIDGET);
    Gtk::VBox *vb = Gtk::manage(new Gtk::VBox());
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

    pack_start(*child);

    pack_start(*Gtk::manage(new Gtk::HSeparator()));

    CurveEditorGroup *EditorG = Gtk::manage(new CurveEditorGroup("", M("TP_COLORTONING_LABREGION_MASK"), 0.7));
    EditorG->setCurveListener(this);

    rtengine::LabCorrectionMask default_params;

    hueMask = static_cast<FlatCurveEditor *>(EditorG->addCurve(CT_Flat, M("TP_COLORTONING_LABREGION_HUEMASK"), nullptr, false, true));
    hueMask->setIdentityValue(0.);
    hueMask->setResetCurve(FlatCurveType(default_params.hueMask[0]), default_params.hueMask);
    hueMask->setCurveColorProvider(this, ID_LABREGION_HUE);
    hueMask->setBottomBarColorProvider(this, ID_LABREGION_HUE);
    hueMask->setEditID(EUID_Lab_HHCurve, BT_SINGLEPLANE_FLOAT);

    chromaticityMask = static_cast<FlatCurveEditor *>(EditorG->addCurve(CT_Flat, M("TP_COLORTONING_LABREGION_CHROMATICITYMASK"), nullptr, false, false));
    chromaticityMask->setIdentityValue(0.);
    chromaticityMask->setResetCurve(FlatCurveType(default_params.chromaticityMask[0]), default_params.chromaticityMask);
    chromaticityMask->setBottomBarColorProvider(this, ID_LABREGION_HUE+1);
    chromaticityMask->setEditID(EUID_Lab_CCurve, BT_SINGLEPLANE_FLOAT);

    lightnessMask = static_cast<FlatCurveEditor *>(EditorG->addCurve(CT_Flat, M("TP_COLORTONING_LABREGION_LIGHTNESSMASK"), nullptr, false, false));
    lightnessMask->setIdentityValue(0.);
    lightnessMask->setResetCurve(FlatCurveType(default_params.lightnessMask[0]), default_params.lightnessMask);
    lightnessMask->setBottomBarBgGradient(milestones);
    lightnessMask->setEditID(EUID_Lab_LCurve, BT_SINGLEPLANE_FLOAT);

    EditorG->curveListComplete();
    EditorG->show();
    pack_start(*EditorG, Gtk::PACK_SHRINK, 2);

    maskBlur = Gtk::manage(new Adjuster(M("TP_COLORTONING_LABREGION_MASKBLUR"), -10, 100, 0.1, 0));
    maskBlur->setLogScale(10, 0);
    maskBlur->setAdjusterListener(this);
    pack_start(*maskBlur);

    // TODO

    labAreaMask = Gtk::manage(new MyExpander(true, M("TP_COLORTONING_LABREGION_AREAMASK")));
    labAreaMask->signal_enabled_toggled().connect(sigc::mem_fun(*this, &LabMasksPanel::labAreaMaskEnableToggled));
    ToolParamBlock *area = Gtk::manage(new ToolParamBlock());
    hb = Gtk::manage(new Gtk::HBox());

    labAreaMaskInverted = Gtk::manage(new Gtk::CheckButton(M("TP_COLORTONING_LABREGION_AREAMASK_INVERTED")));
    labAreaMaskInverted->signal_toggled().connect(sigc::mem_fun(*this, &LabMasksPanel::labAreaMaskInvertedChanged));
    hb->pack_start(*labAreaMaskInverted, Gtk::PACK_EXPAND_WIDGET);

    labAreaMaskToggle = Gtk::manage(new Gtk::ToggleButton());
    labAreaMaskToggle->get_style_context()->add_class("independent");
    labAreaMaskToggle->add(*Gtk::manage(new RTImage("crosshair-adjust.png")));
    labAreaMaskToggle->set_tooltip_text(M("EDIT_OBJECT_TOOLTIP"));
    labAreaMaskToggle->signal_toggled().connect(sigc::mem_fun(*this, &LabMasksPanel::labAreaMaskToggleChanged));
    hb->pack_start(*labAreaMaskToggle, Gtk::PACK_SHRINK, 0);
    area->pack_start(*hb);

    labAreaMaskX = Gtk::manage(new Adjuster(M("TP_COLORTONING_LABREGION_AREAMASK_X"), -100, 100, 0.1, 0));
    labAreaMaskAdjusters.push_back(labAreaMaskX);
    labAreaMaskY = Gtk::manage(new Adjuster(M("TP_COLORTONING_LABREGION_AREAMASK_Y"), -100, 100, 0.1, 0));
    labAreaMaskAdjusters.push_back(labAreaMaskY);
    labAreaMaskWidth = Gtk::manage(new Adjuster(M("TP_COLORTONING_LABREGION_AREAMASK_WIDTH"), 1, 200, 0.1, 100));
    labAreaMaskWidth->setLogScale(10, 1);
    labAreaMaskAdjusters.push_back(labAreaMaskWidth);
    labAreaMaskHeight = Gtk::manage(new Adjuster(M("TP_COLORTONING_LABREGION_AREAMASK_HEIGHT"), 1, 200, 0.1, 100));
    labAreaMaskHeight->setLogScale(10, 1);
    labAreaMaskAdjusters.push_back(labAreaMaskHeight);
    labAreaMaskAngle = Gtk::manage(new Adjuster(M("TP_COLORTONING_LABREGION_AREAMASK_ANGLE"), 0, 180, 0.1, 0));
    labAreaMaskAdjusters.push_back(labAreaMaskAngle);
    labAreaMaskFeather = Gtk::manage(new Adjuster(M("TP_COLORTONING_LABREGION_AREAMASK_FEATHER"), 0, 100, 0.1, 0));
    labAreaMaskAdjusters.push_back(labAreaMaskFeather);
    labAreaMaskRoundness = Gtk::manage(new Adjuster(M("TP_COLORTONING_LABREGION_AREAMASK_ROUNDNESS"), 0, 100, 0.1, 0));
    labAreaMaskAdjusters.push_back(labAreaMaskRoundness);

    for (auto a : labAreaMaskAdjusters) {
        a->setAdjusterListener(this);
        area->pack_start(*a);
    }

    labAreaMask->add(*area, false);
    labAreaMask->setLevel(2);
    pack_start(*labAreaMask);

    ShowMask = Gtk::manage(new Gtk::CheckButton(M("TP_COLORTONING_LABREGION_SHOWMASK")));
    ShowMask->signal_toggled().connect(sigc::mem_fun(*this, &LabMasksPanel::ShowMaskChanged));
    pack_start(*ShowMask, Gtk::PACK_SHRINK, 4);

    maskBlur->delay = options.adjusterMaxDelay;
}


LabMasksPanel::~LabMasksPanel()
{
}
