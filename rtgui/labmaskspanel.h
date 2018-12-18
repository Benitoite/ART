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

#pragma once

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "areamask.h"
#include "colorprovider.h"
#include "curvelistener.h"
#include "curveeditorgroup.h"
#include "curveeditor.h"


class LabMasksContentProvider {
public:
    virtual ~LabMasksContentProvider() {}

    virtual Gtk::Widget *getWidget() = 0;
    virtual void getEvents(rtengine::ProcEvent &mask_list, rtengine::ProcEvent &h_mask, rtengine::ProcEvent &c_mask, rtengine::ProcEvent &l_mask, rtengine::ProcEvent &blur, rtengine::ProcEvent &show, rtengine::ProcEvent &area_mask) = 0;

    virtual ToolPanelListener *listener() = 0;

    virtual void selectionChanging(int idx) = 0;
    virtual void selectionChanged(int idx) = 0;
    virtual bool addPressed() = 0;
    virtual bool removePressed(int idx) = 0;
    virtual bool copyPressed(int idx) = 0;
    virtual bool moveUpPressed(int idx) = 0;
    virtual bool moveDownPressed(int idx) = 0;
    virtual bool resetPressed() = 0;

    virtual int getColumnCount() = 0;
    virtual Glib::ustring getColumnHeader(int col) = 0;
    virtual Glib::ustring getColumnContent(int col, int row) = 0;

    virtual void getEditIDs(EditUniqueID &hcurve, EditUniqueID &ccurve, EditUniqueID &lcurve) = 0;
};


class LabMasksPanel:
    public Gtk::VBox,
    public AdjusterListener,
    public CurveListener,
    public ColorProvider,
    public AreaMask {
public:
    LabMasksPanel(LabMasksContentProvider *cp);
    ~LabMasksPanel();

    void setMasks(const std::vector<rtengine::LabCorrectionMask> &masks, int show_mask_idx);
    void getMasks(std::vector<rtengine::LabCorrectionMask> &masks, int &show_mask_idx);
    int getSelected();

    void adjusterChanged(Adjuster *a, double newval) override;
    void adjusterAutoToggled(Adjuster *a, bool newval) override;
    void curveChanged(CurveEditor* ce) override;
    float blendPipetteValues(CurveEditor *ce, float chan1, float chan2, float chan3) override;

    void setBatchMode();

    bool button1Released() override;
    void switchOffEditMode() override;
    void setEditProvider(EditDataProvider *provider);
    void colorForValue(double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller *caller) override;

    void updateAreaMaskDefaults(const rtengine::procparams::ProcParams* params);

    void setEdited(bool yes);
    bool getEdited();

    void updateSelected();
    
private:
    ToolPanelListener *getListener();
    void populateList();
    void onSelectionChanged();
    void onAddPressed();
    void onRemovePressed();
    void onUpPressed();
    void onDownPressed();
    void onCopyPressed();
    void onResetPressed();
    void onAreaMaskEnableToggled();
    void onAreaMaskInvertedChanged();
    void onAreaMaskToggleChanged();
    void onShowMaskChanged();
    void updateAreaMask(bool from_mask);
    void maskGet(int idx);
    void maskShow(int idx, bool list_only=false, bool unsub=true);

    void disableListener();
    void enableListener();

    LabMasksContentProvider *cp_;
    std::vector<rtengine::LabCorrectionMask> masks_;
    int selected_;

    rtengine::ProcEvent EvMaskList;
    rtengine::ProcEvent EvHMask;
    rtengine::ProcEvent EvCMask;
    rtengine::ProcEvent EvLMask;
    rtengine::ProcEvent EvMaskBlur;
    rtengine::ProcEvent EvShowMask;
    rtengine::ProcEvent EvAreaMask;

    Gtk::ListViewText *list;
    Gtk::Button *reset;
    Gtk::Button *add;
    Gtk::Button *remove;
    Gtk::Button *up;
    Gtk::Button *down;
    Gtk::Button *copy;
    CurveEditorGroup *maskEditorGroup;
    FlatCurveEditor *hueMask;
    FlatCurveEditor *chromaticityMask;
    FlatCurveEditor *lightnessMask;
    Adjuster *maskBlur;
    Gtk::CheckButton *showMask;
    sigc::connection selectionConn;
    MyExpander *areaMask;
    Gtk::HBox *areaMaskButtonsHb;
    Gtk::ToggleButton *areaMaskToggle;
    Gtk::CheckButton *areaMaskInverted;
    Adjuster *areaMaskX;
    Adjuster *areaMaskY;
    Adjuster *areaMaskWidth;
    Adjuster *areaMaskHeight;
    Adjuster *areaMaskAngle;
    Adjuster *areaMaskFeather;
    Adjuster *areaMaskRoundness;
    Adjuster *areaMaskContrast;
    std::vector<Adjuster *> areaMaskAdjusters;
    std::vector<bool> listenerDisabled;
    rtengine::LabCorrectionMask::AreaMask defaultAreaMask;
    bool listEdited;
};
