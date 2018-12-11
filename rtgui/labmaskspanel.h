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


class LabMasksPanel: public Gtk::VBox {
public:
    LabMasksPanel(Gtk::Widget *child, rtengine::ProcEvent h_mask, rtengine::ProcEvent c_mask, rtengine::ProcEvent l_mask, rtengine::ProcEvent blur, rtengine::ProcEvent area_mask);
    virtual ~LabMasksPanel();

    virtual ToolPanelListener *listener();
    
    virtual bool selectionChanged(int idx);
    virtual bool addPressed();
    virtual bool removePressed(int idx);
    virtual bool copyPressed(int idx);
    virtual bool moveUpPressed(int idx);
    virtual bool moveDownPressed(int idx);

    void setMasks(const std::vector<rtengine::LabCorrectionMask> &masks);
    void getMasks(std::vector<rtengine::LabCorrectionMask> &masks);

private:
    std::vector<rtengine::LabCorrectionMask> masks_;
    int selected_;

    rtengine::ProcEvent EvHMask;
    rtengine::ProcEvent EvCMask;
    rtengine::ProcEvent EvLMask;
    rtengine::ProcEvent EvMaskBlur;
    rtengine::ProcEvent EvAreaMask;

    Gtk::ListViewText *list;
    Gtk::Button *add;
    Gtk::Button *remove;
    Gtk::Button *up;
    Gtk::Button *down;
    Gtk::Button *copy;
    FlatCurveEditor *hueMask;
    FlatCurveEditor *chromaticityMask;
    FlatCurveEditor *lightnessMask;
    Adjuster *maskBlur;
    Gtk::CheckButton *showMask;
    sigc::connection selectionConn;
    MyExpander *areaMask;
    Gtk::ToggleButton *areaMaskToggle;
    Gtk::CheckButton *areaMaskInverted;
    Adjuster *areaMaskX;
    Adjuster *areaMaskY;
    Adjuster *areaMaskWidth;
    Adjuster *areaMaskHeight;
    Adjuster *areaMaskAngle;
    Adjuster *areaMaskFeather;
    Adjuster *areaMaskRoundness;
};
