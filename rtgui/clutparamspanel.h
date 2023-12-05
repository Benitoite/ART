/* -*- C++ -*-
 *
 *  This file is part of ART.
 *
 *  Copyright 2023 Alberto Griggio <alberto.griggio@gmail.com>
 *
 *  ART is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ART is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ART.  If not, see <http://www.gnu.org/licenses/>.
 */
#pragma once

#include <gtkmm.h>
#include "toolpanel.h"
#include "adjuster.h"
#include "../rtengine/clutparams.h"


class CLUTParamsPanel: public Gtk::VBox, public AdjusterListener {
public:
    CLUTParamsPanel();

    void setParams(const std::vector<rtengine::CLUTParamDescriptor> &params);
    void setValue(const std::vector<double> &val);
    std::vector<double> getValue() const;
    
    sigc::signal<void> signal_changed() { return sig_changed_; }

    void adjusterChanged(Adjuster *a, double v) override { emit_signal(); }
    void adjusterAutoToggled(Adjuster *a, bool v) override {}
    
private:
    void emit_signal();
    
    sigc::signal<void> sig_changed_;
    bool sig_blocked_;
    std::vector<rtengine::CLUTParamDescriptor> params_;
    std::vector<Gtk::Widget *> widgets_;
};
