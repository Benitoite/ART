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
#pragma once

#include <gtkmm.h>
#include "edit.h"


class AreaMask: public EditSubscriber {
public:
    AreaMask();
    
    CursorShape getCursor(const int objectID) override;
    bool mouseOver(const int modifierKey) override;
    bool button1Pressed(const int modifierKey) override;
    bool button1Released() override;
    bool drag1(const int modifierKey) override;

    void updateGeometry(const int fullWidth=-1, const int fullHeight=-1);
    
protected:
    int last_object_;
    double dragged_point_old_angle_;
    double dragged_point_adjuster_angle_;
    rtengine::Coord dragged_center_;
    double center_x_;
    double center_y_;
    double width_;
    double height_;
    double angle_;

    int top_id_;
    int bottom_id_;
    int left_id_;
    int right_id_;
    int rotate_w_id_;
    int rotate_h_id_;
    int center_id_;
};
