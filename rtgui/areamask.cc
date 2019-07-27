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

#include "areamask.h"

using rtengine::Coord;
using rtengine::PolarCoord;


AreaMask::AreaMask():
    EditSubscriber(ET_OBJECTS),
    last_object_(-1),
    dragged_point_old_angle_(-1000),
    dragged_point_adjuster_angle_(-1000),
    dragged_center_(0, 0),
    center_x_(0),
    center_y_(0),
    width_(100),
    height_(100),
    angle_(0)
{
    const auto mkline = [](std::vector<Geometry *> &geom) -> int
                        {
                            Line *ret = new Line();
                            ret->innerLineWidth = 2;
                            ret->datum = Geometry::IMAGE;
                            geom.push_back(ret);
                            return geom.size() - 1;
                        };

    const auto mkcircle = [](std::vector<Geometry *> &geom) -> int
                          {
                              Circle *ret = new Circle();
                              ret->datum = Geometry::IMAGE;
                              ret->radiusInImageSpace = false;
                              ret->radius = 6;
                              ret->filled = true;
                              geom.push_back(ret);
                              return geom.size() - 1;
                          };
    
    top_id_ = mkline(visibleGeometry); // top
    bottom_id_ = mkline(visibleGeometry); // bottom
    left_id_ = mkline(visibleGeometry); // left
    right_id_ = mkline(visibleGeometry); // right
    rotate_w_id_ = mkline(visibleGeometry); // rotate_w
    rotate_h_id_ = mkline(visibleGeometry); // rotate_h
    center_id_ = mkcircle(visibleGeometry); // center

    // MouseOver geometry
    mkline(mouseOverGeometry); // top
    mkline(mouseOverGeometry); // bottom
    mkline(mouseOverGeometry); // left
    mkline(mouseOverGeometry); // right
    mkline(mouseOverGeometry); // rotate_w
    mkline(mouseOverGeometry); // rotate_h
    mkcircle(mouseOverGeometry); // center
}


CursorShape AreaMask::getCursor(const int objectID)
{
    if (objectID == top_id_ || objectID == bottom_id_) {
        if (angle_ < -135 || (angle_ >= -45 && angle_ <= 45) || angle_ > 135) {
            return CSMove1DV;
        } else {
            return CSMove1DH;
        }
    } else if (objectID == left_id_ || objectID == right_id_) {
        if (angle_ < -135 || (angle_ >= -45 && angle_ <= 45) || angle_ > 135) {
            return CSMove1DH;
        } else {
            return CSMove1DV;
        }
    } else if (objectID == rotate_w_id_ || objectID == rotate_h_id_) {
        return CSMoveRotate;
    } else if (objectID == center_id_) {
        return CSMove2D;
    } else {
        return CSHandOpen;
    }
}


bool AreaMask::mouseOver(const int modifierKey)
{
    EditDataProvider *provider = getEditProvider();

    if (provider && provider->object != last_object_) {
        if (last_object_ > -1) {
            for (int i = top_id_; i <= center_id_; ++i) {
                visibleGeometry[i]->state = Geometry::NORMAL;
            }
        }

        if (provider->object > -1) {
            for (int i = top_id_; i <= center_id_; ++i) {
                visibleGeometry[i]->state = Geometry::PRELIGHT;
            }
        }

        last_object_ = provider->object;
        return true;
    }

    return false;
}


bool AreaMask::button1Pressed(const int modifierKey)
{
    if (last_object_ < 0) {
        return false;
    }

    EditDataProvider *provider = getEditProvider();

    if (!(modifierKey & GDK_CONTROL_MASK)) {
        // button press is valid (no modifier key)
        PolarCoord pCoord;
        int imW, imH;
        provider->getImageSize(imW, imH);
        double halfSizeW = imW / 2.;
        double halfSizeH = imH / 2.;
        dragged_center_.set(halfSizeW + halfSizeW * (center_x_ / 100), halfSizeH + halfSizeH * center_y_ / 100);

        // trick to get the correct angle (clockwise/counter-clockwise)
        rtengine::Coord p1 = dragged_center_;
        rtengine::Coord p2 = provider->posImage;
        int p = p1.y;
        p1.y = p2.y;
        p2.y = p;
        pCoord = p2 - p1;
        dragged_point_old_angle_ = pCoord.angle;
        dragged_point_adjuster_angle_ = angle_;

        EditSubscriber::action = ES_ACTION_DRAGGING;
        return false;
    } else { // should theoretically always be true
        // this will let this class ignore further drag events
        for (int i = top_id_; i <= center_id_; ++i) {
            visibleGeometry[i]->state = Geometry::NORMAL;
        }

        last_object_ = -1;
        return true;
    }

    return false;
}


bool AreaMask::button1Released()
{
    dragged_point_old_angle_ = -1000.;
    last_object_ = -1;
    EditSubscriber::action = ES_ACTION_NONE;
    for (int i = top_id_; i <= center_id_; ++i) {
        visibleGeometry[i]->state = Geometry::NORMAL;
    }
    return true;
}


bool AreaMask::drag1(const int modifierKey)
{
    // compute the polar coordinate of the mouse position
    EditDataProvider *provider = getEditProvider();
    int imW, imH;
    provider->getImageSize(imW, imH);
    double halfSizeW = imW / 2.;
    double halfSizeH = imH / 2.;

    if (last_object_ == rotate_w_id_ || last_object_ == rotate_h_id_) {
        // Dragging a line to change the angle
        rtengine::Coord currPos = provider->posImage + provider->deltaImage;
        rtengine::Coord centerPos = dragged_center_;

        // trick to get the correct angle (clockwise/counter-clockwise)
        std::swap(centerPos.y, currPos.y);
        
        PolarCoord draggedPoint(currPos - centerPos);
        double deltaAngle = draggedPoint.angle - dragged_point_old_angle_;

        if (deltaAngle > 180.) { // crossing the boundary (0->360)
            deltaAngle -= 360.;
        } else if (deltaAngle < -180.) { // crossing the boundary (360->0)
            deltaAngle += 360.;
        }

        dragged_point_old_angle_ = draggedPoint.angle;
        dragged_point_adjuster_angle_ += deltaAngle;

        if (dragged_point_adjuster_angle_ > 180.) {
            dragged_point_adjuster_angle_ = -360. + dragged_point_adjuster_angle_;
        } else if (dragged_point_adjuster_angle_ < -180.) {
            dragged_point_adjuster_angle_ = 360. - dragged_point_adjuster_angle_;
        }

        if (dragged_point_adjuster_angle_ != angle_) {
            if (dragged_point_adjuster_angle_ < 0) {
                angle_ = 180.0 + dragged_point_adjuster_angle_;
            } else {
                angle_ = dragged_point_adjuster_angle_;
            }
            updateGeometry();
            return true;
        }
    } else if (last_object_ >= top_id_ && last_object_ <= right_id_) {
        // Dragging to resize
        rtengine::Coord currPos = provider->posImage + provider->deltaImage;
        rtengine::Coord centerPos = dragged_center_;

        // trick to get the correct angle (clockwise/counter-clockwise)
        std::swap(centerPos.y, currPos.y);

        PolarCoord draggedPoint(currPos - centerPos);
        double cur_offset_h = draggedPoint.radius * sin((draggedPoint.angle - angle_) / 180.*rtengine::RT_PI);
        double cur_offset_w = draggedPoint.radius * cos((draggedPoint.angle - angle_) / 180.*rtengine::RT_PI);

        double ww = width_, hh = height_;

        if (last_object_ == top_id_ || last_object_ == bottom_id_) {
            double ch = (cur_offset_h - height_ / 2) * (last_object_ == top_id_ ? -1 : 1);
            hh = std::max(ch * 200 / imH, 0.1);
        } else { // left_id_ || right_id_
            double cw = (cur_offset_w - width_ / 2) * (last_object_ == left_id_ ? -1 : 1);
            ww = std::max(cw * 200 / imW, 0.1);
        }

        if (ww != width_ || hh != height_) {
            width_ = ww;
            height_ = hh;
            updateGeometry();
            return true;
        }
    } else if (last_object_ == center_id_) {
        // Dragging the circle to change the center
        dragged_center_ += provider->deltaPrevImage;
        rtengine::Coord currPos = dragged_center_;
        currPos.clip(imW, imH);
        double cx = (double(currPos.x) - halfSizeW) / halfSizeW * 100.;
        double cy = (double(currPos.y) - halfSizeH) / halfSizeH * 100.;

        if (cx != center_x_ || cy != center_y_) {
            center_x_ = cx;
            center_y_ = cy;
            updateGeometry();
            return true;
        }
    }

    return false;
}


void AreaMask::updateGeometry(const int fullWidth, const int fullHeight)
{
    EditDataProvider* provider = getEditProvider();

    if (!provider) {
        return;
    }

    int imW = 0;
    int imH = 0;
    if (fullWidth != -1 && fullHeight != -1) {
        imW = fullWidth;
        imH = fullHeight;
    } else {
        provider->getImageSize(imW, imH);
        if (!imW || !imH) {
            return;
        }
    }

    const auto width = width_ / 100 * imW;
    const auto height = height_ / 100 * imH;
    rtengine::Coord origin(imW / 2 + center_x_ * imW / 200, imH / 2 + center_y_ * imH / 200);

    const auto update_border =
        [&](Geometry *geometry, int direction)
        {
            const auto line = static_cast<Line *>(geometry);
            float w = width / 2;
            float h = height / 2;
            int sx, sy;
            if (direction == top_id_) {
                sx = sy = 1;
            } else if (direction == bottom_id_) {
                sx = 1;
                sy = -1;
            } else if (direction == left_id_) {
                sx = -1;
                sy = 1;
            } else {
                sx = 1;
                sy = 1;
            }
            PolarCoord begin, end;
            if (direction == top_id_ || direction == bottom_id_) {
                begin = Coord(-sx * w, sy * h);
                end = Coord(sx * w, sy * h);
            } else {
                begin = Coord(sx * w, -sy * h);
                end = Coord(sx * w, sy * h);
            }
            double r, a;
            begin.get(r, a);
            begin.set(r, a - angle_);
            end.get(r, a);
            end.set(r, a - angle_);

            line->begin = begin;
            line->end = end;
            line->begin += origin;
            line->end += origin;
        };

    const auto update_inner =
        [&](Geometry *geometry, int direction)
        {
            const auto line = static_cast<Line *>(geometry);
            PolarCoord begin, end;
            float w = width / 4;
            float h = height / 4;
            if (direction == rotate_w_id_) {
                begin = Coord(-w, 0);
                end = Coord(w, 0);
            } else {
                begin = Coord(0, -h);
                end = Coord(0, h);
            }
            double r, a;
            begin.get(r, a);
            begin.set(r, a - angle_);
            end.get(r, a);
            end.set(r, a - angle_);
            
            line->begin = begin;
            line->end = end;
            line->begin += origin;
            line->end += origin;
        };

    const auto update_center =
        [&](Geometry *geometry)
        {
            const auto circle = static_cast<Circle *>(geometry);
            circle->center = origin;
        };

    for (int dir = top_id_; dir <= right_id_; ++dir) {
        update_border(visibleGeometry[dir], dir);
        update_border(mouseOverGeometry[dir], dir);
    }

    for (int dir = rotate_w_id_; dir <= rotate_h_id_; ++dir) {
        update_inner(visibleGeometry[dir], dir);
        update_inner(mouseOverGeometry[dir], dir);
    }

    update_center(visibleGeometry[center_id_]);
    update_center(mouseOverGeometry[center_id_]);
}
