/* -*- C++ -*-
 *  
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
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
#include "guiutils.h"
#include "../rtengine/coord.h"

class InspectorBuffer;
class FileCatalog;

class InspectorArea: public Gtk::DrawingArea {
public:
    InspectorArea();
    ~InspectorArea() override;

    /** @brief Mouse movement to a new position
     * @param pos Location of the mouse, in percentage (i.e. [0;1] range) relative to the full size image ; -1,-1 == out of the image
     * @param transform H/V flip and coarse rotation transformation
     */
    void mouseMove (rtengine::Coord2D pos, int transform);

    /** @brief A new image is being flown over
     * @param fullPath Full path of the image that is being hovered inspect, or an empty string if out of any image.
     */
    void switchImage (const Glib::ustring &fullPath);

    /** @brief Set the new coarse rotation transformation
     * @param transform A semi-bitfield coarse transformation using #defines from iimage.h
     */
    // void setTransformation (int transform);

    /** @brief Use this method to flush all image buffer whenever the Inspector panel is hidden
     */
    void flushBuffers ();

    /** @brief Set the inspector on/off
     * @param state true if to activate the Inspector, false to disable it and flush the buffers
     */
    void setActive(bool state);

    /** @brief Get the on/off state
     */
    bool isActive() const
    {
        return active;
    }

    void setInfoText(const Glib::ustring &txt);
    void infoEnabled(bool yes);
    void setZoomFit(bool yes);

    Gtk::SizeRequestMode get_request_mode_vfunc () const override;
    void get_preferred_height_vfunc (int& minimum_height, int& natural_height) const override;
    void get_preferred_width_vfunc (int &minimum_width, int &natural_width) const override;
    void get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const override;
    void get_preferred_width_for_height_vfunc (int height, int &minimum_width, int &natural_width) const override;

    sigc::signal<void> signal_ready() { return sig_ready_; }

private:
    bool on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr) override;
    void deleteBuffers();
    bool doSwitchImage();

    rtengine::Coord center;
    std::vector<InspectorBuffer*> images;
    InspectorBuffer* currImage;
    //double zoom;
    bool active;
    bool first_active_;

    sigc::connection delayconn;
    Glib::ustring next_image_path;

    Glib::ustring info_text_;
    BackBuffer info_bb_;

    sigc::signal<void> sig_ready_;
};


class Inspector: public Gtk::VBox {
public:
    Inspector(FileCatalog *filecatalog);

    void mouseMove(rtengine::Coord2D pos, int transform)
    { ins_.mouseMove(pos, transform); }
    
    void switchImage(const Glib::ustring &fullPath);    
    // void setTransformation(int transform);
    
    void flushBuffers() { ins_.flushBuffers(); }
    
    void setActive(bool state) { ins_.setActive(state); }
    bool isActive() const { return ins_.isActive(); };

    sigc::signal<void> signal_ready() { return ins_.signal_ready(); }
    
private:
    Gtk::HBox *get_toolbar();
    Glib::ustring get_info_text();
    void info_toggled();
    void mode_toggled(Gtk::ToggleButton *b);
    void zoom_toggled(Gtk::ToggleButton *b);
    void cms_toggled();
    bool keyPressed(GdkEventKey *evt);

    FileCatalog *filecatalog_;

    Glib::ustring cur_image_;
    InspectorArea ins_;

    Gtk::ToggleButton *info_;
    Gtk::ToggleButton *jpg_;
    Gtk::ToggleButton *rawlinear_;
    Gtk::ToggleButton *rawfilm_;
    Gtk::ToggleButton *rawshadow_;
    Gtk::ToggleButton *rawclip_;
    Gtk::ToggleButton *zoomfit_;
    Gtk::ToggleButton *zoom11_;
    Gtk::ToggleButton *cms_;

    sigc::connection jpgconn_;
    sigc::connection rawlinearconn_;
    sigc::connection rawfilmconn_;
    sigc::connection rawshadowconn_;
    sigc::connection rawclipconn_;
    sigc::connection zoomfitconn_;
    sigc::connection zoom11conn_;
};
