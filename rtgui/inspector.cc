/*
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

#include <gtkmm.h>
#include <iomanip>

#include "inspector.h"
#include "guiutils.h"
#include "cursormanager.h"
#include "guiutils.h"
#include "options.h"
#include "multilangmgr.h"
#include "filecatalog.h"
#include "../rtengine/previewimage.h"
#include "../rtengine/imagedata.h"

extern Options options;

//-----------------------------------------------------------------------------
// InspectorBuffer
//-----------------------------------------------------------------------------

class InspectorBuffer {
//private:
//    int infoFromImage (const Glib::ustring& fname);

public:
    BackBuffer imgBuffer;
    Glib::ustring imgPath;
    int currTransform;  // coarse rotation from RT, not from shot orientation
    bool fromRaw;

    explicit InspectorBuffer(const Glib::ustring &imgagePath, int width=-1, int height=-1);
    //~InspectorBuffer();
};

InspectorBuffer::InspectorBuffer(const Glib::ustring &imagePath, int width, int height) : currTransform(0), fromRaw(false)
{
    if (!imagePath.empty() && Glib::file_test(imagePath, Glib::FILE_TEST_EXISTS) && !Glib::file_test(imagePath, Glib::FILE_TEST_IS_DIR)) {
        imgPath = imagePath;

        // generate thumbnail image
        Glib::ustring ext = getExtension (imagePath);

        if (ext == "") {
            imgPath.clear();
            return;
        }

        rtengine::PreviewImage pi(imagePath, ext, width, height, options.thumbnail_inspector_enable_cms);
        Cairo::RefPtr<Cairo::ImageSurface> imageSurface = pi.getImage();

        if (imageSurface) {
            imgBuffer.setSurface(imageSurface);
            fromRaw = true;
        } else {
            imgPath.clear();
        }
    }
}


//-----------------------------------------------------------------------------
// InspectorArea
//-----------------------------------------------------------------------------

InspectorArea::InspectorArea():
    currImage(nullptr),
    active(false),
    first_active_(true),
    info_text_("")
{
    Glib::RefPtr<Gtk::StyleContext> style = get_style_context();
    set_name("Inspector");
}


InspectorArea::~InspectorArea()
{
    deleteBuffers();
}


bool InspectorArea::on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr)
{
    Glib::RefPtr<Gdk::Window> win = get_window();

    if (!win) {
        return false;
    }

    if (!active) {
        active = true;
    }


    // cleanup the region


    if (currImage && currImage->imgBuffer.surfaceCreated()) {
        // this will eventually create/update the off-screen pixmap

        // compute the displayed area
        rtengine::Coord availableSize;
        rtengine::Coord topLeft;
        rtengine::Coord displayedSize;
        rtengine::Coord dest(0, 0);
        availableSize.x = win->get_width();
        availableSize.y = win->get_height();
        int imW = currImage->imgBuffer.getWidth();
        int imH = currImage->imgBuffer.getHeight();

        if (imW < availableSize.x) {
            // center the image in the available space along X
            topLeft.x = 0;
            displayedSize.x = availableSize.x;
            dest.x = (availableSize.x - imW) / 2;
        } else {
            // partial image display
            // double clamp
            topLeft.x = center.x + availableSize.x / 2;
            topLeft.x = rtengine::min<int>(topLeft.x, imW);
            topLeft.x -= availableSize.x;
            topLeft.x = rtengine::max<int>(topLeft.x, 0);
        }

        if (imH < availableSize.y) {
            // center the image in the available space along Y
            topLeft.y = 0;
            displayedSize.y = availableSize.y;
            dest.y = (availableSize.y - imH) / 2;
        } else {
            // partial image display
            // double clamp
            topLeft.y = center.y + availableSize.y / 2;
            topLeft.y = rtengine::min<int>(topLeft.y, imH);
            topLeft.y -= availableSize.y;
            topLeft.y = rtengine::max<int>(topLeft.y, 0);
        }

        //printf("center: %d, %d   (img: %d, %d)  (availableSize: %d, %d)  (topLeft: %d, %d)\n", center.x, center.y, imW, imH, availableSize.x, availableSize.y, topLeft.x, topLeft.y);

        // define the destination area
        currImage->imgBuffer.setDrawRectangle(win, dest.x, dest.y, rtengine::min<int>(availableSize.x - dest.x, imW), rtengine::min<int>(availableSize.y - dest.y, imH), false);
        currImage->imgBuffer.setSrcOffset(topLeft.x, topLeft.y);

        if (!currImage->imgBuffer.surfaceCreated()) {
            return false;
        }

        // Draw!

        Gdk::RGBA c;
        Glib::RefPtr<Gtk::StyleContext> style = get_style_context();

        // draw the background
        style->render_background(cr, 0, 0, get_width(), get_height());

        /* --- old method
        c = style->get_background_color (Gtk::STATE_FLAG_NORMAL);
        cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
        cr->set_line_width (0);
        cr->rectangle (0, 0, availableSize.x, availableSize.y);
        cr->fill ();
        */

        currImage->imgBuffer.copySurface(win);

        // draw the frame
        c = style->get_border_color (Gtk::STATE_FLAG_NORMAL);
        cr->set_source_rgb (c.get_red(), c.get_green(), c.get_blue());
        cr->set_line_width (1);
        cr->rectangle (0.5, 0.5, availableSize.x - 1, availableSize.y - 1);
        cr->stroke ();

        if (options.thumbnail_inspector_show_info && info_text_ != "") {
            info_bb_.copySurface(cr);
        }
    }

    if (first_active_) {
        first_active_ = false;
        sig_ready_.emit();
    }

    return true;
}


void InspectorArea::mouseMove (rtengine::Coord2D pos, int transform)
{
    if (!active) {
        return;
    }

    if (currImage) {
        center.set(int(rtengine::LIM01(pos.x)*double(currImage->imgBuffer.getWidth())), int(rtengine::LIM01(pos.y)*double(currImage->imgBuffer.getHeight())));
    } else {
        center.set(0, 0);
    }

    queue_draw();
}


void InspectorArea::switchImage(const Glib::ustring &fullPath)
{
    if (!active) {
        return;
    }

    if (delayconn.connected()) {
        delayconn.disconnect();
    }

    next_image_path = fullPath;
    if (!options.inspectorDelay) {
        doSwitchImage();
    } else {
        delayconn = Glib::signal_timeout().connect(sigc::mem_fun(*this, &InspectorArea::doSwitchImage), options.inspectorDelay);
    }
}


bool InspectorArea::doSwitchImage()
{
    Glib::ustring fullPath = next_image_path;
    
    // we first check the size of the list, it may have been changed in Preference
    if (images.size() > size_t(options.maxInspectorBuffers)) {
        // deleting the last entries
        for (size_t i = images.size() - 1; i > size_t(options.maxInspectorBuffers - 1); --i) {
            delete images.at(i);
            images.at(i) = nullptr;
        }

        // resizing down
        images.resize(options.maxInspectorBuffers);
    }

    if (fullPath.empty()) {
        currImage = nullptr;
//        queue_draw();
    } else {
        bool found = false;

        for (size_t i = 0; i < images.size(); ++i) {
            if (images.at(i) != nullptr && images.at(i)->imgPath == fullPath) {
                currImage = images.at(i);

                // rolling the list 1 step to the beginning
                for (size_t j = i; j < images.size() - 1; ++j) {
                    images.at(j) = images.at(j + 1);
                }

                images.at(images.size() - 1) = currImage; // move the last used image to the tail
                found = true;
                break;
            }
        }

        if (!found) {
            if (images.size() == size_t(options.maxInspectorBuffers)) {
                // The list is full, delete the first entry
                delete images.at(0);
                images.erase(images.begin());
            }

            Glib::RefPtr<Gdk::Window> win = get_window();
            int width = -1, height = -1;
            if (win && options.thumbnail_inspector_zoom_fit) {
                width = win->get_width();
                height = win->get_height();
            }
            
            // Loading a new image
            InspectorBuffer *iBuffer = new InspectorBuffer(fullPath, width, height);

            // and add it to the tail
            if (!iBuffer->imgPath.empty()) {
                images.push_back(iBuffer);
                currImage = images.at(images.size() - 1);
            } else {
                delete iBuffer;
                currImage = nullptr;
            }
        }
    }

    queue_draw();

    return true;
}

void InspectorArea::deleteBuffers ()
{
    for (size_t i = 0; i < images.size(); ++i) {
        if (images.at(i) != nullptr) {
            delete images.at(i);
            images.at(i) = nullptr;
        }
    }

    images.resize(0);
    currImage = nullptr;
}

void InspectorArea::flushBuffers ()
{
    if (!active) {
        return;
    }

    deleteBuffers();
}

void InspectorArea::setActive(bool state)
{
    if (!state) {
        flushBuffers();
    }

    active = state;
    if (!active) {
        first_active_ = true;
    }
}


Gtk::SizeRequestMode InspectorArea::get_request_mode_vfunc () const
{
    return Gtk::SIZE_REQUEST_CONSTANT_SIZE;
}

void InspectorArea::get_preferred_height_vfunc (int &minimum_height, int &natural_height) const
{
    minimum_height= 50 * RTScalable::getScale();
    natural_height = 300 * RTScalable::getScale();
}

void InspectorArea::get_preferred_width_vfunc (int &minimum_width, int &natural_width) const
{
    minimum_width = 50 * RTScalable::getScale();
    natural_width = 200 * RTScalable::getScale();
}

void InspectorArea::get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const
{
    get_preferred_height_vfunc(minimum_height, natural_height);
}

void InspectorArea::get_preferred_width_for_height_vfunc (int height, int &minimum_width, int &natural_width) const
{
    get_preferred_width_vfunc (minimum_width, natural_width);
}


void InspectorArea::setInfoText(const Glib::ustring &text)
{
    info_text_ = text;

    Glib::RefPtr<Pango::Context> context = get_pango_context();
    Pango::FontDescription fontd(get_style_context()->get_font());

    // update font
    fontd.set_weight(Pango::WEIGHT_BOLD);
    fontd.set_size(10 * Pango::SCALE);
    context->set_font_description(fontd);

    // create text layout
    Glib::RefPtr<Pango::Layout> ilayout = create_pango_layout("");
    ilayout->set_markup(text);

    // get size of the text block
    int iw, ih;
    ilayout->get_pixel_size(iw, ih);

    // create BackBuffer
    info_bb_.setDrawRectangle(Cairo::FORMAT_ARGB32, 0, 0, iw + 16, ih + 16, true);
    info_bb_.setDestPosition(8, 8);

    Cairo::RefPtr<Cairo::Context> cr = info_bb_.getContext();

    // cleaning the back buffer (make it full transparent)
    cr->set_source_rgba(0., 0., 0., 0.);
    cr->set_operator(Cairo::OPERATOR_CLEAR);
    cr->paint();
    cr->set_operator(Cairo::OPERATOR_OVER);

    // paint transparent black background
    cr->set_source_rgba(0., 0., 0., 0.5);
    cr->paint ();

    // paint text
    cr->set_source_rgb(1.0, 1.0, 1.0);
    cr->move_to(8, 8);
    ilayout->add_to_cairo_context(cr);
    cr->fill();
}


void InspectorArea::infoEnabled(bool yes)
{
    if (options.thumbnail_inspector_show_info != yes) {
        options.thumbnail_inspector_show_info = yes;
        queue_draw();
    }
}


//-----------------------------------------------------------------------------
// Inspector
//-----------------------------------------------------------------------------

Inspector::Inspector(FileCatalog *filecatalog):
    filecatalog_(filecatalog)
{
    pack_start(ins_);
    pack_start(*get_toolbar(), Gtk::PACK_SHRINK, 2);
    show_all_children();

    signal_key_press_event().connect(sigc::mem_fun(*this, &Inspector::keyPressed));
}


bool Inspector::keyPressed(GdkEventKey *evt)
{
    if (filecatalog_) {
        return filecatalog_->handleShortcutKey(evt);
    }
    return false;    
}


void Inspector::switchImage(const Glib::ustring &fullPath)
{
    cur_image_ = fullPath;
    if (info_->get_active()) {
        ins_.setInfoText(get_info_text());
    }
    ins_.switchImage(fullPath);
}


Gtk::HBox *Inspector::get_toolbar()
{
    Gtk::HBox *tb = Gtk::manage(new Gtk::HBox());
    tb->pack_start(*Gtk::manage(new Gtk::Label("")), Gtk::PACK_EXPAND_WIDGET, 2);

    const auto add_tool =
        [&](const char *icon, const char *tip=nullptr) -> Gtk::ToggleButton *
        {
            Gtk::ToggleButton *ret = Gtk::manage(new Gtk::ToggleButton());
            ret->add(*Gtk::manage(new RTImage(icon)));
            ret->set_relief(Gtk::RELIEF_NONE);
            if (tip) {
                ret->set_tooltip_markup(M(tip));
            }
            tb->pack_start(*ret, Gtk::PACK_SHRINK, 2);
            return ret;
        };
       
    info_ = add_tool("info.png", "INSPECTOR_INFO");

    tb->pack_start(*Gtk::manage(new Gtk::VSeparator()), Gtk::PACK_SHRINK, 4);

    jpg_ = add_tool("wb-camera.png", "INSPECTOR_PREVIEW");
    rawlinear_ = add_tool("raw-linear-curve.png", "INSPECTOR_RAW_LINEAR");
    rawfilm_ = add_tool("raw-film-curve.png", "INSPECTOR_RAW_FILM");
    rawshadow_ = add_tool("raw-shadow-curve.png", "INSPECTOR_RAW_SHADOW");
    rawclip_ = add_tool("raw-clip-curve.png", "INSPECTOR_RAW_CLIP");

    tb->pack_start(*Gtk::manage(new Gtk::VSeparator()), Gtk::PACK_SHRINK, 4);

    zoomfit_ = add_tool("magnifier-fit.png", "INSPECTOR_ZOOM_FIT");
    zoom11_ = add_tool("magnifier-1to1.png", "INSPECTOR_ZOOM_11");

    tb->pack_start(*Gtk::manage(new Gtk::VSeparator()), Gtk::PACK_SHRINK, 4);

    cms_ = add_tool("gamut-softproof.png", "INSPECTOR_ENABLE_CMS");

    //------------------------------------------------------------------------
    
    info_->set_active(options.thumbnail_inspector_show_info);
    info_->signal_toggled().connect(sigc::mem_fun(*this, &Inspector::info_toggled));
    bool use_jpg = options.rtSettings.thumbnail_inspector_mode == rtengine::Settings::ThumbnailInspectorMode::JPEG;
    jpg_->set_active(use_jpg);
    rawlinear_->set_active(!use_jpg && options.rtSettings.thumbnail_inspector_raw_curve == rtengine::Settings::ThumbnailInspectorRawCurve::LINEAR);
    rawfilm_->set_active(!use_jpg && options.rtSettings.thumbnail_inspector_raw_curve == rtengine::Settings::ThumbnailInspectorRawCurve::FILM);
    rawshadow_->set_active(!use_jpg && options.rtSettings.thumbnail_inspector_raw_curve == rtengine::Settings::ThumbnailInspectorRawCurve::SHADOW_BOOST);
    rawclip_->set_active(!use_jpg && options.rtSettings.thumbnail_inspector_raw_curve == rtengine::Settings::ThumbnailInspectorRawCurve::RAW_CLIPPING);

    zoomfit_->set_active(options.thumbnail_inspector_zoom_fit);
    zoom11_->set_active(!options.thumbnail_inspector_zoom_fit);

    cms_->set_active(options.thumbnail_inspector_enable_cms);

    jpgconn_ = jpg_->signal_toggled().connect(sigc::bind(sigc::mem_fun(*this, &Inspector::mode_toggled), jpg_));
    rawlinearconn_ = rawlinear_->signal_toggled().connect(sigc::bind(sigc::mem_fun(*this, &Inspector::mode_toggled), rawlinear_));
    rawfilmconn_ = rawfilm_->signal_toggled().connect(sigc::bind(sigc::mem_fun(*this, &Inspector::mode_toggled), rawfilm_));
    rawshadowconn_ = rawshadow_->signal_toggled().connect(sigc::bind(sigc::mem_fun(*this, &Inspector::mode_toggled), rawshadow_));
    rawclipconn_ = rawclip_->signal_toggled().connect(sigc::bind(sigc::mem_fun(*this, &Inspector::mode_toggled), rawclip_));

    zoomfitconn_ = zoomfit_->signal_toggled().connect(sigc::bind(sigc::mem_fun(*this, &Inspector::zoom_toggled), zoomfit_));
    zoom11conn_ = zoom11_->signal_toggled().connect(sigc::bind(sigc::mem_fun(*this, &Inspector::zoom_toggled), zoom11_));

    cms_->signal_toggled().connect(sigc::mem_fun(*this, &Inspector::cms_toggled));

    return tb;
}


void Inspector::info_toggled()
{
    if (!info_->get_active()) {
        ins_.infoEnabled(false);
        return;
    }

    ins_.setInfoText(get_info_text());
    ins_.infoEnabled(true);
}


Glib::ustring Inspector::get_info_text()
{    
    rtengine::FramesData meta(cur_image_);

    Glib::ustring infoString;
    Glib::ustring expcomp;

    if (meta.hasExif()) {
        infoString = Glib::ustring::compose ("%1 + %2\n<span size=\"small\">f/</span><span size=\"large\">%3</span>  <span size=\"large\">%4</span><span size=\"small\">s</span>  <span size=\"small\">%5</span><span size=\"large\">%6</span>  <span size=\"large\">%7</span><span size=\"small\">mm</span>",
                                              Glib::ustring(meta.getMake() + " " + meta.getModel()),
                                              Glib::ustring(meta.getLens()),
                                              Glib::ustring(meta.apertureToString (meta.getFNumber())),
                                              Glib::ustring(meta.shutterToString (meta.getShutterSpeed())),
                                              M("QINFO_ISO"), meta.getISOSpeed(),
                                              Glib::ustring::format(std::setw (3), std::fixed, std::setprecision (2), meta.getFocalLen()));

        expcomp = Glib::ustring(meta.expcompToString(meta.getExpComp(), true)); // maskZeroexpcomp

        if (!expcomp.empty ()) {
            infoString = Glib::ustring::compose ("%1  <span size=\"large\">%2</span><span size=\"small\">EV</span>",
                                                  infoString,
                                                  expcomp /*Glib::ustring(meta.expcompToString(meta.getExpComp()))*/);
        }

        infoString = Glib::ustring::compose("%1\n<span size=\"small\">%2</span><span>%3</span>",
                                              infoString,
                                              escapeHtmlChars(Glib::path_get_dirname(cur_image_)) + G_DIR_SEPARATOR_S,
                                              escapeHtmlChars(Glib::path_get_basename(cur_image_)));

        // int ww = ipc->getFullWidth();
        // int hh = ipc->getFullHeight();
        // //megapixels
        // infoString = Glib::ustring::compose ("%1\n<span size=\"small\">%2 MP (%3x%4)</span>",
        //                                      infoString,
        //                                      Glib::ustring::format (std::setw (4), std::fixed, std::setprecision (1), (float)ww * hh / 1000000),
        //                                      ww, hh);

        // //adding special characteristics
        // bool isHDR = meta.getHDR();
        // bool isPixelShift = meta.getPixelShift();
        // unsigned int numFrames = meta.getFrameCount();
        // if (isHDR) {
        //     infoString = Glib::ustring::compose ("%1\n" + M("QINFO_HDR"), infoString, numFrames);
        //     if (numFrames == 1) {
        //         int sampleFormat = meta.getSampleFormat();
        //         infoString = Glib::ustring::compose ("%1 / %2", infoString, M(Glib::ustring::compose("SAMPLEFORMAT_%1", sampleFormat)));
        //     }
        // } else if (isPixelShift) {
        //     infoString = Glib::ustring::compose ("%1\n" + M("QINFO_PIXELSHIFT"), infoString, numFrames);
        // } else if (numFrames > 1) {
        //     infoString = Glib::ustring::compose ("%1\n" + M("QINFO_FRAMECOUNT"), infoString, numFrames);
        // }
    } else {
        infoString = M("QINFO_NOEXIF");
    }
    return infoString;
}


void Inspector::mode_toggled(Gtk::ToggleButton *b)
{
    ConnectionBlocker blockj(jpgconn_);
    ConnectionBlocker blockl(rawlinearconn_);
    ConnectionBlocker blockf(rawfilmconn_);
    ConnectionBlocker blocks(rawshadowconn_);
    ConnectionBlocker blockc(rawclipconn_);

    if (!b->get_active()) {
        b->set_active(true);
    } else {
        jpg_->set_active(false);
        rawlinear_->set_active(false);
        rawfilm_->set_active(false);
        rawshadow_->set_active(false);
        rawclip_->set_active(false);
        b->set_active(true);

        if (jpg_->get_active()) {
            options.rtSettings.thumbnail_inspector_mode = rtengine::Settings::ThumbnailInspectorMode::JPEG;
        } else if (rawlinear_->get_active()) {
            options.rtSettings.thumbnail_inspector_mode = rtengine::Settings::ThumbnailInspectorMode::RAW;
            options.rtSettings.thumbnail_inspector_raw_curve = rtengine::Settings::ThumbnailInspectorRawCurve::LINEAR;
        } else if (rawfilm_->get_active()) {
            options.rtSettings.thumbnail_inspector_mode = rtengine::Settings::ThumbnailInspectorMode::RAW;
            options.rtSettings.thumbnail_inspector_raw_curve = rtengine::Settings::ThumbnailInspectorRawCurve::FILM;
        } else if (rawshadow_->get_active()) {
            options.rtSettings.thumbnail_inspector_mode = rtengine::Settings::ThumbnailInspectorMode::RAW;
            options.rtSettings.thumbnail_inspector_raw_curve = rtengine::Settings::ThumbnailInspectorRawCurve::SHADOW_BOOST;
        } else if (rawclip_->get_active()) {
            options.rtSettings.thumbnail_inspector_mode = rtengine::Settings::ThumbnailInspectorMode::RAW;
            options.rtSettings.thumbnail_inspector_raw_curve = rtengine::Settings::ThumbnailInspectorRawCurve::RAW_CLIPPING;
        }

        ins_.flushBuffers();
        ins_.switchImage(cur_image_);
    }
}


void Inspector::zoom_toggled(Gtk::ToggleButton *b)
{
    ConnectionBlocker blockf(zoomfitconn_);
    ConnectionBlocker block1(zoom11conn_);

    if (!b->get_active()) {
        b->set_active(true);
    } else {
        zoomfit_->set_active(false);
        zoom11_->set_active(false);
        b->set_active(true);

        options.thumbnail_inspector_zoom_fit = zoomfit_->get_active();
        ins_.flushBuffers();
        ins_.switchImage(cur_image_);
    }
}


void Inspector::cms_toggled()
{
    options.thumbnail_inspector_enable_cms = cms_->get_active();
    ins_.flushBuffers();
    ins_.switchImage(cur_image_);
}
