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
 *  along with RawTherapee.  If not, see <https://www.gnu.org/licenses/>.
 */
#include "histogrampanel.h"
#include "multilangmgr.h"
#include "guiutils.h"
#include "options.h"
#include <cstring>
#include <cmath>
#include "../rtengine/array2D.h"
#include "../rtengine/LUT.h"
#include "rtimage.h"
#include "../rtengine/color.h"

using namespace rtengine;

constexpr float HistogramArea::MAX_BRIGHT;
constexpr float HistogramArea::MIN_BRIGHT;

using ScopeType = Options::ScopeType;

namespace {

constexpr float rgb_R[3] = { 1.f, 0.f, 0.f };
constexpr float rgb_G[3] = { 0.f, 1.f, 0.f };
constexpr float rgb_B[3] = { 0.f, 0.25f, 1.f };

void set_arr(array2D<int> &dst, const array2D<int> &src)
{
    dst(src.width(), src.height());
    for (int y = 0; y < src.height(); ++y) {
        for (int x = 0; x < src.width(); ++x) {
            dst[y][x] = src[y][x];
        }
    }
}

constexpr double padding = 0;

} // namespace

//
//
// HistogramPanel
HistogramPanel::HistogramPanel () :
    panel_listener(nullptr)
{
    setExpandAlignProperties(this, true, true, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    set_name("HistogramPanel");

    histogramArea = Gtk::manage (new HistogramArea (this));
    setExpandAlignProperties(histogramArea, true, true, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);

    histogramRGBAreaHori.reset(new HistogramRGBAreaHori());
    setExpandAlignProperties(histogramRGBAreaHori.get(), true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_END);

    histogramRGBAreaVert.reset(new HistogramRGBAreaVert());
    setExpandAlignProperties(histogramRGBAreaVert.get(), false, true, Gtk::ALIGN_END, Gtk::ALIGN_FILL);

    switch (options.histogramScopeType) {
        case ScopeType::NONE:
        case ScopeType::HISTOGRAM_RAW:
        case ScopeType::VECTORSCOPE_HC:
        case ScopeType::VECTORSCOPE_HS:
            histogramRGBArea = nullptr;
            break;
        case ScopeType::PARADE:
        case ScopeType::WAVEFORM:
            histogramRGBArea = histogramRGBAreaVert.get();
            break;
        case ScopeType::HISTOGRAM:
            histogramRGBArea = histogramRGBAreaHori.get();
            break;
    }

    // connecting the two childs
    histogramArea->signal_factor_changed().connect( sigc::mem_fun(*histogramRGBAreaHori, &HistogramRGBArea::factorChanged) );
    histogramArea->signal_factor_changed().connect( sigc::mem_fun(*histogramRGBAreaVert, &HistogramRGBArea::factorChanged) );

    gfxGrid = Gtk::manage (new Gtk::Grid ());
    gfxGrid->set_row_spacing(1);
    gfxGrid->set_column_spacing(1);
    gfxGrid->add(*histogramArea);
    gfxGrid->attach_next_to(
        *histogramRGBAreaVert, *histogramArea,
        options.histogramPosition == Options::HISTOGRAM_POS_LEFT ? Gtk::POS_RIGHT : Gtk::POS_LEFT,
        1,
        1
    );
    gfxGrid->attach_next_to(*histogramRGBAreaHori, *histogramArea, Gtk::POS_BOTTOM, 1, 1);
    histogramRGBAreaHori->set_no_show_all();
    histogramRGBAreaVert->set_no_show_all();

    redImage   = new RTImage ("histogram-red-on-small.png");
    greenImage = new RTImage ("histogram-green-on-small.png");
    blueImage  = new RTImage ("histogram-blue-on-small.png");
    valueImage = new RTImage ("histogram-silver-on-small.png");
    chroImage  = new RTImage ("histogram-gold-on-small.png");
    barImage   = new RTImage ("histogram-bar-on-small.png");

    redImage_g   = new RTImage ("histogram-red-off-small.png");
    greenImage_g = new RTImage ("histogram-green-off-small.png");
    blueImage_g  = new RTImage ("histogram-blue-off-small.png");
    valueImage_g = new RTImage ("histogram-silver-off-small.png");
    chroImage_g  = new RTImage ("histogram-gold-off-small.png");
    barImage_g   = new RTImage ("histogram-bar-off-small.png");

    mode_images_[0] = new RTImage ("histogram-mode-linear-small.png");
    mode_images_[1] = new RTImage ("histogram-mode-logx-small.png");
    mode_images_[2] = new RTImage ("histogram-mode-logxy-small.png");
    mode_tips_[0] = M("HISTOGRAM_TOOLTIP_MODE_LINEAR");
    mode_tips_[1] = M("HISTOGRAM_TOOLTIP_MODE_LOG_X");
    mode_tips_[2] = M("HISTOGRAM_TOOLTIP_MODE_LOG_XY");

    Gtk::Image* histImage = Gtk::manage(new RTImage("histogram-type-histogram-small.png"));
    Gtk::Image* histRawImage = Gtk::manage(new RTImage("histogram-type-histogram-raw-small.png"));
    Gtk::Image* paradeImage = Gtk::manage(new RTImage("histogram-type-parade-small.png"));
    Gtk::Image* waveImage = Gtk::manage(new RTImage("histogram-type-waveform-small.png"));
    Gtk::Image* vectHcImage = Gtk::manage(new RTImage("histogram-type-vectorscope-hc-small.png"));
    Gtk::Image* vectHsImage = Gtk::manage(new RTImage("histogram-type-vectorscope-hs-small.png"));

    showRed   = Gtk::manage (new Gtk::ToggleButton ());
    showGreen = Gtk::manage (new Gtk::ToggleButton ());
    showBlue  = Gtk::manage (new Gtk::ToggleButton ());
    showValue = Gtk::manage (new Gtk::ToggleButton ());
    showChro  = Gtk::manage (new Gtk::ToggleButton ());
    showMode  = Gtk::manage (new Gtk::Button ());
    showBAR   = Gtk::manage (new Gtk::ToggleButton ());
    scopeOptions = Gtk::manage (new Gtk::ToggleButton ());

    Gtk::RadioButtonGroup scopeTypeGroup;
    scopeHistBtn = Gtk::manage(new Gtk::RadioButton(scopeTypeGroup));
    scopeHistRawBtn = Gtk::manage(new Gtk::RadioButton(scopeTypeGroup));
    scopeParadeBtn = Gtk::manage(new Gtk::RadioButton(scopeTypeGroup));
    scopeWaveBtn = Gtk::manage(new Gtk::RadioButton(scopeTypeGroup));
    scopeVectHcBtn = Gtk::manage(new Gtk::RadioButton(scopeTypeGroup));
    scopeVectHsBtn = Gtk::manage(new Gtk::RadioButton(scopeTypeGroup));
    scopeHistBtn->set_mode(false);
    scopeHistRawBtn->set_mode(false);
    scopeParadeBtn->set_mode(false);
    scopeWaveBtn->set_mode(false);
    scopeVectHcBtn->set_mode(false);
    scopeVectHsBtn->set_mode(false);

    showRed->set_name("histButton");
    showRed->set_can_focus(false);
    showGreen->set_name("histButton");
    showGreen->set_can_focus(false);
    showBlue->set_name("histButton");
    showBlue->set_can_focus(false);
    showValue->set_name("histButton");
    showValue->set_can_focus(false);
    showChro->set_name("histButton");
    showChro->set_can_focus(false);
    showMode->set_name("histButton");
    showMode->set_can_focus(false);
    scopeOptions->set_name("histButton");
    scopeOptions->set_can_focus(false);
    showBAR->set_name("histButton");
    showBAR->set_can_focus(false);
    scopeHistBtn->set_name("histButton");
    scopeHistBtn->set_can_focus(false);
    scopeHistRawBtn->set_name("histButton");
    scopeHistRawBtn->set_can_focus(false);
    scopeParadeBtn->set_name("histButton");
    scopeParadeBtn->set_can_focus(false);
    scopeWaveBtn->set_name("histButton");
    scopeWaveBtn->set_can_focus(false);
    scopeVectHcBtn->set_name("histButton");
    scopeVectHcBtn->set_can_focus(false);
    scopeVectHsBtn->set_name("histButton");
    scopeVectHsBtn->set_can_focus(false);

    showRed->set_relief (Gtk::RELIEF_NONE);
    showGreen->set_relief (Gtk::RELIEF_NONE);
    showBlue->set_relief (Gtk::RELIEF_NONE);
    showValue->set_relief (Gtk::RELIEF_NONE);
    showChro->set_relief (Gtk::RELIEF_NONE);
    showMode->set_relief (Gtk::RELIEF_NONE);
    scopeOptions->set_relief (Gtk::RELIEF_NONE);
    showBAR->set_relief (Gtk::RELIEF_NONE);
    scopeHistBtn->set_relief (Gtk::RELIEF_NONE);
    scopeHistRawBtn->set_relief (Gtk::RELIEF_NONE);
    scopeParadeBtn->set_relief (Gtk::RELIEF_NONE);
    scopeWaveBtn->set_relief (Gtk::RELIEF_NONE);
    scopeVectHcBtn->set_relief (Gtk::RELIEF_NONE);
    scopeVectHsBtn->set_relief (Gtk::RELIEF_NONE);

    showRed->set_tooltip_text   (M("HISTOGRAM_TOOLTIP_R"));
    showGreen->set_tooltip_text (M("HISTOGRAM_TOOLTIP_G"));
    showBlue->set_tooltip_text  (M("HISTOGRAM_TOOLTIP_B"));
    showValue->set_tooltip_text (M("HISTOGRAM_TOOLTIP_L"));
    showChro->set_tooltip_text  (M("HISTOGRAM_TOOLTIP_CHRO"));
    showMode->set_tooltip_text  (M("HISTOGRAM_TOOLTIP_MODE"));
    scopeOptions->set_tooltip_text(M("HISTOGRAM_TOOLTIP_SHOW_OPTIONS"));
    scopeHistBtn->set_tooltip_text(M("HISTOGRAM_TOOLTIP_TYPE_HISTOGRAM"));
    scopeHistRawBtn->set_tooltip_text(M("HISTOGRAM_TOOLTIP_TYPE_HISTOGRAM_RAW"));
    scopeParadeBtn->set_tooltip_text(M("HISTOGRAM_TOOLTIP_TYPE_PARADE"));
    scopeWaveBtn->set_tooltip_text(M("HISTOGRAM_TOOLTIP_TYPE_WAVEFORM"));
    scopeVectHcBtn->set_tooltip_text(M("HISTOGRAM_TOOLTIP_TYPE_VECTORSCOPE_HC"));
    scopeVectHsBtn->set_tooltip_text(M("HISTOGRAM_TOOLTIP_TYPE_VECTORSCOPE_HS"));

    buttonGrid = Gtk::manage (new Gtk::Grid ());
    buttonGrid->set_orientation(Gtk::ORIENTATION_HORIZONTAL);
    persistentButtons = Gtk::manage(new Gtk::Box());
    persistentButtons->set_orientation(Gtk::ORIENTATION_VERTICAL);
    optionButtons = Gtk::manage(new Gtk::Box());
    optionButtons->set_orientation(Gtk::ORIENTATION_VERTICAL);

    showRed->set_active   (options.histogramRed);
    showGreen->set_active (options.histogramGreen);
    showBlue->set_active  (options.histogramBlue);
    showValue->set_active (options.histogramLuma);
    showChro->set_active  (options.histogramChroma);
    // no showMode->set_active(), as it's not a ToggleButton
    scopeOptions->set_active(options.histogramShowOptionButtons);
    showBAR->set_active   (options.histogramBar);

    showRed->set_image   (showRed->get_active()   ? *redImage   : *redImage_g);
    showGreen->set_image (showGreen->get_active() ? *greenImage : *greenImage_g);
    showBlue->set_image  (showBlue->get_active()  ? *blueImage  : *blueImage_g);
    showValue->set_image (showValue->get_active() ? *valueImage : *valueImage_g);
    showChro->set_image  (showChro->get_active()  ? *chroImage  : *chroImage_g);
    toggleButtonMode();
    scopeHistBtn->set_image(*histImage);
    scopeHistRawBtn->set_image(*histRawImage);
    scopeParadeBtn->set_image(*paradeImage);
    scopeWaveBtn->set_image(*waveImage);
    scopeVectHcBtn->set_image(*vectHcImage);
    scopeVectHsBtn->set_image(*vectHsImage);
    switch(options.histogramScopeType) {
        case ScopeType::HISTOGRAM:
            scopeHistBtn->set_active();
            break;
        case ScopeType::HISTOGRAM_RAW:
            scopeHistRawBtn->set_active();
            break;
        case ScopeType::PARADE:
            scopeParadeBtn->set_active();
            break;
        case ScopeType::WAVEFORM:
            scopeWaveBtn->set_active();
            break;
        case ScopeType::VECTORSCOPE_HS:
            scopeVectHsBtn->set_active();
            break;
        case ScopeType::VECTORSCOPE_HC:
            scopeVectHcBtn->set_active();
            break;
        case ScopeType::NONE:
            break;
    }
    scopeOptions->set_image(*Gtk::manage(new RTImage("histogram-ellipsis-small.png")));
    showBAR->set_image   (showBAR->get_active()   ? *barImage   : *barImage_g);

    setExpandAlignProperties(showRed  , false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_CENTER);
    setExpandAlignProperties(showGreen, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_CENTER);
    setExpandAlignProperties(showBlue , false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_CENTER);
    setExpandAlignProperties(showValue, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_CENTER);
    setExpandAlignProperties(showChro , false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_CENTER);
    setExpandAlignProperties(showMode , false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_CENTER);
    setExpandAlignProperties(scopeOptions, false, false, Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
    setExpandAlignProperties(showBAR  , false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_CENTER);
    setExpandAlignProperties(scopeOptions, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    setExpandAlignProperties(scopeHistBtn, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    setExpandAlignProperties(scopeHistRawBtn, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    setExpandAlignProperties(scopeParadeBtn, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    setExpandAlignProperties(scopeWaveBtn, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    setExpandAlignProperties(scopeVectHcBtn, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    setExpandAlignProperties(scopeVectHsBtn, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_FILL);
    setExpandAlignProperties(persistentButtons, false, true, Gtk::ALIGN_START, Gtk::ALIGN_FILL);
    setExpandAlignProperties(optionButtons, false, true, Gtk::ALIGN_START, Gtk::ALIGN_FILL);

    showRed->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::red_toggled), showRed );
    showGreen->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::green_toggled), showGreen );
    showBlue->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::blue_toggled), showBlue );
    showValue->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::value_toggled), showValue );
    showChro->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::chro_toggled), showChro );
    showMode->signal_released().connect( sigc::mem_fun(*this, &HistogramPanel::mode_released), showMode );
    scopeOptions->signal_toggled().connect(sigc::mem_fun(*this, &HistogramPanel::scopeOptionsToggled));
    showBAR->signal_toggled().connect( sigc::mem_fun(*this, &HistogramPanel::bar_toggled), showBAR );
    scopeHistBtn->signal_toggled().connect(sigc::bind(sigc::mem_fun(*this, &HistogramPanel::type_selected), scopeHistBtn));
    scopeHistRawBtn->signal_toggled().connect(sigc::bind(sigc::mem_fun(*this, &HistogramPanel::type_selected), scopeHistRawBtn));
    scopeParadeBtn->signal_toggled().connect(sigc::bind(sigc::mem_fun(*this, &HistogramPanel::type_selected), scopeParadeBtn));
    scopeWaveBtn->signal_toggled().connect(sigc::bind(sigc::mem_fun(*this, &HistogramPanel::type_selected), scopeWaveBtn));
    scopeVectHcBtn->signal_toggled().connect(sigc::bind(sigc::mem_fun(*this, &HistogramPanel::type_selected), scopeVectHcBtn));
    scopeVectHsBtn->signal_toggled().connect(sigc::bind(sigc::mem_fun(*this, &HistogramPanel::type_selected), scopeVectHsBtn));

    brightnessWidget = Gtk::manage(new Gtk::Scale(Gtk::ORIENTATION_VERTICAL));
    brightnessWidget->set_inverted();
    brightnessWidget->set_range(log(HistogramArea::MIN_BRIGHT), log(HistogramArea::MAX_BRIGHT));
    brightnessWidget->set_draw_value(false);
    brightnessWidget->signal_value_changed().connect(sigc::mem_fun(*this, &HistogramPanel::brightnessWidgetValueChanged));
    brightnessWidget->set_name("histScale");
    brightnessWidget->set_tooltip_text(M("HISTOGRAM_TOOLTIP_TRACE_BRIGHTNESS"));
    setExpandAlignProperties(brightnessWidget, true, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_START);

    optionButtons->add(*showRed);
    optionButtons->add(*showGreen);
    optionButtons->add(*showBlue);
    optionButtons->add(*showValue);
    optionButtons->add(*showChro);
    optionButtons->add(*showMode);
    optionButtons->add(*showBAR);
    optionButtons->add(*brightnessWidget);

    Gtk::VSeparator* separator = Gtk::manage(new Gtk::VSeparator());
    setExpandAlignProperties(separator, true, false, Gtk::ALIGN_FILL, Gtk::ALIGN_CENTER);
    persistentButtons->add(*scopeHistBtn);
    persistentButtons->add(*scopeHistRawBtn);
    persistentButtons->add(*scopeParadeBtn);
    persistentButtons->add(*scopeWaveBtn);
    persistentButtons->add(*scopeVectHsBtn);
    persistentButtons->add(*scopeVectHcBtn);
    persistentButtons->add(*separator);
    persistentButtons->add(*scopeOptions);

    // Put the button vbox next to the window's border to be less disturbing
    if (options.histogramPosition == Options::HISTOGRAM_POS_LEFT) {
        buttonGrid->add(*persistentButtons);
        buttonGrid->add(*optionButtons);

        add (*buttonGrid);
        add (*gfxGrid);
    } else {
        buttonGrid->add(*optionButtons);
        buttonGrid->add(*persistentButtons);

        add (*gfxGrid);
        add (*buttonGrid);
    }

    show_all ();
    optionButtons->set_no_show_all();
    optionButtons->set_visible(options.histogramShowOptionButtons);

    type_changed();
    updateHistAreaOptions();
    if (histogramRGBArea) {
        updateHistRGBAreaOptions();
    }

    brightness_changed_connection = histogramArea->getBrighnessChangedSignal().connect(sigc::mem_fun(*this, &HistogramPanel::brightnessUpdated));
    rconn = signal_size_allocate().connect( sigc::mem_fun(*this, &HistogramPanel::resized) );

    histogramArea->setBrightness(options.histogramTraceBrightness);
}

HistogramPanel::~HistogramPanel ()
{
    delete redImage;
    delete greenImage;
    delete blueImage;
    delete valueImage;
    delete chroImage;
    for (int i = 0; i < 3; ++i) {
        delete mode_images_[i];
    }
    delete barImage;

    delete redImage_g;
    delete greenImage_g;
    delete blueImage_g;
    delete valueImage_g;
    delete chroImage_g;
    delete barImage_g;

}

void HistogramPanel::showRGBBar()
{
    histogramRGBAreaHori->set_visible(
        histogramRGBArea == histogramRGBAreaHori.get() && showBAR->get_active());
    histogramRGBAreaVert->set_visible(
        histogramRGBArea == histogramRGBAreaVert.get() && showBAR->get_active());
    histogramRGBAreaHori->setShow(false);
    histogramRGBAreaVert->setShow(false);

    if (!histogramRGBArea) {
        return;
    }

    setHistRGBInvalid();
    histogramRGBArea->setShow(showBAR->get_active());
}

void HistogramPanel::resized (Gtk::Allocation& req)
{
    static int old_height = 0;
    static int old_width = 0;

    bool size_changed =
        old_height != req.get_height() || old_width != req.get_width();

    if (!histogramArea->updatePending() && size_changed) {
        histogramArea->updateBackBuffer ();
        histogramArea->queue_draw ();
    }

    // set histogramRGBArea invalid;
    if (histogramRGBArea && size_changed) {
        histogramRGBArea->updateBackBuffer(-1, -1, -1);
        histogramRGBArea->queue_draw ();
    }

    // Store current height of the histogram
    options.histogramHeight = get_height();

    old_height = req.get_height();
    old_width = req.get_width();
}

void HistogramPanel::red_toggled ()
{
    showRed->set_image(showRed->get_active() ? *redImage : *redImage_g);
    rgbv_toggled();
}
void HistogramPanel::green_toggled ()
{
    showGreen->set_image(showGreen->get_active() ? *greenImage : *greenImage_g);
    rgbv_toggled();
}
void HistogramPanel::blue_toggled ()
{
    showBlue->set_image(showBlue->get_active() ? *blueImage : *blueImage_g);
    rgbv_toggled();
}
void HistogramPanel::value_toggled ()
{
    removeIfThere(showValue, valueImage, false);
    removeIfThere(showValue, valueImage_g, false);
    showValue->set_image(showValue->get_active() ? *valueImage : *valueImage_g);
    rgbv_toggled();
}
void HistogramPanel::chro_toggled ()
{
    removeIfThere(showChro, chroImage, false);
    removeIfThere(showChro, chroImage_g, false);
    showChro->set_image(showChro->get_active() ? *chroImage : *chroImage_g);
    rgbv_toggled();
}

void HistogramPanel::mode_released ()
{
    options.histogramDrawMode = (options.histogramDrawMode + 1) % 3;
    toggleButtonMode();
    rgbv_toggled();
}

void HistogramPanel::brightnessWidgetValueChanged(void)
{
    ConnectionBlocker blocker(brightness_changed_connection);
    histogramArea->setBrightness(exp(brightnessWidget->get_value()));
    options.histogramTraceBrightness = histogramArea->getBrightness();
}

void HistogramPanel::brightnessUpdated(float brightness)
{
    brightnessWidget->set_value(log(brightness));
    options.histogramTraceBrightness = histogramArea->getBrightness();
}

void HistogramPanel::scopeOptionsToggled()
{
    options.histogramShowOptionButtons = scopeOptions->get_active();
    optionButtons->set_visible(scopeOptions->get_active());
}

void HistogramPanel::type_selected(Gtk::RadioButton* button)
{
    ScopeType new_type = ScopeType::NONE;

    if (button == scopeHistBtn) {
        new_type = ScopeType::HISTOGRAM;
    } else if (button == scopeHistRawBtn) {
        new_type = ScopeType::HISTOGRAM_RAW;
    } else if (button == scopeParadeBtn) {
        new_type = ScopeType::PARADE;
    } else if (button == scopeWaveBtn) {
        new_type = ScopeType::WAVEFORM;
    } else if (button == scopeVectHcBtn) {
        new_type = ScopeType::VECTORSCOPE_HC;
    } else if (button == scopeVectHsBtn) {
        new_type = ScopeType::VECTORSCOPE_HS;
    }

    if (new_type == options.histogramScopeType) {
        return;
    }

    options.histogramScopeType = new_type;

    type_changed();
    updateHistAreaOptions();
    if (histogramRGBArea) {
        updateHistRGBAreaOptions();
    }
    histogramArea->setDirty(true);
    histogramArea->queue_draw();
}

void HistogramPanel::type_changed()
{
    switch (options.histogramScopeType) {
        case ScopeType::HISTOGRAM:
            showRed->show();
            showGreen->show();
            showBlue->show();
            showValue->show();
            showChro->show();
            showMode->show();
            showBAR->show();
            showBAR->set_tooltip_text(M("HISTOGRAM_TOOLTIP_BAR"));
            brightnessWidget->hide();
            histogramRGBArea = histogramRGBAreaHori.get();
            break;
        case ScopeType::HISTOGRAM_RAW:
            showRed->show();
            showGreen->show();
            showBlue->show();
            showValue->hide();
            showChro->hide();
            showMode->show();
            showBAR->hide();
            brightnessWidget->hide();
            histogramRGBArea = nullptr;
            break;
        case ScopeType::PARADE:
        case ScopeType::WAVEFORM:
            showRed->show();
            showGreen->show();
            showBlue->show();
            showValue->show();
            showChro->hide();
            showMode->hide();
            showBAR->show();
            showBAR->set_tooltip_text(M("HISTOGRAM_TOOLTIP_BAR"));
            brightnessWidget->show();
            histogramRGBArea = histogramRGBAreaVert.get();
            break;
        case ScopeType::VECTORSCOPE_HC:
        case ScopeType::VECTORSCOPE_HS:
            showRed->hide();
            showGreen->hide();
            showBlue->hide();
            showValue->hide();
            showChro->hide();
            showMode->hide();
            showBAR->show();
            showBAR->set_tooltip_text(M("HISTOGRAM_TOOLTIP_CROSSHAIR"));
            brightnessWidget->show();
            histogramRGBArea = nullptr;
            break;
        case ScopeType::NONE:
            break;
    }

    if (panel_listener) {
        updateHistAreaOptions();
        panel_listener->scopeTypeChanged(options.histogramScopeType);
    }

    showRGBBar();
}

void HistogramPanel::bar_toggled ()
{
    showBAR->set_image(showBAR->get_active() ? *barImage : *barImage_g);
    rgbv_toggled();
    showRGBBar();
}

void HistogramPanel::rgbv_toggled ()
{
    // Update Display
    updateHistAreaOptions();
    histogramArea->updateBackBuffer ();
    histogramArea->queue_draw ();

    if (histogramRGBArea) {
        updateHistRGBAreaOptions();
        histogramRGBArea->updateBackBuffer(-1, -1, -1);
        histogramRGBArea->queue_draw ();
    }
}

void HistogramPanel::setHistRGBInvalid ()
{
    // do something to un-show vertical bars
    histogramRGBArea->updateBackBuffer(-1, -1, -1);
    histogramRGBArea->queue_draw ();
}

void HistogramPanel::pointerMoved (bool validPos, const Glib::ustring &profile, const Glib::ustring &profileW, int x, int y, int r, int g, int b, bool isRaw)
{
    bool update_hist_area = false;

    if (!validPos) {
        // do something to un-show vertical bars
        if (histogramRGBArea) {
            histogramRGBArea->updateBackBuffer(-1, -1, -1);
        }
        update_hist_area = histogramArea->updatePointer(-1, -1, -1);
    } else {
        // do something to show vertical bars
        if (histogramRGBArea) {
            histogramRGBArea->updateBackBuffer(r, g, b, profile, profileW);
        }
        update_hist_area = histogramArea->updatePointer(r, g, b, profile, profileW);
    }
    if (histogramRGBArea) {
        histogramRGBArea->queue_draw();
    }
    if (update_hist_area) {
        histogramArea->queue_draw();
    }
}

/*
 * Move the vertical button bar to the right side
 * only allowed values for align are Gtk::POS_LEFT and Gtk::POS_RIGHT
 */
void HistogramPanel::reorder (Gtk::PositionType align)
{
    if (align == Gtk::POS_LEFT) {
        gfxGrid->reference();
        removeIfThere(this, gfxGrid, false);
        add (*gfxGrid);
        gfxGrid->unreference();

        gfxGrid->remove(*histogramRGBAreaVert);
        gfxGrid->add(*histogramRGBAreaVert);

        optionButtons->reference();
        removeIfThere(buttonGrid, optionButtons, false);
        buttonGrid->add(*optionButtons);
        optionButtons->unreference();
    } else {
        buttonGrid->reference();
        removeIfThere(this, buttonGrid, false);
        add (*buttonGrid);
        buttonGrid->unreference();

        gfxGrid->remove(*histogramRGBAreaVert);
        gfxGrid->attach_next_to(*histogramRGBAreaVert, *histogramArea, Gtk::POS_LEFT, 1, 1);

        persistentButtons->reference();
        removeIfThere(buttonGrid, persistentButtons, false);
        buttonGrid->add(*persistentButtons);
        persistentButtons->unreference();
    }
}

// DrawModeListener interface:
void HistogramPanel::toggleButtonMode ()
{
    int m = LIM(options.histogramDrawMode, 0, 2);
    showMode->set_image(*mode_images_[m]);
    showMode->set_tooltip_text(Glib::ustring::compose(M("HISTOGRAM_TOOLTIP_MODE"), mode_tips_[m]));
}

void HistogramPanel::setPanelListener(HistogramPanelListener* listener)
{
    panel_listener = listener;

    if (listener) {
        listener->scopeTypeChanged(options.histogramScopeType);
    }
}

void HistogramPanel::updateHistAreaOptions()
{
    histogramArea->updateOptions(
        showRed->get_active(),
        showGreen->get_active(),
        showBlue->get_active(),
        showValue->get_active(),
        showChro->get_active(),
        options.histogramDrawMode,
        options.histogramScopeType,
        showBAR->get_active()
    );
}

void HistogramPanel::updateHistRGBAreaOptions()
{
    histogramRGBArea->updateOptions(
        showRed->get_active(),
        showGreen->get_active(),
        showBlue->get_active(),
        showValue->get_active(),
        showChro->get_active(),
        showBAR->get_active()
    );
}

//
//
//
// HistogramScaling

HistogramScaling::HistogramScaling():
    factor(10.0)
{
    factor = options.histogram_scaling_factor;
}


double HistogramScaling::log(double vsize, double val)
{
    //double factor = 10.0; // can be tuned if necessary - higher is flatter curve
    return vsize * std::log(factor / (factor + val)) / std::log(factor / (factor + vsize));
}

//
//
//
// HistogramRGBArea
HistogramRGBArea::HistogramRGBArea () :
    val(0), r(0), g(0), b(0), valid(false),
    needRed(options.histogramRed), needGreen(options.histogramGreen), needBlue(options.histogramBlue),
    needLuma(options.histogramLuma), needChroma(options.histogramChroma),
    showMode(options.histogramBar), barDisplayed(options.histogramBar), parent(nullptr)
{
    get_style_context()->add_class("drawingarea");
    set_name("HistogramRGBArea");

    harih = new HistogramRGBAreaIdleHelper;
    harih->harea = this;
    harih->destroyed = false;
    harih->pending = 0;
}

HistogramRGBArea::~HistogramRGBArea ()
{
    idle_register.destroy();

    if (harih->pending) {
        harih->destroyed = true;
    } else {
        delete harih;
    }
}


void HistogramRGBArea::getPreferredThickness(int& min_thickness, int& natural_thickness) const
{
    int minimumLength = 0;
    int naturalLength = 0;
    getPreferredLength(minimumLength, naturalLength);
    getPreferredThicknessForLength(minimumLength, min_thickness, natural_thickness);
}

void HistogramRGBArea::getPreferredLength(int& min_length, int& natural_length) const
{
    int s = RTScalable::getScale();
    min_length = 60 * s;
    natural_length = 200 * s;
}

void HistogramRGBArea::getPreferredThicknessForLength(int length, int& min_thickness, int& natural_thickness) const
{
    int bThickness = length / 30;

    int s = RTScalable::getScale();

    if (bThickness > (10 * s)) {
        bThickness = 10 * s;
    } else if (bThickness < (5 * s)) {
        bThickness = 5 * s;
    }

    min_thickness = bThickness;
    natural_thickness = bThickness;
}

// unused?
void HistogramRGBArea::getPreferredLengthForThickness(int thickness, int& min_length, int& natural_length) const
{
    getPreferredLength(min_length, natural_length);
}

bool HistogramRGBArea::getShow()
{
    return(showMode);
}

void HistogramRGBArea::setShow(bool show)
{
    showMode = show;
}

void HistogramRGBArea::updateBackBuffer (int r, int g, int b, const Glib::ustring &profile, const Glib::ustring &profileW)
{
    if (!get_realized () || !showMode || !(
        options.histogramScopeType == ScopeType::HISTOGRAM
        || options.histogramScopeType == ScopeType::PARADE
        || options.histogramScopeType == ScopeType::WAVEFORM
    )) {
        return;
    }

    // Mostly not necessary, but should be in some case
    GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected

    Glib::RefPtr<Gdk::Window> window = get_window();
    int winx, winy, winw, winh;
    window->get_geometry(winx, winy, winw, winh);

    double s = RTScalable::getScale();

    // This will create or update the size of the BackBuffer::surface
    setDrawRectangle(Cairo::FORMAT_ARGB32, 0, 0, winw, winh, true);

    if (surface)  {
        Cairo::RefPtr<Cairo::Context> cc = Cairo::Context::create(surface);

        cc->set_source_rgba (0., 0., 0., 0.);
        cc->set_operator (Cairo::OPERATOR_CLEAR);
        cc->paint ();
        cc->set_operator (Cairo::OPERATOR_OVER);

        cc->set_antialias(Cairo::ANTIALIAS_NONE);
        cc->set_line_width (1.0 * s);

        if ( r != -1 && g != -1 && b != -1 ) {
            if (needRed) {
                // Red
                cc->set_source_rgb(rgb_R[0], rgb_R[1], rgb_R[2]);
                drawBar(cc, r, 255.0, winw, winh, s);
            }

            if (needGreen) {
                // Green
                cc->set_source_rgb(rgb_G[0], rgb_G[1], rgb_G[2]);
                drawBar(cc, g, 255.0, winw, winh, s);
            }

            if (needBlue) {
                // Blue
                cc->set_source_rgb(rgb_B[0], rgb_B[1], rgb_B[2]);
                drawBar(cc, b, 255.0, winw, winh, s);
            }

            if(
                (needLuma || needChroma)
                && (options.histogramScopeType == ScopeType::HISTOGRAM
                    || options.histogramScopeType == ScopeType::PARADE
                    || options.histogramScopeType == ScopeType::WAVEFORM)
            ) {
                float Lab_L, Lab_a, Lab_b;
                rtengine::Color::rgb2lab01(profile, profileW, r / 255.f, g / 255.f, b / 255.f, Lab_L, Lab_a, Lab_b, options.rtSettings.HistogramWorking);

                if (needLuma) {
                    // Luma
                    cc->set_source_rgb(1.0, 1.0, 1.0);
                    drawBar(cc, Lab_L, 100.0, winw, winh, s);
                }

                if (needChroma && options.histogramScopeType == ScopeType::HISTOGRAM) {
                    // Chroma
                    double chromaval = sqrt(Lab_a * Lab_a + Lab_b * Lab_b) / 1.8;
                    cc->set_source_rgb(0.9, 0.9, 0.0);
                    drawBar(cc, chromaval, 100.0, winw, winh, s);
                }
            }
        }
    }

    setDirty(false);
}

void HistogramRGBArea::update (int valh, int rh, int  gh, int bh)
{

    if (valh) {
        val = valh;
        r = rh;
        g = gh;
        b = bh;
        valid = true;
    } else {
        valid = false;
    }

    harih->pending++;

    idle_register.add(
        [this]() -> bool
        {
            if (harih->destroyed) {
                if (harih->pending == 1) {
                    delete harih;
                } else {
                    --harih->pending;
                }

                return false;
            }

            harih->harea->updateBackBuffer(-1, -1, -1);
            harih->harea->queue_draw ();

            --harih->pending;

            return false;
        }
    );
}

void HistogramRGBArea::updateOptions (bool r, bool g, bool b, bool l, bool c, bool bar)
{

    options.histogramRed    = needRed    = r;
    options.histogramGreen  = needGreen  = g;
    options.histogramBlue   = needBlue   = b;
    options.histogramLuma   = needLuma   = l;
    options.histogramChroma = needChroma = c;
    options.histogramBar    = showMode   = bar;

}

void HistogramRGBArea::on_realize ()
{

    Gtk::DrawingArea::on_realize();
    add_events(Gdk::BUTTON_PRESS_MASK);
}

bool HistogramRGBArea::on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr)
{

    const Glib::RefPtr<Gtk::StyleContext> style = get_style_context();
    style->render_background(cr, 0, 0, get_width(), get_height());

    // on_realize & updateBackBuffer have to be called before
    if (surface) {
        if (isDirty()) { // not sure this could happen...
            updateBackBuffer(-1, -1, -1);
        }

        copySurface(cr, NULL);
    }

    style->render_frame (cr, 0, 0, get_width(), get_height());

    return true;
}

bool HistogramRGBArea::on_button_press_event (GdkEventButton* event)
{

    if (event->type == GDK_2BUTTON_PRESS && event->button == 1) {
        // do something?
    }

    return true;
}

void HistogramRGBArea::factorChanged (double newFactor)
{
    factor = newFactor;
    options.histogram_scaling_factor = factor;
}

void HistogramRGBAreaHori::drawBar(Cairo::RefPtr<Cairo::Context> cc, double value, double max_value, int winw, int winh, double scale)
{
    double pos;
    if (options.histogramDrawMode == 0) {
        pos = padding + value * (winw - padding * 2.0) / max_value + 0.5 * scale;
    } else {
        pos = padding + HistogramScaling::log (max_value, value) * (winw - padding * 2.0) / max_value + 0.5 * scale;
    }
    cc->move_to(pos, 0.0);
    cc->line_to(pos, winh - 0.0);
    cc->stroke();
}

Gtk::SizeRequestMode HistogramRGBAreaHori::get_request_mode_vfunc () const
{
    return Gtk::SIZE_REQUEST_HEIGHT_FOR_WIDTH;
}

void HistogramRGBAreaHori::get_preferred_height_vfunc (int &minimum_height, int &natural_height) const
{
    getPreferredThickness(minimum_height, natural_height);
}

void HistogramRGBAreaHori::get_preferred_width_vfunc (int &minimum_width, int &natural_width) const
{
    getPreferredLength(minimum_width, natural_width);
}

void HistogramRGBAreaHori::get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const
{
    getPreferredThicknessForLength(width, minimum_height, natural_height);
}

void HistogramRGBAreaHori::get_preferred_width_for_height_vfunc (int height, int &minimum_width, int &natural_width) const
{
    getPreferredLengthForThickness(height, minimum_width, natural_width);
}

void HistogramRGBAreaVert::drawBar(Cairo::RefPtr<Cairo::Context> cc, double value, double max_value, int winw, int winh, double scale)
{
    double pos;
    if (options.histogramDrawMode == 0 || options.histogramScopeType == ScopeType::PARADE || options.histogramScopeType == ScopeType::WAVEFORM) {
        pos = padding + value * (winh - padding * 2.0 - 1) / max_value + 0.5 * scale;
    } else {
        pos = padding + HistogramScaling::log (max_value, value) * (winh - padding * 2.0) / max_value + 0.5 * scale;
    }
    cc->move_to(0.0, winh - pos);
    cc->line_to(winw, winh - pos);
    cc->stroke();
}

Gtk::SizeRequestMode HistogramRGBAreaVert::get_request_mode_vfunc () const
{
    return Gtk::SIZE_REQUEST_WIDTH_FOR_HEIGHT;
}

void HistogramRGBAreaVert::get_preferred_height_vfunc (int &minimum_height, int &natural_height) const
{
    getPreferredLength(minimum_height, natural_height);
}

void HistogramRGBAreaVert::get_preferred_width_vfunc (int &minimum_width, int &natural_width) const
{
    minimum_width = 10 * RTScalable::getScale();
    natural_width = minimum_width;
}

void HistogramRGBAreaVert::get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const
{
    getPreferredLengthForThickness(width, minimum_height, natural_height);
}

void HistogramRGBAreaVert::get_preferred_width_for_height_vfunc (int height, int &minimum_width, int &natural_width) const
{
    get_preferred_width_vfunc(minimum_width, natural_width);
}

//
//
//
// HistogramArea
HistogramArea::HistogramArea(DrawModeListener *fml, bool is_main):
    vectorscope_scale(0),
    vect_hc(0, 0), vect_hs(0, 0),
    vect_hc_buffer_dirty(true), vect_hs_buffer_dirty(true),
    waveform_scale(0),
    rwave(0, 0), gwave(0, 0), bwave(0, 0), lwave(0, 0),
    parade_buffer_r_dirty(true), parade_buffer_g_dirty(true), parade_buffer_b_dirty(true),
    wave_buffer_dirty(true), wave_buffer_luma_dirty(true),
    valid(false), drawMode(options.histogramDrawMode), myDrawModeListener(fml),
    scopeType(options.histogramScopeType),
    oldwidth(-1), oldheight(-1),
    trace_brightness(1.0),
    needRed(options.histogramRed), needGreen(options.histogramGreen), needBlue(options.histogramBlue),
    needLuma(options.histogramLuma), needChroma(options.histogramChroma),
    isPressed(false), movingPosition(0.0),
    needPointer(options.histogramBar),
    pointer_red(-1), pointer_green(-1), pointer_blue(-1),
    pointer_a(0), pointer_b(0),
    is_main_(is_main)
{

    rhist(256);
    ghist(256);
    bhist(256);
    lhist(256);
    chist(256);

    get_style_context()->add_class("drawingarea");
    set_name("HistogramArea");

    haih = new HistogramAreaIdleHelper;
    haih->harea = this;
    haih->destroyed = false;
    haih->pending = 0;
}

HistogramArea::~HistogramArea ()
{
    idle_register.destroy();

    if (haih->pending) {
        haih->destroyed = true;
    } else {
        delete haih;
    }
}

Gtk::SizeRequestMode HistogramArea::get_request_mode_vfunc () const
{
    return Gtk::SIZE_REQUEST_CONSTANT_SIZE;
}

void HistogramArea::get_preferred_height_vfunc (int &minimum_height, int &natural_height) const
{
    int s = RTScalable::getScale();
    minimum_height = 100 * s;
    natural_height = 200 * s;
}

void HistogramArea::get_preferred_width_vfunc (int &minimum_width, int &natural_width) const
{

    int s = RTScalable::getScale();
    minimum_width = 200 * s;
    natural_width = 400 * s;
}

void HistogramArea::get_preferred_height_for_width_vfunc (int width, int &minimum_height, int &natural_height) const
{

    minimum_height = 0;
    natural_height = 0;
}

void HistogramArea::get_preferred_width_for_height_vfunc (int height, int &minimum_width, int &natural_width) const
{
    get_preferred_width_vfunc (minimum_width, natural_width);
}

void HistogramArea::updateOptions (bool r, bool g, bool b, bool l, bool c, int mode, ScopeType type, bool pointer)
{
    wave_buffer_dirty = wave_buffer_dirty || needRed != r || needGreen != g || needBlue != b;

    needRed = r;
    needGreen = g;
    needBlue = b;
    needLuma = l;
    needChroma = c;
    drawMode = mode;
    scopeType = type;
    needPointer = pointer;

    if (is_main_) {
        options.histogramRed = needRed;
        options.histogramGreen = needGreen;
        options.histogramBlue = needBlue;
        options.histogramLuma = needLuma;
        options.histogramChroma = needChroma;
        options.histogramDrawMode = drawMode;
        options.histogramScopeType = scopeType;
        options.histogramBar = needPointer;
    }
}

bool HistogramArea::updatePending(void)
{
    return haih->pending > 0 && !haih->destroyed;
}

void HistogramArea::update(
    const LUTu& histRed,
    const LUTu& histGreen,
    const LUTu& histBlue,
    const LUTu& histLuma,
    const LUTu& histChroma,
    const LUTu& histRedRaw,
    const LUTu& histGreenRaw,
    const LUTu& histBlueRaw,
    int vectorscopeScale,
    const array2D<int>& vectorscopeHC,
    const array2D<int>& vectorscopeHS,
    int waveformScale,
    const array2D<int>& waveformRed,
    const array2D<int>& waveformGreen,
    const array2D<int>& waveformBlue,
    const array2D<int>& waveformLuma
)
{
    if (histRed) {
        switch (scopeType) {
            case ScopeType::HISTOGRAM:
                rhist = histRed;
                ghist = histGreen;
                bhist = histBlue;
                lhist = histLuma;
                chist = histChroma;
                break;
            case ScopeType::HISTOGRAM_RAW:
                rhistRaw = histRedRaw;
                ghistRaw = histGreenRaw;
                bhistRaw = histBlueRaw;
                break;
            case ScopeType::PARADE:
            case ScopeType::WAVEFORM:
                waveform_scale = waveformScale;
                set_arr(rwave, waveformRed);
                set_arr(gwave, waveformGreen);
                set_arr(bwave, waveformBlue);
                set_arr(lwave, waveformLuma);
                parade_buffer_r_dirty = parade_buffer_g_dirty = parade_buffer_b_dirty = wave_buffer_dirty = wave_buffer_luma_dirty = true;
                break;
            case ScopeType::VECTORSCOPE_HS:
                vectorscope_scale = vectorscopeScale;
                set_arr(vect_hs, vectorscopeHS);
                vect_hs_buffer_dirty = true;
                break;
            case ScopeType::VECTORSCOPE_HC:
                vectorscope_scale = vectorscopeScale;
                set_arr(vect_hc, vectorscopeHC);
                vect_hc_buffer_dirty = true;
                break;
            case ScopeType::NONE:
                break;
        }
        valid = true;
    } else {
        valid = false;
    }

    haih->pending++;

    // Can be done outside of the GUI thread
    idle_register.add(
        [this]() -> bool
        {
            if (haih->destroyed) {
                if (haih->pending == 1) {
                    delete haih;
                } else {
                    --haih->pending;
                }

                return false;
            }

            haih->harea->setDirty(true);
            haih->harea->updateBackBuffer();
            haih->harea->queue_draw();

            --haih->pending;

            return false;
        }
    );
}

void HistogramArea::updateBackBuffer(int custom_w, int custom_h)
{
    int winw = custom_w;
    int winh = custom_h;
    if (winw < 0) {
        if (!get_realized ()) {
            return;
        }

        Glib::RefPtr<Gdk::Window> window = get_window();
        int winx, winy;
        window->get_geometry(winx, winy, winw, winh);
    }

    // This will create or update the size of the BackBuffer::surface
    setDrawRectangle(Cairo::FORMAT_ARGB32, 0, 0, winw, winh, true);

    Cairo::RefPtr<Cairo::Context> cr = Cairo::Context::create(surface);
    const Glib::RefPtr<Gtk::StyleContext> style = get_style_context();

    double s = RTScalable::getScale();

    // Setup drawing
    cr->set_source_rgba (0., 0., 0., 0.);
    cr->set_operator (Cairo::OPERATOR_CLEAR);
    cr->paint ();

    cr->set_operator (Cairo::OPERATOR_SOURCE);

    // Prepare drawing gridlines first
    cr->set_source_rgba (1., 1., 1., 0.25);
    cr->set_line_width (1.0 * s);
    cr->set_antialias(Cairo::ANTIALIAS_NONE);
    cr->set_line_join(Cairo::LINE_JOIN_MITER);
    cr->set_line_cap(Cairo::LINE_CAP_BUTT);
    std::valarray<double> ch_ds (1);
    ch_ds[0] = 4;
    cr->set_dash (ch_ds, 0);

    bool rawMode = (scopeType == ScopeType::HISTOGRAM_RAW);
    if (rawMode) {
        updateRaw(cr);
    } else {
        updateNonRaw(cr);
    }

    // Draw the frame's border
    style->render_frame(cr, 0, 0, surface->get_width(), surface->get_height());

    oldwidth = w;
    oldheight = h;

    setDirty(false);
}


void HistogramArea::updateNonRaw(Cairo::RefPtr<Cairo::Context> cr)
{
    set_has_tooltip(false);
    
    const Glib::RefPtr<Gtk::StyleContext> style = get_style_context();

    double s = RTScalable::getScale();

    // determine the number of h-gridlines based on current h
    int nrOfHGridPartitions = (int)rtengine::min (16.0, pow (2.0, floor ((h - 100) / 250) + 2));
    int nrOfVGridPartitions = 8; // always show 8 stops (lines at 1,3,7,15,31,63,127)

    // draw vertical gridlines
    if (scopeType == ScopeType::HISTOGRAM) {
        if (drawMode == 0) {
            for (int i = 1; i < nrOfVGridPartitions; i++) {
                cr->move_to((pow(2.0,i) - 1) / 255.0 * w + 0.5, 0.);
                cr->line_to((pow(2.0,i) - 1) / 255.0 * w + 0.5, h);
                cr->stroke();
            }
        } else {
            for (int i = 1; i < nrOfVGridPartitions; i++) {
                cr->move_to(HistogramScaling::log(255, pow(2.0,i) - 1) / 255.0 * w + 0.5, 0.);
                cr->line_to(HistogramScaling::log(255, pow(2.0,i) - 1) / 255.0 * w + 0.5, h);
                cr->stroke();
            }
        }
    }

    // draw horizontal gridlines
    if (scopeType == ScopeType::PARADE || scopeType == ScopeType::WAVEFORM) {
        for (int i = 0; i <= nrOfVGridPartitions; i++) {
            const double ypos = h - padding - (pow(2.0,i) - 1) * (h - 2 * padding - 1) / 255.0;
            cr->move_to(0, ypos);
            cr->line_to(w, ypos);
            cr->stroke();
        }
    } else if (scopeType == ScopeType::VECTORSCOPE_HC || scopeType == ScopeType::VECTORSCOPE_HS) {
        // Vectorscope has no gridlines.
    } else if (drawMode != 2) {
        for (int i = 1; i < nrOfHGridPartitions; i++) {
            cr->move_to (0., i * (double)h / nrOfHGridPartitions + 0.5);
            cr->line_to (w, i * (double)h / nrOfHGridPartitions + 0.5);
            cr->stroke ();
        }
    } else {
        for (int i = 1; i < nrOfHGridPartitions; i++) {
            cr->move_to (0., h - HistogramScaling::log (h, i * (double)h / nrOfHGridPartitions) + 0.5*s);
            cr->line_to (w, h - HistogramScaling::log (h, i * (double)h / nrOfHGridPartitions) + 0.5*s);
            cr->stroke ();
        }
    }

    cr->unset_dash();

    if (valid && scopeType == ScopeType::HISTOGRAM) {
        // For RAW mode use the other hists
        LUTu& rh = rhist;
        LUTu& gh = ghist;
        LUTu& bh = bhist;

        // make copies of LUT for faster access
        unsigned int lhisttemp[256] ALIGNED16 {0}, chisttemp[256] ALIGNED16 {0}, rhtemp[256] ALIGNED16 {0}, ghtemp[256] ALIGNED16 {0}, bhtemp[256] ALIGNED16 {0};
        const int scale = 1;

        for(int i = 0; i < 256; i++) {
            if(needLuma) {
                lhisttemp[i] = lhist[i];
            }

            if(needChroma) {
                chisttemp[i] = chist[i];
            }

            if(needRed) {
                rhtemp[i] = rh[i];
            }

            if(needGreen) {
                ghtemp[i] = gh[i];
            }

            if(needBlue) {
                bhtemp[i] = bh[i];
            }
        }

        // Compute the highest point of the histogram for scaling
        // Values at far left and right end (0 and 255) are handled differently

        unsigned int histheight = 0;

        for (int i = 1; i < 255; i++) {
            if (needLuma && lhisttemp[i] > histheight) {
                histheight = lhisttemp[i];
            }

            if (needChroma && chisttemp[i] > histheight) {
                histheight = chisttemp[i];
            }

            if (needRed && rhtemp[i] > histheight) {
                histheight = rhtemp[i];
            }

            if (needGreen && ghtemp[i] > histheight) {
                histheight = ghtemp[i];
            }

            if (needBlue && bhtemp[i] > histheight) {
                histheight = bhtemp[i];
            }
        }

        int realhistheight = histheight;

        if (realhistheight < h - 2) {
            realhistheight = h - 2;
        }

        cr->set_antialias (Cairo::ANTIALIAS_SUBPIXEL);
        cr->set_line_width (1.0 * s);
        cr->set_operator (Cairo::OPERATOR_OVER);

        int ui = 0, oi = 0;

        const bool rawMode = (scopeType == ScopeType::HISTOGRAM_RAW);
        if (needLuma && !rawMode) {
            drawCurve(cr, lhist, realhistheight, w, h);
            cr->set_source_rgba (0.65, 0.65, 0.65, 0.65);
            cr->fill ();
            drawMarks(cr, lhist, scale, w, ui, oi);
        }

        if (needChroma && !rawMode) {
            drawCurve(cr, chist, realhistheight, w, h);
            cr->set_source_rgb (0.9, 0.9, 0.);
            cr->stroke ();
            drawMarks(cr, chist, scale, w, ui, oi);
        }

        if (needRed) {
            drawCurve(cr, rh, realhistheight, w, h);
            cr->set_source_rgb(rgb_R[0], rgb_R[1], rgb_R[2]);
            cr->stroke ();
            drawMarks(cr, rh, scale, w, ui, oi);
        }

        if (needGreen) {
            drawCurve(cr, gh, realhistheight, w, h);
            cr->set_source_rgb(rgb_G[0], rgb_G[1], rgb_G[2]);
            cr->stroke ();
            drawMarks(cr, gh, scale, w, ui, oi);
        }

        if (needBlue) {
            drawCurve(cr, bh, realhistheight, w, h);
            cr->set_source_rgb(rgb_B[0], rgb_B[1], rgb_B[2]);
            cr->stroke ();
            drawMarks(cr, bh, scale, w, ui, oi);
        }

    } else if (scopeType == ScopeType::PARADE && rwave.width() > 0) {
        drawParade(cr, w, h);
    } else if (scopeType == ScopeType::WAVEFORM && rwave.width() > 0) {
        drawWaveform(cr, w, h);
    } else if (scopeType == ScopeType::VECTORSCOPE_HC || scopeType == ScopeType::VECTORSCOPE_HS) {
        drawVectorscope(cr, w, h);
    }
}


namespace {

class RawIdxHelper {
public:
    RawIdxHelper(bool logscale, unsigned int ub, unsigned int width):
        logscale_(logscale),
        ub_(ub)
    {
        incr_ = std::max(ub_ / width, 1u);
    }

    void operator()(unsigned int &cur, unsigned int &i) const
    {
        unsigned int j = i;
        i = cur;
        unsigned int next = cur;
        if (logscale_) {
            next = std::max(static_cast<unsigned int>(j * 1.3f), cur + 1u);
        } else {
            next = cur + incr_;
        }
        if (ub_ - next < next - cur) {
            next = ub_ + 1;
        }
        if (cur < ub_ && next > ub_) {
            cur = ub_;
        } else {
            cur = std::min(next, ub_+1);
        }
    }

private:
    bool logscale_;
    unsigned int ub_;
    unsigned int incr_;
};

} // namespace

bool HistogramArea::updatePointer(int r, int g, int b, const Glib::ustring &profile, const Glib::ustring &profileW)
{
    if (!needPointer || !(scopeType == ScopeType::VECTORSCOPE_HC || scopeType == ScopeType::VECTORSCOPE_HS)) {
        return false;
    }
    if (pointer_red == r && pointer_green == g && pointer_blue == b) {
        return false;
    }

    float L;
    pointer_red = r;
    pointer_green = g;
    pointer_blue = b;
    Color::rgb2lab01(profile, profileW, r / 255.f, g / 255.f, b / 255.f, L, pointer_a, pointer_b, options.rtSettings.HistogramWorking);
    updateBackBuffer();
    return true;
}

void HistogramArea::updateRaw(Cairo::RefPtr<Cairo::Context> cr)
{
    const Glib::RefPtr<Gtk::StyleContext> style = get_style_context();
    double s = RTScalable::getScale();

    // determine the number of h-gridlines based on current h
    int nrOfHGridPartitions = (int)rtengine::min (16.0, pow (2.0, floor ((h - 100) / 250) + 2));
    double sz = 2;
    if (valid) {
        sz = rtengine::max(rhistRaw.getUpperBound(), ghistRaw.getUpperBound(), bhistRaw.getUpperBound());
    }

    int minval[3];
    int mincount[3];
    int maxval[3];
    int maxcount[3];
    int peakval[3];
    int peakcount[3];
    int totcount[3];
    int totsum[3];
    const char *ch[3] = {"#FF0000", "#00FF00", "#0000FF"};

    // draw vertical gridlines
    const double logmax = std::log2(sz);
    const bool logscale = drawMode > 0;
    
    for (double i = sz / 2.0; ; i /= 2.0) {
        double x = i / sz;
        if (logscale) {
            x = logmax + std::log2(x);
            if (x <= 0) {
                break;
            }
            x /= logmax;
        }
        x *= w;
        cr->move_to(x, 0.);
        cr->line_to(x, h);
        cr->stroke();
        if (i <= 1.0) {
            break;
        }
    }

    // draw horizontal gridlines
    if (drawMode == 2) {
        for (int i = 1; i < nrOfHGridPartitions; i++) {
            double y = double(i) / nrOfHGridPartitions;
            y = rtengine::log2lin(y, 10.0);
            cr->move_to(0., y * h);
            cr->line_to(w, y * h);
            cr->stroke();
        }
    } else {
        for (int i = 1; i < nrOfHGridPartitions; i++) {
            cr->move_to(0., i * (double)h / nrOfHGridPartitions + 0.5);
            cr->line_to(w, i * (double)h / nrOfHGridPartitions + 0.5);
            cr->stroke();
        }
    }

    cr->unset_dash();

    if (valid) {
        // For RAW mode use the other hists
        LUTu &rh = rhistRaw;
        LUTu &gh = ghistRaw;
        LUTu &bh = bhistRaw;

        auto ub = rtengine::max(rh.getUpperBound(), gh.getUpperBound(), bh.getUpperBound());

        const auto bin =
            [](LUTu &data, int i, int next) -> double
            {
                double val = 0.0;
                for (int j = i; j < next; ++j) {
                    val += data[j];
                }
                val /= (next - i);
                return val;
            };
        
        //const int delta = std::max(int(ub / std::max(w / 2, 1)), 1);
        //const float delta = 1.05f;
        RawIdxHelper next_raw_idx(logscale, ub, w);
        unsigned int next = 1;
        unsigned int histheight = 0;

        const unsigned int off = drawMode == 2 ? 0 : 1;
        for (unsigned int i = off; i <= ub-off;) {
            if (i < rh.getSize()) {
                double val = bin(rh, i, next);
                if (histheight < val) {
                    histheight = val;
                }
            }
            if (i < gh.getSize()) {
                double val = bin(gh, i, next);
                if (histheight < val) {
                    histheight = val;
                }
            }
            if (i < bh.getSize()) {
                double val = bin(bh, i, next);
                if (histheight < val) {
                    histheight = val;
                }
            }

            next_raw_idx(next, i);
        }

        int realhistheight = histheight;

        if (realhistheight < h - 2) {
            realhistheight = h - 2;
        }

        cr->set_antialias (Cairo::ANTIALIAS_SUBPIXEL);
        cr->set_line_width (1.0 * s);
        cr->set_operator (Cairo::OPERATOR_OVER);

        int ui = 0, oi = 0;

        if (needRed) {
            drawRawCurve(cr, rh, ub, realhistheight, w, h);
            cr->set_source_rgb(rgb_R[0], rgb_R[1], rgb_R[2]);
            cr->stroke();
            drawMarks(cr, rh, 1.0, w, ui, oi);
        }

        if (needGreen) {
            drawRawCurve(cr, gh, ub, realhistheight, w, h);
            cr->set_source_rgb(rgb_G[0], rgb_G[1], rgb_G[2]);
            cr->stroke();
            drawMarks(cr, gh, 1.0, w, ui, oi);
        }

        if (needBlue) {
            drawRawCurve(cr, bh, ub, realhistheight, w, h);
            cr->set_source_rgb(rgb_B[0], rgb_B[1], rgb_B[2]);
            cr->stroke();
            drawMarks(cr, bh, 1.0, w, ui, oi);
        }

        // update stats
        for (int c = 0; c < 3; ++c) {
            minval[c] = -1;
            mincount[c] = 0;
            maxval[c] = -1;
            maxcount[c] = 0;
            peakval[c] = 0;
            peakcount[c] = 0;
            totcount[c] = 0;
            totsum[c] = 0;
        }

        const auto update =
            [&](int c, LUTu &data, unsigned int i) -> void
            {
                if (i < data.getSize()) {
                    int cnt = data[i];
                    if (cnt > 0) {
                        if (minval[c] < 0) {
                            minval[c] = i;
                            mincount[c] = cnt;
                        }
                        maxval[c] = i;
                        maxcount[c] = cnt;
                        if (peakcount[c] < cnt) {
                            peakval[c] = i;
                            peakcount[c] = cnt;
                        }
                        ++totcount[c];
                    }
                    totsum[c] += cnt;
                }
            };
        
        for (unsigned int i = 0; i <= ub; ++i) {
            update(0, rh, i);
            update(1, gh, i);
            update(2, bh, i);
        }
    }


    if (valid) {
        const auto pct =
            [](int n, int d) -> Glib::ustring
            {
                double p = double(int(double(n)/double(d) * 10000)) / 100.0;
                if (p == 0 && n > 0) {
                    return "&lt;0.01%";
                } else {
                    return Glib::ustring::compose("%1%%", p);
                }
            };
        
        set_has_tooltip(true);
        Glib::ustring tip = "";
        for (int i = 0; i < 3; ++i) {
            double ev = minval[i] >= 0 ? std::log2(maxval[i] - minval[i]) : 0;
            ev = double(int(ev * 100)) / 100.0;
            Glib::ustring m = Glib::ustring::compose(M("HISTOGRAM_RAW_STATS_TOOLTIP"), minval[i], pct(mincount[i], totsum[i]), maxval[i], pct(maxcount[i], totsum[i]), peakval[i], pct(peakcount[i], totsum[i]), ev, totcount[i], totsum[i]);
            Glib::ustring s = Glib::ustring::compose("<span font_family=\"Arial\" size=\"larger\" foreground=\"%1\">&#9632;</span> %2", ch[i], m);
            if (i > 0) {
                tip += "\n";
            }
            tip += s;
        }
        set_tooltip_markup(tip);
    } else {
        set_has_tooltip(false);
    }
}


void HistogramArea::drawRawCurve(Cairo::RefPtr<Cairo::Context> &cr,
                                 LUTu &data, unsigned int ub, double scale,
                                 int hsize, int vsize)
{
    double s = RTScalable::getScale();

    cr->set_line_width(s);
    cr->move_to(0, vsize - 1);
    scale = scale <= 0.f ? 0.001f : scale; // avoid division by zero and negative values

    unsigned int next = 1;
    const bool logscale = drawMode > 0;
    const double logmax = std::log2(ub);
    //const double ybase = std::pow(10, std::max(std::floor(std::log(scale) / std::log(10))-1, 1.0));
    const double ylogmax = std::log10(scale);

    RawIdxHelper next_raw_idx(logscale, data.getUpperBound(), hsize);

    for (unsigned int i = logscale ? 1 : 0; i < data.getSize(); ) {
        double val = 0.0;
        for (unsigned int j = i; j < next; ++j) {
            val += data[j];
        }
        val /= (next - i);
        
        val = std::min(val / scale, 1.0);

        if (drawMode == 2 && val > 0) { // scale y for log-scale
            //val = rtengine::lin2log(val, ybase);
            val = std::max((ylogmax + std::log10(val)) / ylogmax, 0.0);
        }

        double iscaled = std::min(double(i) / ub, 1.0);
        if (logscale) { // scale x for log-scale
            iscaled = logmax + std::log2(iscaled);
            if (iscaled < 0) {
                next_raw_idx(next, i);
                continue;
            }
            iscaled /= logmax;
        }

        double posX = iscaled * hsize;
        double posY = vsize - 2 + val * (4 - vsize);

        cr->line_to(posX, posY);

        next_raw_idx(next, i);
    }

    cr->line_to(hsize - 1, vsize - 1);
}


void HistogramArea::on_realize ()
{

    Gtk::DrawingArea::on_realize();
    add_events(Gdk::POINTER_MOTION_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK);
}

void HistogramArea::drawCurve(Cairo::RefPtr<Cairo::Context> &cr,
                              const LUTu &data, double scale, int hsize, int vsize)
{
    double s = RTScalable::getScale();

    cr->set_line_width(s);
    cr->move_to (0, vsize - 1);
    scale = scale <= 0.f ? 0.001f : scale; // avoid division by zero and negative values

    for (int i = 0; i < 256; i++) {
        double val = data[i] * (double)vsize / scale;

        if (drawMode == 2) { // scale y for log-scale
            val = HistogramScaling::log ((double)vsize, val);
        }

        double iscaled = i;
        if (drawMode > 0) { // scale x for log-scale
            iscaled = HistogramScaling::log (255.0, (double)i);
        }

        double posX = (iscaled / 255.0) * (hsize - 1);
        double posY = vsize - 2 + val * (4 - vsize) / vsize;

        cr->line_to (posX, posY);
    }

    cr->line_to (hsize - 1, vsize - 1);
}

void HistogramArea::drawMarks(Cairo::RefPtr<Cairo::Context> &cr,
                              const LUTu &data, double scale, int hsize, int & ui, int & oi)
{
    int s = 8 * RTScalable::getScale();

    if(data[0] > scale) {
        cr->rectangle(0, (ui++)*s, s, s);
    }

    if(data[data.getUpperBound()] > scale) {
        cr->rectangle(hsize - s, (oi++)*s, s, s);
    }

    cr->fill();
}

void HistogramArea::drawParade(Cairo::RefPtr<Cairo::Context> &cr, int w, int h)
{
    // Arbitrary scale factor divided by current scale.
    const float scale = trace_brightness * 32.f * 255.f / waveform_scale;
    const int wave_width = rwave.width();
    const int wave_height = rwave.height();

    // See Cairo documentation on stride.
    const int cairo_stride = Cairo::ImageSurface::format_stride_for_width(Cairo::FORMAT_ARGB32, rwave.width());
    const auto buffer_size = static_cast<std::vector<unsigned char>::size_type>(wave_height) * cairo_stride;

    if (parade_buffer_r_dirty && needRed) {
        parade_buffer_r.assign(buffer_size, 0);
        assert(parade_buffer_r.size() % 4 == 0);

        for (int val = 0; val < wave_height; val++) {
            const int* const r_row = rwave[val];
            std::uint32_t* const buffer_r_row = reinterpret_cast<uint32_t*>(parade_buffer_r.data() + (255 - val) * cairo_stride);
            for (int col = 0; col < wave_width; col++) {
                const unsigned char r = std::min<float>(scale * r_row[col], 0xff);
                if (r != 0) {
                    buffer_r_row[col] = (r << 16) | (r << 24);
                }
            }
        }

        parade_buffer_r_dirty = false;
    }

    if (parade_buffer_g_dirty && needGreen) {
        parade_buffer_g.assign(buffer_size, 0);
        assert(parade_buffer_g.size() % 4 == 0);

        for (int val = 0; val < wave_height; val++) {
            const int* const g_row = gwave[val];
            std::uint32_t* const buffer_g_row = reinterpret_cast<uint32_t*>(parade_buffer_g.data() + (255 - val) * cairo_stride);
            for (int col = 0; col < wave_width; col++) {
                const unsigned char g = std::min<float>(scale * g_row[col], 0xff);
                if (g != 0) {
                    buffer_g_row[col] = (g << 8) | (g << 24);
                }
            }
        }

        parade_buffer_g_dirty = false;
    }

    if (parade_buffer_b_dirty && needBlue) {
        parade_buffer_b.assign(buffer_size, 0);
        assert(parade_buffer_b.size() % 4 == 0);

        for (int val = 0; val < wave_height; val++) {
            const int* const b_row = bwave[val];
            std::uint32_t* const buffer_b_row = reinterpret_cast<uint32_t*>(parade_buffer_b.data() + (255 - val) * cairo_stride);
            for (int col = 0; col < wave_width; col++) {
                const unsigned char b = std::min<float>(scale * b_row[col], 0xff);
                if (b != 0) {
                    const unsigned char green = b * rgb_B[1]; // Make blue easier to see.
                    buffer_b_row[col] = b | (green << 8) | (b << 24);
                }
            }
        }

        parade_buffer_b_dirty = false;
    }

    if (wave_buffer_luma_dirty && needLuma) {
        wave_buffer_luma.assign(buffer_size, 0);
        assert(wave_buffer_luma.size() % 4 == 0);

        for (int val = 0; val < wave_height; val++) {
            const int* const l_row = lwave[val];
            std::uint32_t* const buffer_row =
                reinterpret_cast<uint32_t*>(wave_buffer_luma.data() + (255 - val) * cairo_stride);
            for (int col = 0; col < wave_width; col++) {
                const unsigned char l = std::min<float>(scale * l_row[col], 0xff);
                buffer_row[col] = l | (l << 8) | (l << 16) | (l << 24);
            }
        }

        wave_buffer_luma_dirty = false;
    }

    std::vector<unsigned char*> buffers;
    if (needLuma) {
        buffers.push_back(wave_buffer_luma.data());
    }
    if (needRed) {
        buffers.push_back(parade_buffer_r.data());
    }
    if (needGreen) {
        buffers.push_back(parade_buffer_g.data());
    }
    if (needBlue) {
        buffers.push_back(parade_buffer_b.data());
    }

    auto orig_matrix = cr->get_matrix();
    const double display_wave_width = static_cast<double>(w) / buffers.size();
    for (unsigned i = 0; i < buffers.size(); i++) {
        Cairo::RefPtr<Cairo::ImageSurface> surface;
        cr->translate(i * display_wave_width, padding);
        cr->scale(display_wave_width / wave_width, (h - 2 * padding) / wave_height);
        surface = Cairo::ImageSurface::create(
            buffers[i], Cairo::FORMAT_ARGB32, wave_width, wave_height, cairo_stride);
        cr->set_source(surface, 0, 0);
        cr->set_operator(Cairo::OPERATOR_OVER);
        cr->paint();
        surface->finish();
        cr->set_matrix(orig_matrix);
    }
}

void HistogramArea::drawVectorscope(Cairo::RefPtr<Cairo::Context> &cr, int w, int h)
{
    if (scopeType != ScopeType::VECTORSCOPE_HC && scopeType != ScopeType::VECTORSCOPE_HS) {
        return;
    }

    const auto& vect = (scopeType == ScopeType::VECTORSCOPE_HC) ? vect_hc : vect_hs;
    auto& vect_buffer = (scopeType == ScopeType::VECTORSCOPE_HC) ? vect_hc_buffer : vect_hs_buffer;
    auto& vect_buffer_dirty = (scopeType == ScopeType::VECTORSCOPE_HC) ? vect_hc_buffer_dirty : vect_hs_buffer_dirty;

    const int vect_width = vect.width();
    const int vect_height = vect.height();
    // Arbitrary scale factor multiplied by vectorscope area and divided by
    // current scale.
    const float scale = trace_brightness * 8.f * vect_width * vect_height / vectorscope_scale;

    // See Cairo documentation on stride.
    const int cairo_stride = Cairo::ImageSurface::format_stride_for_width(Cairo::FORMAT_ARGB32, vect_width);

    if (vect_buffer_dirty && vectorscope_scale > 0) {
        if (vect_buffer.size() != static_cast<std::size_t>(cairo_stride) * vect_height) {
            vect_buffer.resize(static_cast<std::size_t>(cairo_stride) * vect_height);
        }

        assert(vect_buffer.size() % 4 == 0);

        for (int y = 0; y < vect_height; y++) {
            const int* const vect_row = vect[y];
            std::uint32_t* const buffer_row =
                reinterpret_cast<uint32_t*>(vect_buffer.data() + (vect_height - 1 - y) * cairo_stride);
            for (int x = 0; x < vect_width; x++) {
                const unsigned char value = std::min<float>(scale * vect_row[x], 0xff);
                buffer_row[x] = value | (value << 8) | (value << 16) | (value << 24);
            }
        }

        vect_buffer_dirty = false;
    }

    const bool fit_width =
        vect_width * (h - 2 * padding) > vect_height * (w - 2 * padding);
    const float scope_scale = fit_width ?
        (w - 2 * padding) / vect_width : (h - 2 * padding) / vect_height;
    const float scope_size = (vectorscope_scale > 0) ?
        scope_scale * std::max<double>(vect_width, vect_height) : std::min<float>(w, h) - 2 * padding;
    const float o_x = (w - scope_scale * vect_width) / 2;
    const float o_y = (h - scope_scale * vect_height) / 2;
    const double s = RTScalable::getScale();
    auto orig_matrix = cr->get_matrix();
    const double line_length = scope_size / 2.0;
    std::valarray<double> ch_ds(1);

    cr->translate(w / 2.0, h / 2.0);
    cr->set_line_width (1.0 * s);
    cr->set_antialias(Cairo::ANTIALIAS_SUBPIXEL);
    ch_ds[0] = 4;

    if (scopeType == ScopeType::VECTORSCOPE_HS) { // Hue-Saturation.
        // RYGCBM lines.
        cr->set_line_width (2.0 * s);
        constexpr double color_labels[6][3] = {
            {1, 0, 0}, // R
            {0, 1, 0}, // G
            {0, 0, 1}, // B
            {0, 1, 1}, // C
            {1, 0, 1}, // M
            {1, 1, 0}, // Y
        };
        for (int i = 0; i < 3; i++) {
            auto gradient = Cairo::LinearGradient::create(-line_length, 0, line_length, 0);
            const double (&color_1)[3] = color_labels[i];
            const double (&color_2)[3] = color_labels[i + 3];
            cr->set_source(gradient);
            gradient->add_color_stop_rgba(0, color_2[0], color_2[1], color_2[2], 0.5);
            gradient->add_color_stop_rgba(0.5, 1, 1, 1, 0.25);
            gradient->add_color_stop_rgba(1, color_1[0], color_1[1], color_1[2], 0.5);
            cr->move_to(-line_length, 0);
            cr->line_to(line_length, 0);
            cr->rotate_degrees(-120);
            cr->stroke();
        }
        cr->set_line_width (1.0 * s);
        cr->set_source_rgba (1, 1, 1, 0.25);
        // 100% saturation circle.
        cr->arc(0, 0, scope_size / 2.0, 0, 2 * RT_PI);
        cr->stroke();
        // 25%, 50%, and 75% saturation.
        cr->set_dash(ch_ds, 0);
        for (int i = 1; i < 4; i++) {
            cr->arc(0, 0, i * scope_size / 8.0, 0, 2 * RT_PI);
            cr->stroke();
        }
        // HSV skin tone line derived from -I axis of YIQ.
        cr->rotate(-0.134900 * RT_PI);
        cr->move_to(0, 0);
        cr->line_to(line_length, 0);
        cr->stroke();
        cr->unset_dash();
    } else if (scopeType == ScopeType::VECTORSCOPE_HC) { // Hue-Chroma.
        // a and b axes.
        Cairo::RefPtr<Cairo::LinearGradient> gradient;
        cr->set_line_width (2.0 * s);
        gradient = Cairo::LinearGradient::create(0, -line_length, 0, line_length);
        cr->set_source(gradient);
        gradient->add_color_stop_rgba(0, 1, 1, 0, 0.5); // "yellow"
        gradient->add_color_stop_rgba(0.5, 1, 1, 1, 0.25); // neutral
        gradient->add_color_stop_rgba(1, 0, 0, 1, 0.5); // "blue"
        cr->move_to(0, 0);
        cr->line_to(0, line_length);
        cr->move_to(0, 0);
        cr->line_to(0, -line_length);
        cr->stroke();
        gradient = Cairo::LinearGradient::create(-line_length, 0, line_length, 0);
        cr->set_source(gradient);
        gradient->add_color_stop_rgba(0, 0, 1, 0, 0.5); // "green"
        gradient->add_color_stop_rgba(0.5, 1, 1, 1, 0.25); // neutral
        gradient->add_color_stop_rgba(1, 1, 0, 1, 0.5); // "magenta"
        cr->move_to(0, 0);
        cr->line_to(line_length, 0);
        cr->move_to(0, 0);
        cr->line_to(-line_length, 0);
        cr->stroke();
        cr->set_source_rgba (1, 1, 1, 0.25);
        cr->set_line_width (1.0 * s);
        // 25%, 50%, 75%, and 100% of standard chroma range.
        cr->set_dash(ch_ds, 0);
        for (int i = 1; i <= 4; i++) {
            cr->arc(0, 0, i * scope_size / 8.0, 0, 2 * RT_PI);
            cr->stroke();
        }
        // CIELAB skin tone line, approximated by 50% saturation and
        // value along the HSV skin tone line.
        cr->rotate(-0.321713 * RT_PI);
        cr->move_to(0, 0);
        cr->line_to(line_length, 0);
        cr->stroke();
        cr->unset_dash();
    }
    cr->set_matrix(orig_matrix);

    // Vectorscope trace.
    if (vectorscope_scale > 0) {
        Cairo::RefPtr<Cairo::ImageSurface> surface = Cairo::ImageSurface::create(
            vect_buffer.data(), Cairo::FORMAT_ARGB32, vect_width, vect_height, cairo_stride);
        cr->translate(o_x, o_y);
        cr->scale(scope_scale, scope_scale);
        cr->set_source(surface, 0, 0);
        cr->set_operator(Cairo::OPERATOR_OVER);
        cr->paint();
        surface->finish();
        cr->set_matrix(orig_matrix);

        if (needPointer && pointer_red >= 0 && pointer_green >= 0 && pointer_blue >= 0) {
            float cx, cy;
            if (scopeType == ScopeType::VECTORSCOPE_HS) {
                float H, S, L;
                Color::rgb2hslfloat(pointer_red * 257.f, pointer_green * 257.f, pointer_blue * 257.f, H, S, L);
                cx = (w + scope_size * S * std::cos(H * 2 * RT_PI_F)) / 2;
                cy = (h - scope_size * S * std::sin(H * 2 * RT_PI_F)) / 2;
            } else {
                constexpr float ab_factor = 1.f / 256.f;
                cx = w / 2.f + scope_size * pointer_a * ab_factor;
                cy = h / 2.f - scope_size * pointer_b * ab_factor;
            }
            const float crosshair_size = 20.f * s;
            cr->set_source_rgba(1, 1, 1, 0.5);
            cr->move_to(cx - crosshair_size, cy);
            cr->line_to(cx + crosshair_size, cy);
            cr->move_to(cx, cy - crosshair_size);
            cr->line_to(cx, cy + crosshair_size);
            cr->stroke();
            cr->arc(cx, cy, 3 * s, 0, 2 * RT_PI);
            cr->set_source_rgb(1, 1, 1);
            cr->fill_preserve();
            cr->set_source_rgb(0, 0, 0);
            cr->set_line_width (1.0 * s);
            cr->stroke();
        }
    }
}

void HistogramArea::drawWaveform(Cairo::RefPtr<Cairo::Context> &cr, int w, int h)
{
    // Arbitrary scale factor divided by current scale.
    const float scale = trace_brightness * 32.f * 255.f / waveform_scale;
    const int wave_width = rwave.width();
    const int wave_height = rwave.height();

    // See Cairo documentation on stride.
    const int cairo_stride = Cairo::ImageSurface::format_stride_for_width(Cairo::FORMAT_ARGB32, rwave.width());
    const auto buffer_size = static_cast<std::vector<unsigned char>::size_type>(wave_height) * cairo_stride;

    if (wave_buffer_dirty && (needRed || needGreen || needBlue)) {
        wave_buffer.assign(buffer_size, 0);
        assert(wave_buffer.size() % 4 == 0);

        for (int val = 0; val < wave_height; val++) {
            const int* const r_row = rwave[val];
            const int* const g_row = gwave[val];
            const int* const b_row = bwave[val];
            std::uint32_t* const buffer_row = reinterpret_cast<uint32_t*>(wave_buffer.data() + (255 - val) * cairo_stride);
            for (int col = 0; col < wave_width; col++) {
                const unsigned char r = needRed ? std::min<float>(scale * r_row[col], 0xff) : 0;
                const unsigned char b = needBlue ? std::min<float>(scale * b_row[col], 0xff) : 0;
                const unsigned char g = needGreen ? rtengine::LIM(std::min<float>(scale * g_row[col], 0xff) + float(needBlue ? b * rgb_B[1] : 0), 0.f, float(0xff)) : 0;
                const unsigned char value = rtengine::max(r, g, b);
                if (value != 0) {
                    // Ensures correct order regardless of endianness.
                    buffer_row[col] = b | (g << 8) | (r << 16) | (value << 24);
                }
            }
        }

        wave_buffer_dirty = false;
    }

    if (wave_buffer_luma_dirty && needLuma) {
        wave_buffer_luma.assign(buffer_size, 0);
        assert(wave_buffer_luma.size() % 4 == 0);

        for (int val = 0; val < wave_height; val++) {
            const int* const l_row = lwave[val];
            std::uint32_t* const buffer_row =
                reinterpret_cast<uint32_t*>(wave_buffer_luma.data() + (255 - val) * cairo_stride);
            for (int col = 0; col < wave_width; col++) {
                const unsigned char l = std::min<float>(scale * l_row[col], 0xff);
                buffer_row[col] = l | (l << 8) | (l << 16) | (l << 24);
            }
        }

        wave_buffer_luma_dirty = false;
    }

    Cairo::RefPtr<Cairo::ImageSurface> surface;
    auto orig_matrix = cr->get_matrix();
    cr->translate(0, padding);
    cr->scale(static_cast<double>(w) / wave_width, (h - 2 * padding) / wave_height);
    if (needLuma) {
        surface = Cairo::ImageSurface::create(
            wave_buffer_luma.data(), Cairo::FORMAT_ARGB32, wave_width, wave_height, cairo_stride);
        cr->set_source(surface, 0, 0);
        cr->set_operator(Cairo::OPERATOR_OVER);
        cr->paint();
        surface->finish();
    }
    if (needRed || needGreen || needBlue) {
        surface = Cairo::ImageSurface::create(
            wave_buffer.data(), Cairo::FORMAT_ARGB32, wave_width, wave_height, cairo_stride);
        cr->set_source(surface, 0, 0);
        cr->set_operator(Cairo::OPERATOR_OVER);
        cr->paint();
        surface->finish();
    }
    cr->set_matrix(orig_matrix);
}

bool HistogramArea::on_draw(const ::Cairo::RefPtr< Cairo::Context> &cr)
{

    if (!updatePending() && (get_width() != oldwidth || get_height() != oldheight || isDirty())) {
        updateBackBuffer ();
    }

    const Glib::RefPtr<Gtk::StyleContext> style = get_style_context();
    style->render_background(cr, 0, 0, get_width(), get_height());
    copySurface(cr, NULL);
    style->render_frame (cr, 0, 0, get_width(), get_height());

    return true;
}

bool HistogramArea::on_button_press_event (GdkEventButton* event)
{
    if (!is_main_) {
        return true;
    }
    
    isPressed = true;
    movingPosition = event->x;

    if (event->type == GDK_2BUTTON_PRESS && event->button == 1
        && (scopeType == ScopeType::HISTOGRAM || scopeType == ScopeType::HISTOGRAM_RAW)
    ) {

        drawMode = (drawMode + 1) % 3;
        options.histogramDrawMode = (options.histogramDrawMode + 1) % 3;

        if (myDrawModeListener) {
            myDrawModeListener->toggleButtonMode ();
        }

        updateBackBuffer ();
        queue_draw ();
    }

    return true;
}

bool HistogramArea::on_button_release_event (GdkEventButton* event)
{
    if (!is_main_) {
        return true;
    }
    
    isPressed = false;
    return true;
}

bool HistogramArea::on_motion_notify_event (GdkEventMotion* event)
{
    if (!is_main_) {
        return true;
    }
    
    if (
        drawMode == 0
        && (scopeType == ScopeType::HISTOGRAM || scopeType == ScopeType::HISTOGRAM_RAW)
    ) {
        return false;
    }

    if (!isPressed) {
        return true;
    }

    if (scopeType == ScopeType::HISTOGRAM || scopeType == ScopeType::HISTOGRAM_RAW) { // Adjust log scale.
        double mod = 1 + (event->x - movingPosition) / get_width();

        factor /= mod;
        if (factor < 1.0)
            factor = 1.0;
        if (factor > 100.0)
            factor = 100.0;

        sigFactorChanged.emit(factor);

        setDirty(true);
        queue_draw ();
    } else if (
        scopeType == ScopeType::PARADE
        || scopeType == ScopeType::WAVEFORM
        || scopeType == ScopeType::VECTORSCOPE_HC
        || scopeType == ScopeType::VECTORSCOPE_HS
    ) { // Adjust brightness.
        constexpr float RANGE = MAX_BRIGHT / MIN_BRIGHT;
        double dx = (event->x - movingPosition) / get_width();
        float new_brightness = LIM<float>(trace_brightness * pow(RANGE, dx), MIN_BRIGHT, MAX_BRIGHT);
        setBrightness(new_brightness);
        movingPosition = event->x;
    }

    return true;
}

float HistogramArea::getBrightness(void)
{
    return trace_brightness;
}

void HistogramArea::setBrightness(float brightness)
{
    brightness = LIM<float>(brightness, MIN_BRIGHT, MAX_BRIGHT);
    if (brightness != trace_brightness) {
        parade_buffer_r_dirty = parade_buffer_g_dirty = parade_buffer_b_dirty = wave_buffer_dirty = wave_buffer_luma_dirty = vect_hc_buffer_dirty = vect_hs_buffer_dirty = true;
        trace_brightness = brightness;
        setDirty(true);
        queue_draw();

        signal_brightness_changed.emit(trace_brightness);
    }
}

HistogramArea::SignalBrightnessChanged HistogramArea::getBrighnessChangedSignal(void)
{
    return signal_brightness_changed;
}

HistogramArea::type_signal_factor_changed HistogramArea::signal_factor_changed()
{
    return sigFactorChanged;
}
