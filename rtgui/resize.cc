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
#include <iomanip>
#include "resize.h"
#include "guiutils.h"
#include "eventmapper.h"

using namespace rtengine;
using namespace rtengine::procparams;

Resize::Resize () : FoldableToolPanel(this, "resize", M("TP_RESIZE_LABEL"), false, true), maxw(100000), maxh(100000)
{
    auto m = ProcEventMapper::getInstance();
    EvResizeAllowUpscaling = m->newEvent(RESIZE, "HISTORY_MSG_RESIZE_ALLOWUPSCALING");

    cropw = 0;
    croph = 0;

    Gtk::Table* combos = Gtk::manage (new Gtk::Table (2, 2));

    appliesTo = Gtk::manage (new MyComboBoxText ());
    appliesTo->append (M("TP_RESIZE_CROPPEDAREA"));
    appliesTo->append (M("TP_RESIZE_FULLIMAGE"));
    appliesTo->set_active (0);

    Gtk::Label *label = Gtk::manage (new Gtk::Label (M("TP_RESIZE_APPLIESTO")));
    label->set_alignment(0., 0.);
    combos->attach (*label, 0, 1, 0, 1, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    combos->attach (*appliesTo, 1, 2, 0, 1, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    // See Resize::methodChanged() when adding a new method.
    method = Gtk::manage (new MyComboBoxText ());
    method->append (M("TP_RESIZE_LANCZOS"));
    method->append (M("TP_RESIZE_NEAREST"));
    method->set_active (0);

    label = Gtk::manage (new Gtk::Label (M("TP_RESIZE_METHOD")));
    label->set_alignment(0., 0.);
    combos->attach (*label, 0, 1, 1, 2, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    combos->attach (*method, 1, 2, 1, 2, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    spec = Gtk::manage (new MyComboBoxText ());
    spec->append (M("TP_RESIZE_SCALE"));
    spec->append (M("TP_RESIZE_WIDTH"));
    spec->append (M("TP_RESIZE_HEIGHT"));
    spec->append (M("TP_RESIZE_FITBOX"));
    spec->set_active (0);

    label = Gtk::manage (new Gtk::Label (M("TP_RESIZE_SPECIFY")));
    label->set_alignment(0., 0.);
    combos->attach (*label, 0, 1, 2, 3, Gtk::SHRINK, Gtk::SHRINK, 2, 2);
    combos->attach (*spec, 1, 2, 2, 3, Gtk::EXPAND | Gtk::FILL, Gtk::SHRINK, 2, 2);

    pack_start (*combos, Gtk::PACK_SHRINK, 4);

    scale = new Adjuster (M("TP_RESIZE_SCALE"), 0.01, MAX_SCALE, 0.01, 1.);
    scale->setAdjusterListener (this);
    scale->setLogScale(MAX_SCALE * 10.f, 1.f, true);

    pack_start (*scale, Gtk::PACK_SHRINK, 4);

    sizeBox = Gtk::manage (new Gtk::VBox ());

    Gtk::HBox* sbox = Gtk::manage (new Gtk::HBox ());
    Gtk::HBox* wbox = Gtk::manage (new Gtk::HBox ());
    Gtk::HBox* hbox = Gtk::manage (new Gtk::HBox ());
    w = Gtk::manage (new MySpinButton ());
    h = Gtk::manage (new MySpinButton ());
    wbox->set_spacing(3);
    wbox->pack_start (*Gtk::manage (new Gtk::Label (M("TP_RESIZE_W"))), Gtk::PACK_SHRINK, 0);
    wbox->pack_start (*w);
    hbox->set_spacing(3);
    hbox->pack_start (*Gtk::manage (new Gtk::Label (M("TP_RESIZE_H"))), Gtk::PACK_SHRINK, 0);
    hbox->pack_start (*h);
    sbox->set_spacing(4);
    sbox->pack_start (*wbox);
    sbox->pack_start (*hbox);

    sizeBox->pack_start (*sbox, Gtk::PACK_SHRINK, 0);
    sizeBox->show_all ();
    sizeBox->reference ();

    allowUpscaling = Gtk::manage(new Gtk::CheckButton(M("TP_RESIZE_ALLOW_UPSCALING")));
    pack_start(*allowUpscaling);
    allowUpscaling->signal_toggled().connect(sigc::mem_fun(*this, &Resize::allowUpscalingChanged));

    w->set_digits (0);
    w->set_increments (1, 100);
    w->set_value (800);
    w->set_range (32, MAX_SCALE * maxw);

    h->set_digits (0);
    h->set_increments (1, 100);
    h->set_value (600);
    h->set_range (32, MAX_SCALE * maxh);

    wconn = w->signal_value_changed().connect ( sigc::mem_fun(*this, &Resize::entryWChanged), true);
    hconn = h->signal_value_changed().connect ( sigc::mem_fun(*this, &Resize::entryHChanged), true);
    aconn = appliesTo->signal_changed().connect ( sigc::mem_fun(*this, &Resize::appliesToChanged) );
    method->signal_changed().connect ( sigc::mem_fun(*this, &Resize::methodChanged) );
    sconn = spec->signal_changed().connect ( sigc::mem_fun(*this, &Resize::specChanged) );

    packBox = Gtk::manage (new ToolParamBlock ());
    pack_end (*packBox);
    packBox->hide();
    packBox->set_tooltip_markup (M("TP_PRSHARPENING_TOOLTIP"));

    show_all();
}

Resize::~Resize ()
{
    idle_register.destroy();
    delete scale;
    delete sizeBox;
}

void Resize::read(const ProcParams* pp)
{

    disableListener ();
    aconn.block (true);
    wconn.block (true);
    hconn.block (true);
    sconn.block (true);
    scale->block(true);

    scale->setValue (pp->resize.scale);
    w->set_value (pp->resize.width);
    h->set_value (pp->resize.height);
    setEnabled (pp->resize.enabled);
    spec->set_active (pp->resize.dataspec);
    allowUpscaling->set_active(pp->resize.allowUpscaling);
    updateGUI();

    appliesTo->set_active (0);

    if (pp->resize.appliesTo == "Cropped area") {
        appliesTo->set_active (0);
    } else if (pp->resize.appliesTo == "Full image") {
        appliesTo->set_active (1);
    }

    if (pp->resize.method == "Lanczos") {
        method->set_active (0);
    } else if (pp->resize.method == "Nearest") {
        method->set_active (1);
    } else {
        method->set_active (0);
    }

    wDirty = false;
    hDirty = false;

    scale->block(false);
    sconn.block (false);
    wconn.block (false);
    hconn.block (false);
    aconn.block (false);
    enableListener ();
}

void Resize::write(ProcParams* pp)
{
    int dataSpec = spec->get_active_row_number();

    pp->resize.scale  = scale->getValue();

    pp->resize.appliesTo = "Cropped area";

    if (appliesTo->get_active_row_number() == 0) {
        pp->resize.appliesTo = "Cropped area";
    } else if (appliesTo->get_active_row_number() == 1) {
        pp->resize.appliesTo = "Full image";
    }

    pp->resize.method = "Lanczos";

    if (method->get_active_row_number() == 0) {
        pp->resize.method = "Lanczos";
    } else if (method->get_active_row_number() == 1) {
        pp->resize.method = "Nearest";
    }

    pp->resize.dataspec = dataSpec;
    pp->resize.width = w->get_value_as_int ();
    pp->resize.height = h->get_value_as_int ();
    pp->resize.enabled = getEnabled ();
    //printf("  L:%d   H:%d\n", pp->resize.width, pp->resize.height);

    pp->resize.allowUpscaling = allowUpscaling->get_active();
}

void Resize::setDefaults(const ProcParams* defParams)
{
    scale->setDefault (defParams->resize.scale);
}

void Resize::adjusterChanged(Adjuster* a, double newval)
{
    wconn.block (true);
    hconn.block (true);
    h->set_value ((croph && appliesTo->get_active_row_number() == 0 ? croph : maxh) * a->getValue ());
    w->set_value ((cropw && appliesTo->get_active_row_number() == 0 ? cropw : maxw) * a->getValue ());
    wconn.block (false);
    hconn.block (false);

    if (listener && getEnabled()) {
        listener->panelChanged (EvResizeScale, Glib::ustring::format (std::setw(5), std::fixed, std::setprecision(2), scale->getValue()));
    }
}

void Resize::adjusterAutoToggled(Adjuster* a, bool newval)
{
}

int Resize::getComputedWidth()
{

    if (cropw && appliesTo->get_active_row_number() == 0)
        // we use the crop dimensions
    {
        return (int)((double)(cropw) * (h->get_value() / (double)(croph)) + 0.5);
    } else
        // we use the image dimensions
    {
        return (int)((double)(maxw) * (h->get_value() / (double)(maxh)) + 0.5);
    }
}

int Resize::getComputedHeight()
{

    if (croph && appliesTo->get_active_row_number() == 0)
        // we use the crop dimensions
    {
        return (int)((double)(croph) * (w->get_value() / (double)(cropw)) + 0.5);
    } else
        // we use the image dimensions
    {
        return (int)((double)(maxh) * (w->get_value() / (double)(maxw)) + 0.5);
    }
}

void Resize::appliesToChanged ()
{

    //printf("\nPASSAGE EN MODE \"%s\"\n\n", appliesTo->get_active_text().c_str());
    setDimensions();

    if (listener && getEnabled()) {
        //printf("Appel du listener\n");
        listener->panelChanged (EvResizeAppliesTo, appliesTo->get_active_text());
    }
}

void Resize::methodChanged ()
{

    if (listener && getEnabled()) {
        listener->panelChanged (EvResizeMethod, method->get_active_text());
    }

    // Post-resize Sharpening assumes the image is in Lab space, and currently Lanczos is the only method which uses that space, and Lanczos is on row 0.
    if (method->get_active_row_number() == 0) {
        packBox->set_sensitive(true);
    } else {
        packBox->set_sensitive(false);
    }
}

void Resize::update (bool isCropped, int cw, int ch, int ow, int oh)
{

    // updating crop values now
    if (isCropped) {
        cropw = cw;
        croph = ch;
    } else {
        cropw = 0;
        croph = 0;
    }

    // updating the full image dimensions
    if (ow && oh) {
        maxw = ow;
        maxh = oh;
    }

    // updating the GUI synchronously
    setDimensions();
}

void Resize::sizeChanged(int mw, int mh, int ow, int oh)
{
    // updating max values now
    maxw = ow;
    maxh = oh;

    // updating the GUI synchronously
    setDimensions();
}

void Resize::setDimensions ()
{
    idle_register.add(
        [this]() -> bool
        {
            wconn.block(true);
            hconn.block(true);
            scale->block(true);

            int refw, refh;

            if (appliesTo->get_active_row_number() == 0 && cropw) {
                // Applies to Cropped area
                refw = cropw;
                refh = croph;
            } else {
                // Applies to Full image or crop is disabled
                refw = maxw;
                refh = maxh;
            }

            w->set_range(32, MAX_SCALE * refw);
            h->set_range(32, MAX_SCALE * refh);

            switch (spec->get_active_row_number()) {
                case 0: {
                    // Scale mode
                    w->set_value(static_cast<double>(static_cast<int>(static_cast<double>(refw) * scale->getValue() + 0.5)));
                    h->set_value(static_cast<double>(static_cast<int>(static_cast<double>(refh) * scale->getValue() + 0.5)));
                    break;
                }

                case 1: {
                    // Width mode
                    const double tmp_scale = w->get_value() / static_cast<double>(refw);
                    scale->setValue(tmp_scale);
                    h->set_value(static_cast<double>(static_cast<int>(static_cast<double>(refh) * tmp_scale + 0.5)));
                    break;
                }

                case 2: {
                    // Height mode
                    const double tmp_scale = h->get_value() / static_cast<double>(refh);
                    scale->setValue(tmp_scale);
                    w->set_value(static_cast<double>(static_cast<int>(static_cast<double>(refw) * tmp_scale + 0.5)));
                    break;
                }

                case 3: {
                    // Bounding box mode
                    const double tmp_scale =
                        w->get_value() / h->get_value() < static_cast<double>(refw) / static_cast<double>(refh)
                            ? w->get_value() / static_cast<double>(refw)
                            : h->get_value() / static_cast<double>(refh);

                    scale->setValue(tmp_scale);
                    break;
                }

                default: {
                    break;
                }
            }

            scale->block(false);
            wconn.block(false);
            hconn.block(false);

            return false;
        }
    );
}

void Resize::fitBoxScale()
{
    double tmpScale;
    double neww = w->get_value ();
    double newh = h->get_value ();

    if (cropw && appliesTo->get_active_row_number() == 0) {
        // we use the crop dimensions
        if (((double)(cropw) / (double)(croph)) > (neww / newh)) {
            // the new scale is given by the image width
            tmpScale = neww / (double)(cropw);
        } else {
            // the new scale is given by the image height
            tmpScale = newh / (double)(croph);
        }
    } else {
        // we use the image dimensions
        if (((double)(maxw) / (double)(maxh)) > (neww / newh)) {
            // the new scale is given by the image width
            tmpScale = neww / (double)(maxw);
        } else {
            // the new scale is given by the image height
            tmpScale = newh / (double)(maxh);
        }
    }

    scale->setValue (tmpScale);
}

void Resize::entryWChanged ()
{

    wDirty = true;

    // updating width
    if (spec->get_active_row_number() == 3) {
        // Fit box mode
        fitBoxScale();
    } else {
        // Other modes
        hconn.block (true);
        scale->block (true);

        h->set_value ((double)(getComputedHeight()));
        scale->setValue (w->get_value () / (cropw && appliesTo->get_active_row_number() == 0 ? (double)cropw : (double)maxw));

        scale->block (false);
        hconn.block (false);
    }

    if (listener) {
        if (spec->get_active_row_number() == 3) {
            notifyBBox();
        } else {
            if (getEnabled()) {
                listener->panelChanged (EvResizeWidth, Glib::ustring::format (w->get_value_as_int()));
            }
        }
    }
}

void Resize::entryHChanged ()
{

    hDirty = true;

    if (listener) {
        if (spec->get_active_row_number() == 3) {
            // Fit box mode
            fitBoxScale();
        } else {
            // Other modes
            wconn.block (true);
            scale->block (true);

            w->set_value ((double)(getComputedWidth()));
            scale->setValue (h->get_value () / (croph && appliesTo->get_active_row_number() == 0 ? (double)croph : (double)maxh));

            scale->block (false);
            wconn.block (false);
        }
    }

    if (listener) {
        if (spec->get_active_row_number() == 3) {
            notifyBBox();
        } else {
            if (getEnabled()) {
                listener->panelChanged (EvResizeHeight, Glib::ustring::format (h->get_value_as_int()));
            }
        }
    }
}

void Resize::specChanged ()
{

    switch (spec->get_active_row_number()) {
    case (0):
        // Scale mode
        scale->sliderChanged ();
        break;

    case (1):
        // Width mode
        w->set_value((double)(getComputedWidth()));
        entryWChanged ();
        break;

    case (2):
        // Height mode
        h->set_value((double)(getComputedHeight()));
        entryHChanged ();
        break;

    case (3):
        // Bounding box mode
        notifyBBox();
        break;

    default:
        break;
    }

    updateGUI();
}

void Resize::updateGUI ()
{

    removeIfThere (this, scale, false);
    removeIfThere (this, sizeBox, false);

    switch (spec->get_active_row_number()) {
    case (0):
        // Scale mode
        pack_start (*scale, Gtk::PACK_SHRINK, 4);
        reorder_child(*allowUpscaling, 4);
        break;

    case (1):
        // Width mode
        pack_start (*sizeBox, Gtk::PACK_SHRINK, 4);
        reorder_child(*allowUpscaling, 4);
        w->set_sensitive (true);
        h->set_sensitive (false);
        break;

    case (2):
        // Height mode
        pack_start (*sizeBox, Gtk::PACK_SHRINK, 4);
        reorder_child(*allowUpscaling, 4);
        w->set_sensitive (false);
        h->set_sensitive (true);
        break;

    case (3):
        // Bounding box mode
        pack_start (*sizeBox, Gtk::PACK_SHRINK, 4);
        reorder_child(*allowUpscaling, 4);
        w->set_sensitive (true);
        h->set_sensitive (true);
        break;

    default:
        break;
    }
}

void Resize::notifyBBox()
{
    if (listener && getEnabled()) {
        listener->panelChanged (EvResizeBoundingBox, Glib::ustring::compose("(%1x%2)", (int)w->get_value(), (int)h->get_value() ));
    }
}


void Resize::enabledChanged ()
{

    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvResizeEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvResizeEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvResizeEnabled, M("GENERAL_DISABLED"));
        }
    }
}


void Resize::allowUpscalingChanged()
{

    if (listener) {
        if (allowUpscaling->get_inconsistent()) {
            listener->panelChanged(EvResizeAllowUpscaling, M("GENERAL_UNCHANGED"));
        } else if (allowUpscaling->get_active()) {
            listener->panelChanged(EvResizeAllowUpscaling, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged(EvResizeAllowUpscaling, M("GENERAL_DISABLED"));
        }
    }
}


void Resize::trimValues (rtengine::procparams::ProcParams* pp)
{

    scale->trimValue(pp->resize.scale);
}
