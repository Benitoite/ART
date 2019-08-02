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
#include "exposure.h"
#include "adjuster.h"
#include <sigc++/slot.h>
#include <iomanip>
#include "ppversion.h"
#include "edit.h"
#include "eventmapper.h"

using namespace rtengine;
using namespace rtengine::procparams;

Exposure::Exposure():
    FoldableToolPanel(this, "exposure", M("TP_EXPOSURE_LABEL"), false, true)
{
    auto m = ProcEventMapper::getInstance();
    EvClampOOG = m->newEvent(DARKFRAME, "HISTORY_MSG_CLAMPOOG");
    EvToolEnabled.set_action(DARKFRAME);

//----------- OOG clamping ----------------------------------
    clampOOG = Gtk::manage(new Gtk::CheckButton(M("TP_EXPOSURE_CLAMPOOG")));
    pack_start(*clampOOG);
    pack_start (*Gtk::manage (new  Gtk::HSeparator()));
    clampOOG->signal_toggled().connect(sigc::mem_fun(*this, &Exposure::clampOOGChanged));
    
//----------- Auto Levels ----------------------------------
    abox = Gtk::manage (new Gtk::HBox ());
    abox->set_spacing (4);

    autolevels = Gtk::manage (new Gtk::ToggleButton (M("TP_EXPOSURE_AUTOLEVELS")));
    autolevels->set_tooltip_markup (M("TP_EXPOSURE_AUTOLEVELS_TIP"));
    autoconn = autolevels->signal_toggled().connect( sigc::mem_fun(*this, &Exposure::autolevels_toggled) );

    lclip = Gtk::manage (new Gtk::Label (M("TP_EXPOSURE_CLIP")));
    lclip->set_tooltip_text (M("TP_EXPOSURE_CLIP_TIP"));

    sclip = Gtk::manage (new MySpinButton ());
    sclip->set_range (0.0, 0.99);
    sclip->set_increments (0.01, 0.10);
    sclip->set_value (0.02);
    sclip->set_digits (2);
    sclip->set_width_chars(4);
    sclip->set_max_width_chars(4);
    sclip->signal_value_changed().connect( sigc::mem_fun(*this, &Exposure::clip_changed) );

    neutral = Gtk::manage (new Gtk::Button (M("TP_NEUTRAL")));
    neutral->set_tooltip_text (M("TP_NEUTRAL_TIP"));
    neutralconn = neutral->signal_pressed().connect( sigc::mem_fun(*this, &Exposure::neutral_pressed) );
    neutral->show();

    abox->pack_start (*autolevels, true, true, 0);
    // pack_end is used for these controls as autolevels is replaceable using pack_start in batchmode
    abox->pack_end (*neutral, true, true, 0);
    abox->pack_end (*sclip, false, false, 0);
    abox->pack_end (*lclip, false, false, 0);
    pack_start (*abox);

//-------------- Highlight Reconstruction -----------------
    pack_start (*Gtk::manage (new  Gtk::HSeparator()));

    hrmode = Gtk::manage (new MyComboBoxText ());
    hrmode->append(M("TP_HLREC_OFF"));
    hrmode->append(M("TP_HLREC_BLEND"));
    hrmode->append(M("TP_HLREC_COLOR"));

    hrmode->set_active(ExposureParams::HR_OFF);
    hlrbox = Gtk::manage(new Gtk::HBox());
    Gtk::Label* lab = Gtk::manage(new Gtk::Label(M("TP_HLREC_LABEL") + ": "));
    hlrbox->pack_start(*lab, Gtk::PACK_SHRINK);//, 4);
    hlrbox->pack_start(*hrmode);
    pack_start (*hlrbox);

    hrmode->signal_changed().connect ( sigc::mem_fun(*this, &Exposure::hrmodeChanged) );

    //----------- Exposure Compensation ---------------------
    pack_start (*Gtk::manage (new  Gtk::HSeparator()));

    expcomp   = Gtk::manage (new Adjuster (M("TP_EXPOSURE_EXPCOMP"), -5, 12, 0.05, 0));
    expcomp->setLogScale(64, 0, true);
    pack_start (*expcomp);

    //----------- Highlight recovery & threshold -------------
    hlcompr = Gtk::manage (new Adjuster (M("TP_EXPOSURE_COMPRHIGHLIGHTS"), 0, 500, 1, 0));
    pack_start (*hlcompr);
    hlcomprthresh = Gtk::manage (new Adjuster (M("TP_EXPOSURE_COMPRHIGHLIGHTSTHRESHOLD"), 0, 100, 1, 0));
    pack_start (*hlcomprthresh);

//----------- Black Level & Compression -------------------
    black = Gtk::manage (new Adjuster (M("TP_EXPOSURE_BLACKLEVEL"), -16384, 32768, 50, 0));
    black->setLogScale(500, 0, true);
    pack_start (*black);
    shcompr = Gtk::manage (new Adjuster (M("TP_EXPOSURE_COMPRSHADOWS"), 0, 100, 1, 50));
    pack_start (*shcompr);

    pack_start (*Gtk::manage (new  Gtk::HSeparator()));

// --------- Set Up Listeners -------------
    expcomp->setAdjusterListener (this);
    black->setAdjusterListener (this);
    hlcompr->setAdjusterListener (this);
    hlcomprthresh->setAdjusterListener (this);
    shcompr->setAdjusterListener (this);
}


Exposure::~Exposure ()
{
    idle_register.destroy();
}


void Exposure::read(const ProcParams* pp)
{
    disableListener();

    autoconn.block (true);

    setEnabled(pp->exposure.enabled);
    
    autolevels->set_active (pp->exposure.autoexp);
    lastAuto = pp->exposure.autoexp;
    sclip->set_value (pp->exposure.clip);

    expcomp->setValue (pp->exposure.expcomp);
    black->setValue (pp->exposure.black);
    hlcompr->setValue (pp->exposure.hlcompr);
    hlcomprthresh->setValue (pp->exposure.hlcomprthresh);
    shcompr->setValue (pp->exposure.shcompr);

    if (!black->getAddMode()) {
        shcompr->set_sensitive(!((int)black->getValue () == 0));    //at black=0 shcompr value has no effect
    }
    
    if (!hlcompr->getAddMode()) {
        hlcomprthresh->set_sensitive(!((int)hlcompr->getValue () == 0));    //at hlcompr=0 hlcomprthresh value has no effect
    }

    clampOOG->set_active(pp->exposure.clampOOG);

    hrmode->set_active(pp->exposure.hrmode);
    autoconn.block (false);

    enableListener ();
}


void Exposure::write(ProcParams *pp)
{
    pp->exposure.enabled = getEnabled();
    pp->exposure.autoexp = autolevels->get_active();
    pp->exposure.clip = sclip->get_value ();
    pp->exposure.expcomp = expcomp->getValue ();
    pp->exposure.black = (int)black->getValue ();
    pp->exposure.hlcompr = (int)hlcompr->getValue ();
    pp->exposure.hlcomprthresh = (int)hlcomprthresh->getValue ();
    pp->exposure.shcompr = (int)shcompr->getValue ();
    pp->exposure.clampOOG = clampOOG->get_active();
    pp->exposure.hrmode = ExposureParams::HighlightReconstruction(hrmode->get_active_row_number());
}

void Exposure::hrmodeChanged ()
{
    if (listener && getEnabled()) {
        //setHistmatching(false);
        listener->panelChanged(EvHRMethod, hrmode->get_active_text());
    }
}


void Exposure::clampOOGChanged()
{
    if (listener && getEnabled()) {
        listener->panelChanged(EvClampOOG, clampOOG->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}



void Exposure::setRaw (bool raw)
{
    disableListener ();
    if (raw) {
        hrmode->set_sensitive(true);
    } else {
        hrmode->set_active(0);
        hrmode->set_sensitive(false);
    }
    enableListener ();
}


void Exposure::setDefaults (const ProcParams* defParams)
{
    expcomp->setDefault (defParams->exposure.expcomp);
    black->setDefault (defParams->exposure.black);
    hlcompr->setDefault (defParams->exposure.hlcompr);
    hlcomprthresh->setDefault (defParams->exposure.hlcomprthresh);
    shcompr->setDefault (defParams->exposure.shcompr);
}


void Exposure::adjusterChanged(Adjuster* a, double newval)
{
    // Switch off auto exposure if user changes sliders manually
    if (autolevels->get_active() && (a == expcomp || a == black || a == hlcompr || a == hlcomprthresh)) {
        autoconn.block(true);
        autolevels->set_active (false);
        autoconn.block(false);
        autolevels->set_inconsistent (false);
    }

    if (!listener || !getEnabled()) {
        return;
    }

    Glib::ustring costr;

    if (a == expcomp) {
        costr = Glib::ustring::format (std::setw(3), std::fixed, std::setprecision(2), a->getValue());
    } else {
        costr = Glib::ustring::format ((int)a->getValue());
    }

    if (a == expcomp) {
        listener->panelChanged (EvExpComp, costr);
    } else if (a == black) {
        listener->panelChanged (EvBlack, costr);

        if (!black->getAddMode()) {
            shcompr->set_sensitive(!((int)black->getValue () == 0));    //at black=0 shcompr value has no effect
        }
    } else if (a == hlcompr) {
        listener->panelChanged (EvHLCompr, costr);
        
        if (!hlcompr->getAddMode()) {
            hlcomprthresh->set_sensitive(!((int)hlcompr->getValue () == 0));    //at hlcompr=0 hlcomprthresh value has no effect
        }
    } else if (a == hlcomprthresh) {
        listener->panelChanged (EvHLComprThreshold, costr);
    } else if (a == shcompr) {
        listener->panelChanged (EvSHCompr, costr);
    }
}

void Exposure::adjusterAutoToggled(Adjuster* a, bool newval)
{
}

void Exposure::neutral_pressed ()
{
// This method deselects auto levels and HL reconstruction auto
// and sets neutral values to params in exposure panel
    
    autolevels->set_active (false);
    autolevels->set_inconsistent (false);

    expcomp->setValue(0);
    hlcompr->setValue(0);
    hlcomprthresh->setValue(0);
    black->setValue(0);
    shcompr->setValue(50);
    hrmode->set_active(0);

    if (!black->getAddMode()) {
        shcompr->set_sensitive(!((int)black->getValue () == 0));    //at black=0 shcompr value has no effect
    }
    
    if (!hlcompr->getAddMode()) {
        hlcomprthresh->set_sensitive(!((int)hlcompr->getValue () == 0));    //at hlcompr=0 hlcomprthresh value has no effect
    }

    listener->panelChanged (EvNeutralExp, M("GENERAL_ENABLED"));
}


void Exposure::autolevels_toggled ()
{
    if (listener && getEnabled()) {
        if (autolevels->get_active()) {
            listener->panelChanged (EvAutoExp, M("GENERAL_ENABLED"));
            waitForAutoExp ();

            if (!black->getAddMode()) {
                shcompr->set_sensitive(!((int)black->getValue () == 0));    //at black=0 shcompr value has no effect
            }
            
            if (!hlcompr->getAddMode()) {
                hlcomprthresh->set_sensitive(!((int)hlcompr->getValue () == 0));    //at hlcompr=0 hlcomprthresh value has no effect
            }

        } else {
            listener->panelChanged (EvFixedExp, M("GENERAL_DISABLED"));
        }
    }
}


void Exposure::clip_changed ()
{
    clipDirty = true;

    if (autolevels->get_active() && listener) {
        Glib::signal_idle().connect (sigc::mem_fun(*this, &Exposure::clip_changed_));
    }
}

bool Exposure::clip_changed_ ()
{

    if (listener && getEnabled()) {
        listener->panelChanged (EvClip, Glib::ustring::format (std::setprecision(5), sclip->get_value()));

        waitForAutoExp ();
    }

    return false;
}

void Exposure::waitForAutoExp ()
{
    sclip->set_sensitive (false);
    expcomp->setEnabled (false);
    black->setEnabled (false);
    hlcompr->setEnabled (false);
    hlcomprthresh->setEnabled (false);
    shcompr->setEnabled (false);
    hrmode->set_sensitive(false);
}


void Exposure::enableAll ()
{
    sclip->set_sensitive (true);
    expcomp->setEnabled (true);
    black->setEnabled (true);
    hlcompr->setEnabled (true);
    hlcomprthresh->setEnabled (true);
    shcompr->setEnabled (true);
    hrmode->set_sensitive(true);
}


void Exposure::trimValues (rtengine::procparams::ProcParams* pp)
{
    expcomp->trimValue(pp->exposure.expcomp);
    hlcompr->trimValue(pp->exposure.hlcompr);
    hlcomprthresh->trimValue(pp->exposure.hlcomprthresh);
    black->trimValue(pp->exposure.black);
    shcompr->trimValue(pp->exposure.shcompr);
}


void Exposure::autoExpChanged(double expcomp, int bright, int contr, int black, int hlcompr, int hlcomprthresh, bool hlrecons)
{
    nextBlack = black;
    nextExpcomp = expcomp;
    nextHlcompr = hlcompr;
    nextHlcomprthresh = hlcomprthresh;
    nextHLRecons = hlrecons;

    idle_register.add(
        [this]() -> bool
        {
            GThreadLock lock; // All GUI access from idle_add callbacks or separate thread HAVE to be protected
            disableListener ();
            enableAll ();
            this->expcomp->setValue (nextExpcomp);
            this->black->setValue (nextBlack);
            this->hlcompr->setValue (nextHlcompr);
            this->hlcomprthresh->setValue (nextHlcomprthresh);
            this->hrmode->set_active(nextHLRecons ? ExposureParams::HR_BLEND : ExposureParams::HR_OFF);
            
            if (nextHLRecons) {
                hlrbox->show();
            } else {
                hlrbox->hide();
            }

            if (!this->black->getAddMode()) {
                this->shcompr->set_sensitive(!((int)this->black->getValue () == 0));    //at black=0 shcompr value has no effect
            }
            
            if (!this->hlcompr->getAddMode()) {
                this->hlcomprthresh->set_sensitive(!((int)this->hlcompr->getValue () == 0));    //at hlcompr=0 hlcomprthresh value has no effect
            }
            
            enableListener ();
            
            return false;
        });
}
