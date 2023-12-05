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
#ifndef __TOOLPANEL__
#define __TOOLPANEL__

#include <gtkmm.h>
#include <glibmm.h>
#include "../rtengine/rtengine.h"
#include "../rtengine/procparams.h"
#include "../rtengine/tweakoperator.h"
#include "guiutils.h"
#include "multilangmgr.h"
#include "paramsedited.h"
#include "edit.h"
#include "pparamschangelistener.h"
#include "shortcutmanager.h"

class ToolPanel;
class FoldableToolPanel;

class ToolPanelListener
{
public:
    virtual ~ToolPanelListener() = default;
    /// @brief Ask to refresh the preview not triggered by a parameter change (e.g. 'On Preview' editing).
    virtual void refreshPreview(const rtengine::ProcEvent& event) = 0;
    /// @brief Used to notify all listeners that a parameters has been effectively changed
    virtual void panelChanged(const rtengine::ProcEvent& event, const Glib::ustring& descr) = 0;
    /// @brief Set the TweakOperator to the StagedImageProcessor, to let some tool enter into special modes
    virtual void setTweakOperator (rtengine::TweakOperator *tOperator) = 0;
    /// @brief Unset the TweakOperator to the StagedImageProcessor
    virtual void unsetTweakOperator (rtengine::TweakOperator *tOperator) = 0;
};

/// @brief This class control the space around the group of tools inside a tab, as well as the space separating each tool. */
class ToolVBox: public Gtk::Box {
public:
    ToolVBox();
};

/// @brief This class control the space around a tool's block of parameter. */
class ToolParamBlock: public Gtk::Box {
public:
    ToolParamBlock();
    Gtk::SizeRequestMode get_request_mode_vfunc () const override;
};


class ToolPanel {
protected:
    Glib::ustring toolName;
    ToolPanelListener* listener;
    ToolPanelListener* tmp;
    bool need100Percent;

public:
    ToolPanel(Glib::ustring toolName = "", bool need11 = false) : toolName(toolName), listener(nullptr), tmp(nullptr), need100Percent(need11) {}
    virtual ~ToolPanel() {}

    virtual void setParent(Gtk::Box* parent) {}
    virtual Gtk::Box *getParent() { return nullptr; }
    virtual MyExpander *getExpander() { return nullptr; }
    virtual void setExpanded(bool expanded) {}
    virtual bool getExpanded() { return false; }
    virtual void setListener(ToolPanelListener *tpl) { listener = tpl; }
    virtual void setEditProvider(EditDataProvider *provider) {}
    virtual void read(const rtengine::procparams::ProcParams *pp) {}
    virtual void write(rtengine::procparams::ProcParams *pp) {}
//    virtual void read(const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr) {}
//    virtual void write(rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr) {}
    virtual void trimValues(rtengine::procparams::ProcParams *pp) {}
    virtual void setDefaults(const rtengine::procparams::ProcParams *defParams) {}
    virtual void autoOpenCurve() {}

    /** @brief Disable the event broadcasting mechanism
     *
     * @return Return the previous state of the broadcast (true: enabled ; false: disabled)
     */
    bool disableListener()
    {
        if (tmp == nullptr) {
            tmp = listener;
        }

        bool prevState = listener != nullptr;
        listener = nullptr;
        return prevState;
    }

    /** @brief Enable the event broadcasting mechanism
     */
    void enableListener()
    {
        if (tmp != nullptr) {
            listener = tmp;
        }

        tmp = nullptr;
    }

    virtual Glib::ustring getToolName() { return toolName; }
    virtual PParamsChangeListener *getPParamsChangeListener() { return nullptr; }

    virtual void registerShortcuts(ToolShortcutManager *mgr) {}
};


class FoldableToolPanel: public ToolPanel {

protected:
    Gtk::Box* parentContainer;
    MyExpander* exp;
    bool lastEnabled;
    sigc::connection enaConn;
    void foldThemAll (GdkEventButton* event);
    void enabled_toggled();
    rtengine::ProcEvent EvToolEnabled;
    rtengine::ProcEvent EvToolReset;

    Gtk::EventBox *imageEvBox;

    bool on_enter_leave_reset(GdkEventCrossing *event);
    bool on_reset_change(GdkEventButton *event);

public:

    FoldableToolPanel(Gtk::Box *content, const Glib::ustring &toolName, const Glib::ustring &UILabel, bool need11=false, bool useEnabled=false, bool useReset=false);

    MyExpander* getExpander() override
    {
        return exp;
    }
    void setExpanded (bool expanded) override
    {
        if (exp) {
            exp->set_expanded( expanded );
        }
    }

    void hide() {
        if (exp) {  // conditional hide
            exp->hide();
        }
    }

    void show() {
        if (exp) {                // always show
            exp->show();
        }
    }
    bool getExpanded () override
    {
        if (exp) {
            return exp->get_expanded();
        }

        return false;
    }
    void setParent (Gtk::Box* parent) override
    {
        parentContainer = parent;
    }
    Gtk::Box* getParent () override
    {
        return parentContainer;
    }

    virtual void enabledChanged();

    bool getUseEnabled ()
    {
        if (exp) {
            return exp->getUseEnabled();
        } else {
            return true;
        }
    }
    bool getEnabled();  // related to the enabled/disabled state
    void setEnabled(bool isActive);  // related to the enabled/disabled state
    void setEnabledTooltipMarkup(Glib::ustring tooltipMarkup);
    void setEnabledTooltipText(Glib::ustring tooltipText);
    bool get_inconsistent();  // related to the enabled/disabled state
    void set_inconsistent(bool isInconsistent);  // related to the enabled/disabled state
    void setGrayedOut(bool doGrayOut); // Set whether the tool should be disabled, collapsed and grayed-out.

    void setLevel (int level);

    // Functions that want to receive an enabled/disabled event from this class
    // will have to receive it from MyExpander directly, we do not create
    // a relaying event
    MyExpander::type_signal_enabled_toggled signal_enabled_toggled()
    {
        return exp->signal_enabled_toggled();
    }

    virtual void toolReset(bool to_initial)
    {
    }

    Glib::ustring getUILabel() const { return EvToolReset.get_message(); }
};

#endif
