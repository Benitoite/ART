/* -*- C++ -*-
 *  
 *  This file is part of RawTherapee.
 */
#ifndef _GRADIENT_H_
#define _GRADIENT_H_

#include <gtkmm.h>
#include "adjuster.h"
#include "toolpanel.h"
#include "edit.h"
#include "guiutils.h"

class Gradient : public ToolParamBlock, public AdjusterListener, public FoldableToolPanel, public EditSubscriber
{

private:
    int lastObject;

protected:
    Gtk::HBox *editHBox;
    Gtk::ToggleButton* edit;
    Adjuster* degree;
    Adjuster* feather;
    Adjuster* strength;
    Adjuster* centerX;
    Adjuster* centerY;
    double draggedPointOldAngle;
    double draggedPointAdjusterAngle;
    double draggedFeatherOffset;
    rtengine::Coord draggedCenter;
    sigc::connection editConn;

    void editToggled ();

public:

    Gradient();
    ~Gradient() override;

    void read(const rtengine::procparams::ProcParams* pp) override;
    void write(rtengine::procparams::ProcParams* pp) override;
    void setDefaults(const rtengine::procparams::ProcParams* defParams) override;
    void adjusterChanged (Adjuster* a, double newval) override;
    void adjusterAutoToggled(Adjuster* a, bool newval) override;
    void enabledChanged  () override;
    void trimValues          (rtengine::procparams::ProcParams* pp) override;
    void updateGeometry  (const int centerX, const int centerY, const double feather, const double degree, const int fullWidth=-1, const int fullHeight=-1);

    void setEditProvider (EditDataProvider* provider) override;

    // EditSubscriber interface
    CursorShape getCursor(const int objectID) override;
    bool mouseOver(const int modifierKey) override;
    bool button1Pressed(const int modifierKey) override;
    bool button1Released() override;
    bool drag1(const int modifierKey) override;
    void switchOffEditMode () override;
};

#endif
