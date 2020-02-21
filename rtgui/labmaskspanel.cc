/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright 2018 Alberto Griggio <alberto.griggio@gmail.com>
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

#include "labmaskspanel.h"
#include "eventmapper.h"
#include "../rtengine/iccstore.h"

#include <iostream>

namespace {

using Shape = rtengine::procparams::AreaMask::Shape;

constexpr int ID_HUE_MASK = 5;

inline bool hasMask(const rtengine::procparams::Mask &m, const std::vector<double> &dflt, const std::vector<double> &mask)
{
    return m.parametricMask.enabled && !(mask.empty() || mask[0] == FCT_Linear || mask == dflt);
}


class DeltaEArea: public Gtk::DrawingArea, public BackBuffer, public EditSubscriber {
public:
    DeltaEArea():
        Gtk::DrawingArea(),
        EditSubscriber(ET_OBJECTS),
        L_(0.f), C_(0.f), H_(0.f)
    {
        // Editing geometry; create the spot rectangle
        Rectangle *spotRect = new Rectangle();
        spotRect->filled = false;
    
        visibleGeometry.push_back(spotRect);

        // Stick a dummy rectangle over the whole image in mouseOverGeometry.
        // This is to make sure the getCursor() call is fired everywhere.
        Rectangle *imgRect = new Rectangle();
        imgRect->filled = true;
        mouseOverGeometry.push_back(imgRect);
    }

    ~DeltaEArea()
    {
        for (auto geometry : visibleGeometry) {
            delete geometry;
        }

        for (auto geometry : mouseOverGeometry) {
            delete geometry;
        }        
    }

    void setColor(float L, float C, float H)
    {
        L_ = L;
        C_ = C;
        H_ = H;
        setDirty(true);
        queue_draw();
    }

    bool on_draw(const ::Cairo::RefPtr<Cairo::Context> &crf) override
    {
        Gtk::Allocation allocation = get_allocation();
        allocation.set_x(0);
        allocation.set_y(0);

        // setDrawRectangle will allocate the backbuffer Surface
        if (setDrawRectangle(Cairo::FORMAT_ARGB32, allocation)) {
            setDirty(true);
        }

        if (!isDirty() || !surfaceCreated()) {
            return true;
        }

        Glib::RefPtr<Gtk::StyleContext> style = get_style_context();
        Gtk::Border padding = getPadding(style);  // already scaled
        Cairo::RefPtr<Cairo::Context> cr = getContext();

        if (isDirty()) {
            float h = H_ * 2.f * rtengine::RT_PI / 360.f;
            float a = C_ * std::cos(h);
            float b = C_ * std::sin(h);
            float R, G, B;
            auto iws = rtengine::ICCStore::getInstance()->workingSpaceInverseMatrix("sRGB");
            rtengine::Color::lab2rgb(L_ * 32768.f, a * 42000.f, b * 42000.f, R, G, B, iws);
            cr->set_line_cap(Cairo::LINE_CAP_SQUARE);

            const auto c = [](float v) { return std::pow(v / 65535.f, 1.f/2.2f); };
            cr->set_source_rgb(c(R), c(G), c(B));
            cr->set_operator(Cairo::OPERATOR_OVER);
            cr->paint();
        }

        copySurface(crf);
        return false;
    }

    void on_style_updated() override
    {
        setDirty(true);
        queue_draw();
    }
    
    Gtk::SizeRequestMode get_request_mode_vfunc() const override
    {
        return Gtk::SIZE_REQUEST_HEIGHT_FOR_WIDTH;
    }
        
    void get_preferred_width_vfunc(int &minimum_width, int &natural_width) const override
    {
        Glib::RefPtr<Gtk::StyleContext> style = get_style_context();
        Gtk::Border padding = getPadding(style);  // already scaled
        int s = RTScalable::getScale();
        int p = padding.get_left() + padding.get_right();

        minimum_width = 20 * s + p;
        natural_width = 50 * s + p;
    }
    
    void get_preferred_height_for_width_vfunc(int width, int &minimum_height, int &natural_height) const override
    {
        Glib::RefPtr<Gtk::StyleContext> style = get_style_context();
        Gtk::Border padding = getPadding(style);  // already scaled

        minimum_height = natural_height = width - padding.get_left() - padding.get_right() + padding.get_top() + padding.get_bottom();
    }

    // EditSubscriber interface
    CursorShape getCursor(int objectID) override
    {
        return CSSpotWB;
    }
    
    bool mouseOver(int modifierKey) override
    {
        auto provider = getEditProvider();
        Rectangle *spotRect = static_cast<Rectangle *>(visibleGeometry.at(0));
        spotRect->setXYWH(provider->posImage.x - 16, provider->posImage.y - 16, 32, 32);

        return true;
    }
    
    bool button1Pressed(int modifierKey) override
    { 
        auto provider = getEditProvider();
        EditSubscriber::action = ES_ACTION_NONE;
        sig_spot_requested_.emit(provider->posImage);
        //switchOffEditMode();
        return true;
    }
    
    bool button1Released() override
    {
        EditSubscriber::action = ES_ACTION_NONE;
        return true;
    }
    
    typedef sigc::signal<void, rtengine::Coord> SigSpotRequested;
    SigSpotRequested signal_spot_requested() { return sig_spot_requested_; }

    void activateSpot()
    {
        subscribe();

        int w, h;
        getEditProvider()->getImageSize(w, h);

        // Stick a dummy rectangle over the whole image in mouseOverGeometry.
        // This is to make sure the getCursor() call is fired everywhere.
        Rectangle *imgRect = static_cast<Rectangle *>(mouseOverGeometry.at(0));
        imgRect->setXYWH(0, 0, w, h);
    }

private:
    float L_;
    float C_;
    float H_;

    SigSpotRequested sig_spot_requested_;
};


class DrawnMaskPanel: public MyExpander, public EditSubscriber, public AdjusterListener, public CurveListener {
public:
    DrawnMaskPanel():
        MyExpander(true, M("TP_LABMASKS_DRAWNMASK")),
        EditSubscriber(ET_OBJECTS),
        prev_erase_(false)
    {
        // Editing geometry; create the spot rectangle
        pen_ = new Circle();
        pen_->radiusInImageSpace = true;
        // pen_->setInnerLineColor(255, 255, 255);
        // pen_->setHoverable(false);
        // pen_->setInnerLineWidth(1.f);
        pen_->filled = false;
    
        visibleGeometry.push_back(pen_);

        // Stick a dummy rectangle over the whole image in mouseOverGeometry.
        // This is to make sure the getCursor() call is fired everywhere.
        Rectangle *imgRect = new Rectangle();
        imgRect->filled = true;
        mouseOverGeometry.push_back(imgRect);


        const auto set_btn_style =
            [](Gtk::Button *w) -> void
            {
                w->set_relief(Gtk::RELIEF_NONE);
                w->get_style_context()->add_class(GTK_STYLE_CLASS_FLAT);
                w->set_can_focus(false);
            };

        auto tb = Gtk::manage(new ToolParamBlock());
        add(*tb, false);
        setLevel(1);

        Gtk::VBox *vbg = Gtk::manage(new Gtk::VBox());
        Gtk::HBox *hb = Gtk::manage(new Gtk::HBox());

        Gtk::Button *btncopy = Gtk::manage(new Gtk::Button());
        btncopy->add(*Gtk::manage(new RTImage("copy.png")));
        btncopy->set_tooltip_text(M("TP_LABMASKS_AREA_MASK_COPY_TOOLTIP"));
        btncopy->signal_clicked().connect(sigc::mem_fun(*this, &DrawnMaskPanel::on_copy_pressed));
        set_btn_style(btncopy);
        hb->pack_start(*btncopy, Gtk::PACK_SHRINK);

        Gtk::Button *btnpaste = Gtk::manage(new Gtk::Button());
        btnpaste->add(*Gtk::manage(new RTImage("paste.png")));
        btnpaste->set_tooltip_text(M("TP_LABMASKS_AREA_MASK_PASTE_TOOLTIP"));
        btnpaste->signal_clicked().connect(sigc::mem_fun(*this, &DrawnMaskPanel::on_paste_pressed));
        set_btn_style(btnpaste);
        hb->pack_start(*btnpaste, Gtk::PACK_SHRINK);

        info_ = Gtk::manage(new Gtk::Label(M("TP_LABMASKS_DRAWNMASK_INFO")));
        hb->pack_start(*info_, Gtk::PACK_EXPAND_WIDGET, 4);
        
        const char *img[2] = {
            "area-shape-intersect.png",
            "area-shape-add.png"
        };
        const char *tips[2] = {
            "TP_LABMASKS_DRAWNMASK_INTERSECT_TOOLTIP",
            "TP_LABMASKS_DRAWNMASK_ADD_TOOLTIP"
        };
        for (int i = 0; i < 2; ++i) {
            mode_[i] = Gtk::manage(new Gtk::ToggleButton());
            mode_[i]->add(*Gtk::manage(new RTImage(img[i])));
            modeconn_[i] = mode_[i]->signal_toggled().connect(sigc::bind(sigc::mem_fun(*this, &DrawnMaskPanel::on_mode_changed), i));
            set_btn_style(mode_[i]);
            mode_[i]->set_tooltip_text(M(tips[i]));
        }

        toggle_ = Gtk::manage(new Gtk::ToggleButton());
        toggle_->add(*Gtk::manage(new RTImage("brush.png")));
        set_btn_style(toggle_);
        toggle_->set_tooltip_text(M("TP_LABMASKS_DRAWNMASK_TOGGLE_TIP"));
        toggle_->signal_toggled().connect(sigc::mem_fun(this, &DrawnMaskPanel::on_toggled));

        reset_ = Gtk::manage(new Gtk::Button());
        reset_->add(*Gtk::manage(new RTImage("undo-all.png")));
        set_btn_style(reset_);
        reset_->set_tooltip_text(M("TP_LABMASKS_DRAWNMASK_RESET_TIP"));
        reset_->signal_clicked().connect(sigc::mem_fun(this, &DrawnMaskPanel::on_reset));

        undo_ = Gtk::manage(new Gtk::Button());
        undo_->add(*Gtk::manage(new RTImage("undo.png")));
        set_btn_style(undo_);
        undo_->set_tooltip_text(M("TP_LABMASKS_DRAWNMASK_UNDO_TIP"));
        undo_->signal_clicked().connect(sigc::mem_fun(this, &DrawnMaskPanel::on_undo));
        
        vbg->pack_start(*hb);
        hb = Gtk::manage(new Gtk::HBox());
        
        Gtk::Frame *f = Gtk::manage(new Gtk::Frame(M("TP_LABMASKS_DRAWNMASK_PEN_SETTINGS")));
        Gtk::VBox *vb = Gtk::manage(new Gtk::VBox());
        f->add(*vb);

        radius_ = Gtk::manage(new Adjuster(M("TP_LABMASKS_DRAWNMASK_RADIUS"), 0.1, 100, 0.1, 10));
        radius_->setLogScale(2, 0);
        vb->pack_start(*radius_);

        hardness_ = Gtk::manage(new Adjuster(M("TP_LABMASKS_DRAWNMASK_HARDNESS"), 0, 100, 1, 100));
        vb->pack_start(*hardness_);

        erase_ = Gtk::manage(new Gtk::CheckButton(M("TP_LABMASKS_DRAWNMASK_ERASE")));
        vb->pack_start(*erase_);

        vbg->pack_start(*f);
        hb->pack_start(*vbg);

        Gtk::Grid *grid = Gtk::manage(new Gtk::Grid());
        setExpandAlignProperties(grid, false, true, Gtk::ALIGN_CENTER, Gtk::ALIGN_CENTER);
        grid->attach(*reset_, 0, 1, 1, 1);
        grid->attach(*undo_, 0, 2, 1, 1);
        grid->attach(*toggle_, 0, 3, 1, 1);
        grid->attach(*mode_[0], 0, 4, 1, 1);
        grid->attach(*mode_[1], 0, 5, 1, 1);
        hb->pack_start(*grid, Gtk::PACK_SHRINK, 2);
        
        tb->pack_start(*hb);
        
        feather_ = Gtk::manage(new Adjuster(M("TP_LABMASKS_AREA_FEATHER"), 0, 100, 0.1, 0));
        tb->pack_start(*feather_);
        feather_->setAdjusterListener(this);

        transparency_ = Gtk::manage(new Adjuster(M("TP_LABMASKS_DRAWNMASK_TRANSPARENCY"), 0, 100, 1, 0));
        tb->pack_start(*transparency_);
        transparency_->setAdjusterListener(this);

        smoothness_ = Gtk::manage(new Adjuster(M("TP_LABMASKS_DRAWNMASK_SMOOTHNESS"), 0, 100, 1, 0));
        tb->pack_start(*smoothness_);
        smoothness_->setAdjusterListener(this);

        CurveEditorGroup *cg = Gtk::manage(new CurveEditorGroup(options.lastToneCurvesDir, M("TP_LABMASKS_AREA_CONTRAST")));
        cg->setCurveListener(this);

        contrast_ = static_cast<DiagonalCurveEditor *>(cg->addCurve(CT_Diagonal, ""));
        cg->curveListComplete();
        tb->pack_start(*cg, Gtk::PACK_SHRINK, 2);
        
        signal_enabled_toggled().connect(sigc::mem_fun(*this, &DrawnMaskPanel::on_enabled_toggled));
    }

    ~DrawnMaskPanel()
    {
        for (auto geometry : visibleGeometry) {
            delete geometry;
        }

        for (auto geometry : mouseOverGeometry) {
            delete geometry;
        }        
    }

    void adjusterChanged(Adjuster *a, double newval)
    {
        if (mask_) {
            mask_->feather = feather_->getValue();
            mask_->transparency = transparency_->getValue() / 100.0;
            mask_->smoothness = smoothness_->getValue() / 100.0;
            sig_draw_updated_.emit();
        }
    }

    void curveChanged()
    {
        if (mask_) {
            mask_->contrast = contrast_->getCurve();
            sig_draw_updated_.emit();
        }
    }

    void adjusterAutoToggled(Adjuster *a, bool newval) {}

    // EditSubscriber interface
    CursorShape getCursor(int objectID) override
    {
        return CSEmpty;
    }
    
    bool mouseOver(int modifierKey) override
    {
        if ((modifierKey & GDK_SHIFT_MASK) != prev_erase_) {
            prev_erase_ = (modifierKey & GDK_SHIFT_MASK);
            erase_->set_active(!erase_->get_active());
        }
        update_pen(false);
        return true;
    }
    
    bool button1Pressed(int modifierKey) override
    {
        EditSubscriber::action = ES_ACTION_DRAGGING;
        undo_stack_.push_back(mask_->strokes.size());
        if (modifierKey & GDK_CONTROL_MASK) {
            mask_->strokes.push_back(rtengine::procparams::DrawnMask::Stroke());
        }
        if ((modifierKey & GDK_SHIFT_MASK) != prev_erase_) {
            prev_erase_ = (modifierKey & GDK_SHIFT_MASK);
            erase_->set_active(!erase_->get_active());
        }
        add_stroke(false);
        return true;
    }

    bool drag1(int modifierKey) override
    {
        add_stroke(true);
        update_pen(true);
        return true;
    }
    
    bool button1Released() override
    {
        EditSubscriber::action = ES_ACTION_NONE;
        sig_draw_updated_.emit();
        return true;
    }

    void switchOffEditMode() override
    {
        toggle_->set_active(false);
    }

    void setTargetMask(rtengine::procparams::DrawnMask *mask)
    {
        if (mask != mask_) {
            mask_ = nullptr;
            undo_stack_.clear();
            if (mask) {
                if (toggle_->get_active()) {
                    toggle_->set_active(false);
                }
                info_->set_markup(Glib::ustring::compose(M("TP_LABMASKS_DRAWNMASK_INFO"), mask->strokes.size()));
                setEnabled(mask->enabled);
                feather_->setValue(mask->feather);
                transparency_->setValue(mask->transparency * 100.0);
                smoothness_->setValue(mask->smoothness * 100.0);
                contrast_->setCurve(mask->contrast);
                set_mode(mask->addmode ? 1 : 0);
            }
            mask_ = mask;
        }
    }

    typedef sigc::signal<void> SigDrawUpdated;
    SigDrawUpdated signal_draw_updated() { return sig_draw_updated_; }
    
private:
    void on_toggled()
    {
        if (mask_ && getEnabled() && getEditProvider()) {
            if (toggle_->get_active()) {
                subscribe();
            } else {
                unsubscribe();
            }
        }
    }

    void on_reset()
    {
        if (mask_) {
            mask_->strokes.clear();
            info_->set_markup(Glib::ustring::compose(M("TP_LABMASKS_DRAWNMASK_INFO"), mask_->strokes.size()));
            sig_draw_updated_.emit();
        }
    }

    void on_enabled_toggled()
    {
        if (mask_) {
            mask_->enabled = getEnabled();
            if (!getEnabled()) {
                toggle_->set_active(false);
            }
            sig_draw_updated_.emit();
        }
    }

    void on_undo()
    {
        if (mask_ && !undo_stack_.empty()) {
            size_t n = undo_stack_.back();
            undo_stack_.pop_back();
            mask_->strokes.resize(n);
            info_->set_markup(Glib::ustring::compose(M("TP_LABMASKS_DRAWNMASK_INFO"), n));
            sig_draw_updated_.emit();
        }
    }
    
    void add_stroke(bool dragging)
    {
        auto provider = getEditProvider();
        int w, h;
        provider->getImageSize(w, h);
        mask_->strokes.push_back(rtengine::procparams::DrawnMask::Stroke());
        auto &s = mask_->strokes.back();
        rtengine::Coord p = provider->posImage;
        if (dragging) {
            p += provider->deltaImage;
        }
        s.x = double(p.x) / double(w);
        s.y = double(p.y) / double(h);
        s.radius = radius_->getValue() / 100.0;
        s.erase = erase_->get_active();
        s.hardness = hardness_->getValue() / 100.0;
        info_->set_markup(Glib::ustring::compose(M("TP_LABMASKS_DRAWNMASK_INFO"), mask_->strokes.size()));
    }

    void update_pen(bool dragging)
    {
        auto provider = getEditProvider();
        int w, h;
        provider->getImageSize(w, h);
        pen_->center = provider->posImage;
        if (dragging) {
            pen_->center += provider->deltaImage;
        }
        pen_->radius = radius_->getValue() / 100.0 * std::min(w, h) * 0.25;
    }

    void set_mode(int i)
    {
        for (int j = 0; j < 2; ++j) {
            ConnectionBlocker blocker(modeconn_[j]);
            mode_[j]->set_active(i == j);
        }
    }

    void on_mode_changed(int i)
    {
        if (!mode_[i]->get_active()) {
            ConnectionBlocker blocker(modeconn_[i]);
            mode_[i]->set_active(true);
            return;
        }
        set_mode(i);
        if (mask_) {
            mask_->addmode = mode_[1]->get_active();
            sig_draw_updated_.emit();
        }
    }

    void on_copy_pressed()
    {
        if (mask_) {
            clipboard.setDrawnMask(*mask_);
        }
    }

    void on_paste_pressed()
    {
        if (mask_ && clipboard.hasDrawnMask()) {
            *mask_ = clipboard.getDrawnMask();
            sig_draw_updated_.emit();
        }
    }
    
    rtengine::procparams::DrawnMask *mask_;
    Circle *pen_;
    Gtk::ToggleButton *toggle_;
    Gtk::Label *info_;
    Gtk::Button *reset_;
    Gtk::Button *undo_;
    Adjuster *feather_;
    Adjuster *radius_;
    Adjuster *transparency_;
    Adjuster *smoothness_;
    Adjuster *hardness_;
    Gtk::CheckButton *erase_;
    DiagonalCurveEditor *contrast_;
    Gtk::ToggleButton *mode_[2];
    sigc::connection modeconn_[2];

    SigDrawUpdated sig_draw_updated_;
    bool prev_erase_;
    std::vector<size_t> undo_stack_;
};

bool on_release_event_ignore(GdkEventButton *event)
{
    return true;
}

} // namespace


LabMasksPanel::LabMasksPanel(LabMasksContentProvider *cp):
    Gtk::VBox(),
    cp_(cp),
    masks_(),
    selected_(0),
    area_shape_index_(0),
    listEdited(false),
    adl_(nullptr),
    deltaE_provider_(nullptr)
{
    Gtk::Widget *child = cp_->getWidget();
    cp_->getEvents(EvMaskList, EvParametricMask, EvHMask, EvCMask, EvLMask, EvMaskBlur, EvShowMask, EvAreaMask, EvDeltaEMask, EvContrastThresholdMask, EvDrawnMask);
    EvAreaMaskVoid = ProcEventMapper::getInstance()->newEvent(M_VOID, EvAreaMask.get_message());
    EvDeltaEMaskVoid = ProcEventMapper::getInstance()->newEvent(M_VOID, EvDeltaEMask.get_message());
    EvMaskName = ProcEventMapper::getInstance()->newEvent(M_VOID, "HISTORY_MSG_LABMASKS_MASK_NAME");
    
    CurveListener::setMulti(true);
    
    const auto add_button =
        [](Gtk::Button *btn, Gtk::Box *box, int h=20) -> void
        {
            setExpandAlignProperties(btn, false, false, Gtk::ALIGN_CENTER, Gtk::ALIGN_START);
            btn->set_relief(Gtk::RELIEF_NONE);
            btn->get_style_context()->add_class(GTK_STYLE_CLASS_FLAT);
            btn->set_can_focus(false);
            if (h > 0) {
                btn->set_size_request(-1, h);
            }
            box->pack_start(*btn, false, false);
        };
    
    int n = cp_->getColumnCount();

    list_model_columns_.reset(new ListColumns(n));
    list_model_ = Gtk::ListStore::create(*list_model_columns_);

    list_enabled_renderer_.property_mode() = Gtk::CELL_RENDERER_MODE_ACTIVATABLE;
    list_enabled_renderer_.signal_toggled().connect(sigc::mem_fun(this, &LabMasksPanel::onListEnabledToggled));
    list_enabled_column_.pack_start(list_enabled_renderer_);
    list_enabled_column_.set_cell_data_func(list_enabled_renderer_, sigc::mem_fun(this, &LabMasksPanel::setListEnabled));
    {
        auto img = Gtk::manage(new RTImage("power-on-small-faded.png"));
        img->show();
        list_enabled_column_.set_widget(*img);
    }
    list = Gtk::manage(new Gtk::TreeView());
    list->set_model(list_model_);
    list->set_size_request(-1, 150);
    list->set_can_focus(false);
    list->append_column("#", list_model_columns_->id);
    for (int i = 0; i < n; ++i) {
        int col = list->append_column(cp_->getColumnHeader(i), list_model_columns_->cols[i]);
        list->get_column(col-1)->set_expand(true);
    }
    /*int col =*/ list->append_column(M("TP_LABMASKS_MASK"), list_model_columns_->mask);
    //list->get_column(col-1)->set_expand(true);
    list->append_column(list_enabled_column_);
    list->set_activate_on_single_click(true);
    selectionConn = list->get_selection()->signal_changed().connect(sigc::mem_fun(this, &LabMasksPanel::onSelectionChanged));
    Gtk::HBox *hb = Gtk::manage(new Gtk::HBox());
    Gtk::ScrolledWindow *scroll = Gtk::manage(new Gtk::ScrolledWindow());
    scroll->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_NEVER);
    scroll->add(*list);
    hb->pack_start(*scroll, Gtk::PACK_EXPAND_WIDGET);
    Gtk::VBox *vb = Gtk::manage(new Gtk::VBox());
    reset = Gtk::manage(new Gtk::Button());
    reset->add(*Gtk::manage(new RTImage("undo-small.png")));
    reset->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onResetPressed));
    add_button(reset, vb);
    add = Gtk::manage(new Gtk::Button());
    add->add(*Gtk::manage(new RTImage("add-small.png")));
    add->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onAddPressed));
    add_button(add, vb);
    remove = Gtk::manage(new Gtk::Button());
    remove->add(*Gtk::manage(new RTImage("remove-small.png")));
    remove->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onRemovePressed));
    add_button(remove, vb);
    up = Gtk::manage(new Gtk::Button());
    up->add(*Gtk::manage(new RTImage("arrow-up-small.png")));
    up->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onUpPressed));
    add_button(up, vb);
    down = Gtk::manage(new Gtk::Button());
    down->add(*Gtk::manage(new RTImage("arrow-down-small.png")));
    down->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onDownPressed));
    add_button(down, vb);
    copy = Gtk::manage(new Gtk::Button());
    copy->add(*Gtk::manage(new RTImage("copy-small.png")));
    copy->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onCopyPressed));
    add_button(copy, vb);
    hb->pack_start(*vb, Gtk::PACK_SHRINK);
    pack_start(*hb, true, true);

    pack_start(*Gtk::manage(new Gtk::HSeparator()));
    pack_start(*child);
    pack_start(*Gtk::manage(new Gtk::HSeparator()));

    Gtk::VBox *mask_box = Gtk::manage(new Gtk::VBox());

    showMask = Gtk::manage(new Gtk::CheckButton(M("TP_LABMASKS_SHOW")));
    showMask->signal_toggled().connect(sigc::mem_fun(*this, &LabMasksPanel::onShowMaskChanged));
    hb = Gtk::manage(new Gtk::HBox());
    hb->pack_start(*showMask);

    maskInverted = Gtk::manage(new Gtk::CheckButton(M("TP_LABMASKS_INVERTED")));
    maskInverted->signal_toggled().connect(sigc::mem_fun(*this, &LabMasksPanel::onMaskInvertedChanged));
    hb->pack_start(*maskInverted);
    mask_box->pack_start(*hb);

    parametricMask = Gtk::manage(new MyExpander(true, M("TP_LABMASKS_PARAMETRIC")));
    ToolParamBlock *tb = Gtk::manage(new ToolParamBlock());
    parametricMask->add(*tb, false);
    parametricMask->setLevel(1);
    parametricMask->signal_enabled_toggled().connect(sigc::mem_fun(*this, &LabMasksPanel::onParametricMaskEnableToggled));
        
    maskEditorGroup = Gtk::manage(new CurveEditorGroup(options.lastColorToningCurvesDir, "", 0.7));
    maskEditorGroup->setCurveListener(this);

    rtengine::Mask default_params;

    EditUniqueID eh, ec, el, ede;
    cp_->getEditIDs(eh, ec, el, ede);

    hueMask = static_cast<FlatCurveEditor *>(maskEditorGroup->addCurve(CT_Flat, M("TP_LABMASKS_HUE"), nullptr, false, true));
    hueMask->setIdentityValue(0.);
    hueMask->setResetCurve(FlatCurveType(default_params.parametricMask.hue[0]), default_params.parametricMask.hue);
    hueMask->setCurveColorProvider(this, ID_HUE_MASK);
    hueMask->setBottomBarColorProvider(this, ID_HUE_MASK);
    hueMask->setEditID(eh, BT_SINGLEPLANE_FLOAT);

    chromaticityMask = static_cast<FlatCurveEditor *>(maskEditorGroup->addCurve(CT_Flat, M("TP_LABMASKS_CHROMATICITY"), nullptr, false, false));
    chromaticityMask->setIdentityValue(0.);
    chromaticityMask->setResetCurve(FlatCurveType(default_params.parametricMask.chromaticity[0]), default_params.parametricMask.chromaticity);
    chromaticityMask->setBottomBarColorProvider(this, ID_HUE_MASK+1);
    chromaticityMask->setEditID(ec, BT_SINGLEPLANE_FLOAT);

    lightnessMask = static_cast<FlatCurveEditor *>(maskEditorGroup->addCurve(CT_Flat, M("TP_LABMASKS_LIGHTNESS"), nullptr, false, false));
    lightnessMask->setIdentityValue(0.);
    lightnessMask->setResetCurve(FlatCurveType(default_params.parametricMask.lightness[0]), default_params.parametricMask.lightness);
    std::vector<GradientMilestone> milestones = {
        GradientMilestone(0., 0., 0., 0.),
        GradientMilestone(1., 1., 1., 1.)
    };
    lightnessMask->setBottomBarBgGradient(milestones);
    lightnessMask->setEditID(el, BT_SINGLEPLANE_FLOAT);

    maskEditorGroup->curveListComplete();
    maskEditorGroup->show();
    tb->pack_start(*maskEditorGroup);//, Gtk::PACK_SHRINK, 4);

    contrastThreshold = Gtk::manage(new Adjuster(M("TP_LABMASKS_CONTRASTTHRESHOLDMASK"), -150, 150, 1, 0));
    contrastThreshold->setAdjusterListener(this);
    //mask_box->pack_start(*contrastThreshold);
    tb->pack_start(*contrastThreshold);

    maskBlur = Gtk::manage(new Adjuster(M("TP_LABMASKS_BLUR"), -10, 500, 0.1, 0));
    maskBlur->setLogScale(10, -10);
    maskBlur->setAdjusterListener(this);
    //mask_box->pack_start(*maskBlur);
    tb->pack_start(*maskBlur);

    mask_box->pack_start(*parametricMask);
    
    //-------------------------------------------------------------------------
    deltaEMask = Gtk::manage(new MyExpander(true, M("TP_LABMASKS_DELTAE")));
    deltaEMask->getLabelWidget()->set_tooltip_text(M("TP_LABMASKS_DELTAE_TOOLTIP"));
    DeltaEArea *dE_area = Gtk::manage(new DeltaEArea());
    deltaEColor = dE_area;
    const auto deltaAdj =
        [&](const Glib::ustring &lbl, double vmin, double vmax, double vdflt, unsigned int prec, int i) -> ThresholdAdjuster *
        {
            ThresholdAdjuster *a = Gtk::manage(new ThresholdAdjuster(lbl, 0, 100, 100, M("TP_LABMASKS_DELTAE_CHAN_WEIGHT"), 0, vmin, vmax, vdflt, lbl, prec, nullptr));
            a->setBgColorProvider(this, i);
            a->setAdjusterListener(this);
            a->setUpdatePolicy(RTUP_DYNAMIC);
            return a;
        };
    deltaEL = deltaAdj(M("TP_LABMASKS_DELTAE_L"), 0, 1, 0, 3, ID_HUE_MASK+2);
    deltaEC = deltaAdj(M("TP_LABMASKS_DELTAE_C"), 0, 1, 0, 3, ID_HUE_MASK+3);
    deltaEH = deltaAdj(M("TP_LABMASKS_DELTAE_H"), 0, 360, 0, 1, ID_HUE_MASK+4);
    deltaERange = Gtk::manage(new Adjuster(M("TP_LABMASKS_DELTAE_RANGE"), 1, 100, 0.1, 1));
    deltaERange->setLogScale(10.f, 3.f, true);
    deltaERange->setAdjusterListener(this);
    deltaEDecay = Gtk::manage(new Adjuster(M("TP_LABMASKS_DELTAE_DECAY"), 0, 100, 1, 0));
    deltaEDecay->setLogScale(10.f, 10.f, true);
    deltaEDecay->setAdjusterListener(this);
    vb = Gtk::manage(new Gtk::VBox());
    hb = Gtk::manage(new Gtk::HBox());
    vb->pack_start(*deltaERange);
    vb->pack_start(*deltaEDecay);
    hb->pack_start(*vb);
    vb = Gtk::manage(new Gtk::VBox());
    deltaEPick = Gtk::manage(new Gtk::Button(M("TP_LABMASKS_DELTAE_PICK")));
    vb->pack_start(*deltaEPick);
    vb->pack_start(*deltaEColor, Gtk::PACK_EXPAND_WIDGET, 4);
    hb->pack_start(*vb, Gtk::PACK_SHRINK, 4);
    tb = Gtk::manage(new ToolParamBlock());
    tb->pack_start(*deltaEL);
    tb->pack_start(*deltaEC);
    tb->pack_start(*deltaEH);
    tb->pack_start(*hb);
    deltaEMask->add(*tb, false);
    deltaEMask->setLevel(1);
    mask_box->pack_start(*deltaEMask);
    deltaEMask->signal_enabled_toggled().connect(sigc::mem_fun(*this, &LabMasksPanel::onDeltaEMaskEnableToggled));
    dE_area->setEditID(ede, BT_SINGLEPLANE_FLOAT);
    dE_area->signal_spot_requested().connect(sigc::mem_fun(*this, &LabMasksPanel::onDeltaESpotRequested));
    deltaEPick->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onDeltaEPickClicked));
    //-------------------------------------------------------------------------
        
    areaMask = Gtk::manage(new MyExpander(true, M("TP_LABMASKS_AREA")));
    areaMask->signal_enabled_toggled().connect(sigc::mem_fun(*this, &LabMasksPanel::onAreaMaskEnableToggled));
    ToolParamBlock *area = Gtk::manage(new ToolParamBlock());
    hb = Gtk::manage(new Gtk::HBox());
    areaMaskButtonsHb = hb;

    areaMaskCopy = Gtk::manage(new Gtk::Button());
    areaMaskCopy->add(*Gtk::manage(new RTImage("copy.png")));
    areaMaskCopy->set_tooltip_text(M("TP_LABMASKS_AREA_MASK_COPY_TOOLTIP"));
    areaMaskCopy->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onAreaMaskCopyPressed));
    add_button(areaMaskCopy, hb, 24);
    
    areaMaskPaste = Gtk::manage(new Gtk::Button());
    areaMaskPaste->add(*Gtk::manage(new RTImage("paste.png")));
    areaMaskPaste->set_tooltip_text(M("TP_LABMASKS_AREA_MASK_PASTE_TOOLTIP"));
    areaMaskPaste->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onAreaMaskPastePressed));
    add_button(areaMaskPaste, hb, 24);

    hb->pack_start(*Gtk::manage(new Gtk::Label("")), Gtk::PACK_EXPAND_WIDGET);

    areaMaskDrawAdd = new Gtk::Button();
    areaMaskDrawAdd->add(*Gtk::manage(new RTImage("area-shape-draw-add.png")));
    areaMaskDrawAdd->set_tooltip_text(M("TP_LABMASKS_AREA_MASK_DRAW_ADD_TOOLTIP"));
    areaMaskDrawAdd->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onAreaMaskDrawAddPressed));
    add_button(areaMaskDrawAdd, hb, 24);

    areaMaskDraw = new Gtk::ToggleButton();
    areaMaskDraw->add(*Gtk::manage(new RTImage("area-shape-draw.png")));
    areaMaskDraw->set_tooltip_text(M("TP_LABMASKS_AREA_MASK_DRAW_TOOLTIP"));
    areaMaskDrawConn = areaMaskDraw->signal_toggled().connect(sigc::mem_fun(*this, &LabMasksPanel::onAreaMaskDrawChanged));
    add_button(areaMaskDraw, hb, 24);
    
    areaMaskToggle = new Gtk::ToggleButton();
    areaMaskToggle->add(*Gtk::manage(new RTImage("crosshair-adjust.png")));
    areaMaskToggle->set_tooltip_text(M("TP_LABMASKS_AREA_MASK_TOGGLE_TOOLTIP"));
    areaMaskToggle->signal_toggled().connect(sigc::mem_fun(*this, &LabMasksPanel::onAreaMaskToggleChanged));
    add_button(areaMaskToggle, hb, 24);
    area->pack_start(*hb);

    const auto add_adjuster =
        [&](Adjuster *a, Gtk::Box *b) -> void
        {
            areaMaskAdjusters.push_back(a);
            a->setAdjusterListener(this);
            if (b) {
                b->pack_start(*a);
            }
        };
    
    areaMaskFeather = Gtk::manage(new Adjuster(M("TP_LABMASKS_AREA_FEATHER"), 0, 100, 0.1, 0));
    add_adjuster(areaMaskFeather, area);

    areaMaskBlur = Gtk::manage(new Adjuster(M("TP_LABMASKS_BLUR"), 0, 500, 0.1, 0));
    add_adjuster(areaMaskBlur, area);
    
    CurveEditorGroup *cg = Gtk::manage(new CurveEditorGroup(options.lastToneCurvesDir, M("TP_LABMASKS_AREA_CONTRAST")));
    cg->setCurveListener(this);

    areaMaskContrast = static_cast<DiagonalCurveEditor *>(cg->addCurve(CT_Diagonal, ""));
    cg->curveListComplete();
    if (area) {
        area->pack_start(*cg, Gtk::PACK_SHRINK, 2);
    }
    
    Gtk::Frame *areaFrame = Gtk::manage(new Gtk::Frame(M("TP_LABMASKS_AREA_SHAPES")));
    areaFrame->set_label_align(0.025, 0.5);
    areaFrame->set_border_width(4);

    areaMaskShapes = Gtk::manage(new Gtk::ListViewText(2));
    areaMaskShapes->set_size_request(-1, 80);
    areaMaskShapes->set_can_focus(false);
    areaMaskShapes->set_column_title(0, "#");
    areaMaskShapes->set_column_title(1, M("TP_LABMASKS_AREA_SHAPE"));
    areaMaskShapes->set_activate_on_single_click(true);
    shapeSelectionConn = areaMaskShapes->get_selection()->signal_changed().connect(sigc::mem_fun(this, &LabMasksPanel::onAreaShapeSelectionChanged));
    hb = Gtk::manage(new Gtk::HBox());
    scroll = Gtk::manage(new Gtk::ScrolledWindow());
    scroll->set_policy(Gtk::POLICY_AUTOMATIC, Gtk::POLICY_NEVER);
    scroll->add(*areaMaskShapes);
    hb->pack_start(*scroll, Gtk::PACK_EXPAND_WIDGET);
    vb = Gtk::manage(new Gtk::VBox());
    areaMaskReset = Gtk::manage(new Gtk::Button());
    areaMaskReset->add(*Gtk::manage(new RTImage("undo-small.png")));
    areaMaskReset->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onAreaShapeResetPressed));
    add_button(areaMaskReset, vb);
    areaMaskAdd = Gtk::manage(new Gtk::Button());
    areaMaskAdd->add(*Gtk::manage(new RTImage("add-small.png")));
    areaMaskAdd->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onAreaShapeAddPressed));
    add_button(areaMaskAdd, vb);
    areaMaskRemove = Gtk::manage(new Gtk::Button());
    areaMaskRemove->add(*Gtk::manage(new RTImage("remove-small.png")));
    areaMaskRemove->signal_clicked().connect(sigc::mem_fun(*this, &LabMasksPanel::onAreaShapeRemovePressed));
    add_button(areaMaskRemove, vb);
    hb->pack_start(*vb, Gtk::PACK_SHRINK);
    
    vb = Gtk::manage(new Gtk::VBox());
    vb->set_spacing(2);
    vb->pack_start(*hb);

    hb = Gtk::manage(new Gtk::HBox());
    hb->pack_start(*Gtk::manage(new Gtk::Label(M("TP_LABMASKS_AREA_SHAPE_MODE") + ":")), Gtk::PACK_SHRINK, 4);
    hb->pack_start(*Gtk::manage(new Gtk::Label("")), Gtk::PACK_EXPAND_WIDGET);
    const char *img[3] = {
        "area-shape-add.png",
        "area-shape-subtract.png",
        "area-shape-intersect.png"
    };
    const char *tips[3] = {
        "TP_LABMASKS_AREA_SHAPE_MODE_ADD",
        "TP_LABMASKS_AREA_SHAPE_MODE_SUBTRACT",
        "TP_LABMASKS_AREA_SHAPE_MODE_INTERSECT"
    };
    for (int i = 0; i < 3; ++i) {
        areaMaskMode[i] = Gtk::manage(new Gtk::ToggleButton());
        areaMaskMode[i]->add(*Gtk::manage(new RTImage(img[i])));
        areaMaskMode[i]->set_tooltip_text(M(tips[i]));
        areaMaskModeConn[i] = areaMaskMode[i]->signal_toggled().connect(sigc::bind(sigc::mem_fun(*this, &LabMasksPanel::onAreaShapeModeChanged), i));
        add_button(areaMaskMode[i], hb, 24);
    }
    vb->pack_start(*hb);
    
    areaMaskRoundness = Gtk::manage(new Adjuster(M("TP_LABMASKS_AREA_ROUNDNESS"), 0, 100, 0.1, 0));
    add_adjuster(areaMaskRoundness, vb);
    areaMaskX = Gtk::manage(new Adjuster(M("TP_LABMASKS_AREA_X"), -100, 100, 0.1, 0));
    add_adjuster(areaMaskX, vb);
    areaMaskY = Gtk::manage(new Adjuster(M("TP_LABMASKS_AREA_Y"), -100, 100, 0.1, 0));
    add_adjuster(areaMaskY, vb);
    areaMaskWidth = Gtk::manage(new Adjuster(M("TP_LABMASKS_AREA_WIDTH"), 1, 200, 0.1, 100));
    areaMaskWidth->setLogScale(10, 1);
    add_adjuster(areaMaskWidth, vb);
    areaMaskHeight = Gtk::manage(new Adjuster(M("TP_LABMASKS_AREA_HEIGHT"), 1, 200, 0.1, 100));
    areaMaskHeight->setLogScale(10, 1);
    add_adjuster(areaMaskHeight, vb);
    areaMaskAngle = Gtk::manage(new Adjuster(M("TP_LABMASKS_AREA_ANGLE"), 0, 180, 0.1, 0));
    add_adjuster(areaMaskAngle, vb);

    areaFrame->add(*vb);
    vb = Gtk::manage(new Gtk::VBox());
    vb->pack_start(*area);
    vb->pack_start(*areaFrame);
    areaMask->add(*vb, false);
    areaMask->setLevel(1);
    mask_box->pack_start(*areaMask);

    drawnMask = Gtk::manage(new DrawnMaskPanel());
    mask_box->pack_start(*drawnMask);
    static_cast<DrawnMaskPanel *>(drawnMask)->signal_draw_updated().connect(sigc::mem_fun(this, &LabMasksPanel::onDrawnMaskUpdated));
    
    mask_box->set_border_width(4);

    MyExpander *mask_exp = nullptr;
    {
        Gtk::HBox *hb = Gtk::manage(new Gtk::HBox());
        Gtk::Label *l = Gtk::manage(new Gtk::Label());
        l->set_markup("<b>" + M("TP_LABMASKS_MASK") + "</b>");
        l->set_alignment(Gtk::ALIGN_START, Gtk::ALIGN_CENTER);
        hb->pack_start(*l);
        Gtk::Entry *e = Gtk::manage(new Gtk::Entry());
        hb->pack_start(*e, Gtk::PACK_EXPAND_WIDGET, 2);
        e->get_style_context()->add_class(GTK_STYLE_CLASS_FLAT);
        e->set_alignment(Gtk::ALIGN_END);
        e->set_placeholder_text("(" + M("TP_LABMASKS_MASK_UNNAMED") + ")");
        mask_exp = Gtk::manage(new MyExpander(false, hb));
        e->signal_button_release_event().connect(&on_release_event_ignore);
        e->add_events(Gdk::FOCUS_CHANGE_MASK);
        e->signal_focus_out_event().connect(sigc::mem_fun(*this, &LabMasksPanel::onMaskNameFocusOut));
        maskName = e;
    }
    mask_exp->add(*mask_box, false);
    mask_exp->setLevel(2);
    pack_start(*mask_exp);

    mask_exp_ = mask_exp;
    first_mask_exp_ = true;
    mask_exp_->signal_button_release_event().connect_notify(sigc::mem_fun(this, &LabMasksPanel::onMaskFold));
        
    maskBlur->delay = options.adjusterMaxDelay;
}


LabMasksPanel::~LabMasksPanel()
{
    delete areaMaskToggle;
}


ToolPanelListener *LabMasksPanel::getListener()
{
    if (listenerDisabled.empty() || !listenerDisabled.back()) {
        return cp_->listener();
    }
    return nullptr;
}


void LabMasksPanel::disableListener()
{
    listenerDisabled.push_back(true);
    //static_cast<DrawnMaskPanel *>(drawnMask)->setTargetMask(nullptr);
}


void LabMasksPanel::enableListener()
{
    listenerDisabled.pop_back();
    // if (listenerDisabled.empty() || !listenerDisabled.back()) {
    //     static_cast<DrawnMaskPanel *>(drawnMask)->setTargetMask(&masks_[selected_].drawnMask);
    // }
}


void LabMasksPanel::onSelectionChanged()
{
    auto s = list->get_selection()->get_selected_rows();
    if (!s.empty()) {
        int idx = s[0][0];
        // update the selected values
        maskGet(selected_);
        cp_->selectionChanging(selected_);
        selected_ = idx;
        area_shape_index_ = 0;
        maskShow(selected_);
        if (showMask->get_active()) {
            onShowMaskChanged();
        }
    }
}


void LabMasksPanel::maskGet(int idx)
{
    if (idx < 0 || size_t(idx) >= masks_.size()) {
        return;
    }
    
    auto &r = masks_[idx];
    r.parametricMask.enabled = parametricMask->getEnabled();
    r.parametricMask.hue = hueMask->getCurve();
    r.parametricMask.chromaticity = chromaticityMask->getCurve();
    r.parametricMask.lightness = lightnessMask->getCurve();
    r.parametricMask.blur = maskBlur->getValue();
    r.parametricMask.contrastThreshold = contrastThreshold->getValue();
    r.inverted = maskInverted->get_active();
    r.areaMask.enabled = areaMask->getEnabled();
    r.areaMask.feather = areaMaskFeather->getValue();
    r.areaMask.blur = areaMaskBlur->getValue();
    r.areaMask.contrast = areaMaskContrast->getCurve();
    if (area_shape_index_ < r.areaMask.shapes.size()) {
        auto &a = r.areaMask.shapes[area_shape_index_];
        a.x = areaMaskX->getValue();
        a.y = areaMaskY->getValue();
        a.width = areaMaskWidth->getValue();
        a.height = areaMaskHeight->getValue();
        a.angle = areaMaskAngle->getValue();
        a.roundness = areaMaskRoundness->getValue();
        a.mode = Shape::Mode(getAreaShapeMode());
    }
    r.deltaEMask.enabled = deltaEMask->getEnabled();
    double b, t;
    deltaEL->getValue(b, t);
    r.deltaEMask.L = t;
    r.deltaEMask.weight_L = b;
    deltaEC->getValue(b, t);
    r.deltaEMask.C = t;
    r.deltaEMask.weight_C = b;
    deltaEH->getValue(b, t);
    r.deltaEMask.H = t;
    r.deltaEMask.weight_H = b;
    r.deltaEMask.range = deltaERange->getValue();
    r.deltaEMask.decay = deltaEDecay->getValue();
    r.name = maskName->get_text();
}


void LabMasksPanel::onAddPressed()
{
    if (!cp_->addPressed()) {
        return;
    }

    listEdited = true;
    selected_ = masks_.size();
    masks_.push_back(rtengine::procparams::Mask());
    masks_.back().areaMask.shapes = {defaultAreaShape};
    populateList();
    area_shape_index_ = 0;
    maskShow(selected_);

    auto l = getListener();
    if (l) {
        l->panelChanged(EvMaskList, M("HISTORY_CHANGED"));
    }
}


void LabMasksPanel::onRemovePressed()
{
    if (list_model_->children().size() <= 1 || !cp_->removePressed(selected_)) {
        return;
    }
    
    listEdited = true;
    masks_.erase(masks_.begin() + selected_);
    selected_ = rtengine::LIM(selected_-1, 0u, unsigned(masks_.size()-1));
    populateList();
    area_shape_index_ = 0;
    maskShow(selected_);

    auto l = getListener();
    if (l) {
        l->panelChanged(EvMaskList, M("HISTORY_CHANGED"));
    }
}


void LabMasksPanel::onUpPressed()
{
    if (selected_ > 0 && cp_->moveUpPressed(selected_)) {
        listEdited = true;
        
        auto r = masks_[selected_];
        masks_.erase(masks_.begin() + selected_);
        --selected_;
        masks_.insert(masks_.begin() + selected_, r);
        populateList();

        auto l = getListener();
        if (l) {
            l->panelChanged(EvMaskList, M("HISTORY_CHANGED"));
        }
    }
}


void LabMasksPanel::onDownPressed()
{
    if (selected_ < masks_.size()-1 && cp_->moveDownPressed(selected_)) {
        listEdited = true;
        
        auto r = masks_[selected_];
        masks_.erase(masks_.begin() + selected_);
        ++selected_;
        masks_.insert(masks_.begin() + selected_, r);
        populateList();

        auto l = getListener();
        if (l) {
            l->panelChanged(EvMaskList, M("HISTORY_CHANGED"));
        }
    }
}


void LabMasksPanel::onCopyPressed()
{
    if (selected_ < masks_.size() && cp_->copyPressed(selected_)) {
        listEdited = true;
        
        auto r = masks_[selected_];
        masks_.push_back(r);
        selected_ = masks_.size()-1;
        populateList();

        auto l = getListener();
        if (l) {
            l->panelChanged(EvMaskList, M("HISTORY_CHANGED"));
        }
    }
}


void LabMasksPanel::onResetPressed()
{
    if (cp_->resetPressed()) {
        listEdited = true;
        auto l = getListener();
        if (l) {
            l->panelChanged(EvMaskList, M("HISTORY_CHANGED"));
        }
    }
}


void LabMasksPanel::onShowMaskChanged()
{
    auto l = getListener();
    if (l) {
        l->panelChanged(EvShowMask, showMask->get_active() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}


void LabMasksPanel::populateList()
{
    ConnectionBlocker b(selectionConn);
    list_model_->clear();
    rtengine::procparams::Mask dflt;

    int n = cp_->getColumnCount();
    for (size_t i = 0; i < masks_.size(); ++i) {
        auto &r = masks_[i];
        auto j = list_model_->append();
        auto row = *j;
        row[list_model_columns_->enabled] = r.enabled;
        row[list_model_columns_->id] = i+1;
        for (int c = 0; c < n; ++c) {
            row[list_model_columns_->cols[c]] = cp_->getColumnContent(c, i);
        }
        if (!r.name.empty()) {
            row[list_model_columns_->mask] = r.name;            
        } else {
            Glib::ustring am("");
            if (r.areaMask.enabled && !r.areaMask.isTrivial()) {
                am = Glib::ustring::compose("\n%1 shape%2", r.areaMask.shapes.size(), r.areaMask.shapes.size() > 1 ? "s" : "");
            }
            if (!r.drawnMask.isTrivial()) {
                if (am.empty()) {
                    am = "\n";
                } else {
                    am += " ";
                }
                am += Glib::ustring::compose("%1 stroke%2", r.drawnMask.strokes.size(), r.drawnMask.strokes.size() == 1 ? "" : "s");
            }
            row[list_model_columns_->mask] = 
                Glib::ustring::compose(
                    "%1%2%3%4%7%5%6",
                    hasMask(r, dflt.parametricMask.hue, r.parametricMask.hue) ? "H" : "",
                    hasMask(r, dflt.parametricMask.chromaticity, r.parametricMask.chromaticity) ? "C" : "",
                    hasMask(r, dflt.parametricMask.lightness, r.parametricMask.lightness) ? "L" : "",
                    r.deltaEMask.enabled ? "ΔE" : "",
                    r.parametricMask.blur ? Glib::ustring::compose(" b=%1", r.parametricMask.blur) : "",
                    am, (r.parametricMask.enabled && r.parametricMask.contrastThreshold) ? Glib::ustring::compose(" c=%1", r.parametricMask.contrastThreshold) : "");
        }
    }
}


void LabMasksPanel::maskShow(int idx, bool list_only, bool unsub)
{
    disableListener();
    rtengine::procparams::Mask dflt;
    auto &r = masks_[idx];
    if (!list_only) {
        cp_->selectionChanged(idx);
        parametricMask->setEnabled(r.parametricMask.enabled);
        hueMask->setCurve(r.parametricMask.hue);
        chromaticityMask->setCurve(r.parametricMask.chromaticity);
        lightnessMask->setCurve(r.parametricMask.lightness);
        contrastThreshold->setValue(r.parametricMask.contrastThreshold);
        maskBlur->setValue(r.parametricMask.blur);
        maskInverted->set_active(r.inverted);
        maskName->set_text(r.name);

        if (unsub && isCurrentSubscriber()) {
            unsubscribe();
        }
        areaMaskToggle->set_active(false);
        areaMask->setEnabled(r.areaMask.enabled);
        areaMaskFeather->setValue(r.areaMask.feather);
        areaMaskBlur->setValue(r.areaMask.blur);
        areaMaskContrast->setCurve(r.areaMask.contrast);
        if (area_shape_index_ < r.areaMask.shapes.size()) {
            auto &a = r.areaMask.shapes[area_shape_index_];
            areaMaskX->setValue(a.x);
            areaMaskY->setValue(a.y);
            areaMaskWidth->setValue(a.width);
            areaMaskHeight->setValue(a.height);
            areaMaskAngle->setValue(a.angle);
            areaMaskRoundness->setValue(a.roundness);
            toggleAreaShapeMode(int(a.mode));
        }
        populateShapeList(idx, area_shape_index_);

        deltaEMask->setEnabled(r.deltaEMask.enabled);
        deltaEL->setValue(r.deltaEMask.weight_L, r.deltaEMask.L);
        deltaEC->setValue(r.deltaEMask.weight_C, r.deltaEMask.C);
        deltaEH->setValue(r.deltaEMask.weight_H, r.deltaEMask.H);
        deltaERange->setValue(r.deltaEMask.range);
        deltaEDecay->setValue(r.deltaEMask.decay);
        static_cast<DeltaEArea *>(deltaEColor)->setColor(r.deltaEMask.L, r.deltaEMask.C, r.deltaEMask.H);
        
        updateAreaMask(false);
    }
    static_cast<DrawnMaskPanel *>(drawnMask)->setTargetMask(&r.drawnMask);

    int n = cp_->getColumnCount();
    auto row = list_model_->children()[idx];
    row[list_model_columns_->enabled] = r.enabled;
    for (int c = 0; c < n; ++c) {
        row[list_model_columns_->cols[c]] = cp_->getColumnContent(c, idx);
    }
    if (!r.name.empty()) {
        row[list_model_columns_->mask] = r.name;
    } else {
        Glib::ustring am("");
        if (r.areaMask.enabled && !r.areaMask.isTrivial()) {
            am = Glib::ustring::compose("\n%1 shape%2", r.areaMask.shapes.size(),
                                        r.areaMask.shapes.size() > 1 ? "s" : "");
        }
        if (!r.drawnMask.isTrivial()) {
            if (am.empty()) {
                am = "\n";
            } else {
                am += " ";
            }
            am += Glib::ustring::compose("%1 stroke%2", r.drawnMask.strokes.size(), r.drawnMask.strokes.size() == 1 ? "" : "s");
        }
        row[list_model_columns_->mask] = 
            Glib::ustring::compose(
                "%1%2%3%4%7%5%6",
                hasMask(r, dflt.parametricMask.hue, r.parametricMask.hue) ? "H" : "",
                hasMask(r, dflt.parametricMask.chromaticity, r.parametricMask.chromaticity) ? "C" : "",
                hasMask(r, dflt.parametricMask.lightness, r.parametricMask.lightness) ? "L" : "",
                r.deltaEMask.enabled ? "ΔE" : "",
                r.parametricMask.blur ? Glib::ustring::compose(" b=%1", r.parametricMask.blur) : "", am,
                (r.parametricMask.enabled && r.parametricMask.contrastThreshold) ? Glib::ustring::compose(" c=%1", r.parametricMask.contrastThreshold) : "");
    }
    Gtk::TreePath pth;
    pth.push_back(idx);
    list->get_selection()->select(pth);
    enableListener();
}


void LabMasksPanel::setEditProvider(EditDataProvider *provider)
{
    hueMask->setEditProvider(provider);
    chromaticityMask->setEditProvider(provider);
    lightnessMask->setEditProvider(provider);
    static_cast<DeltaEArea *>(deltaEColor)->setEditProvider(provider);
    AreaMask::setEditProvider(provider);
    static_cast<DrawnMaskPanel *>(drawnMask)->setEditProvider(provider);
}


void LabMasksPanel::onAreaMaskToggleChanged()
{
    if (areaMaskToggle->get_active()) {
        areaMaskDraw->set_active(false);
        subscribe();
    } else {
        unsubscribe();
    }
}


void LabMasksPanel::onMaskInvertedChanged()
{
    auto l = getListener();
    if (l) {
        l->panelChanged(EvAreaMask, M("GENERAL_CHANGED"));
    }
}


void LabMasksPanel::updateAreaMask(bool from_mask)
{
    disableListener();
    if (from_mask) {
        areaMaskX->setValue(center_x_);
        areaMaskY->setValue(center_y_);
        areaMaskWidth->setValue(width_);
        areaMaskHeight->setValue(height_);
        areaMaskAngle->setValue(angle_);
    } else {
        center_x_ = areaMaskX->getValue();
        center_y_ = areaMaskY->getValue();
        width_ = areaMaskWidth->getValue();
        height_ = areaMaskHeight->getValue();
        angle_ = areaMaskAngle->getValue();
        updateGeometry();
    }
    enableListener();
}


inline const rtengine::ProcEvent &LabMasksPanel::areaMaskEvent() const
{
    return areaMask->getEnabled() ? EvAreaMask : EvAreaMaskVoid;
}


bool LabMasksPanel::button1Released()
{
    if (last_object_ != -1) {
        updateAreaMask(true);
        onAreaShapeSelectionChanged();
        populateShapeList(selected_, area_shape_index_);
        auto l = getListener();
        if (l) {
            l->panelChanged(areaMaskEvent(), M("GENERAL_CHANGED"));
        }
    }
    return AreaMask::button1Released();
}


void LabMasksPanel::switchOffEditMode()
{
    areaMaskToggle->set_active(false);
    AreaMask::switchOffEditMode();
}


void LabMasksPanel::onAreaMaskEnableToggled()
{
    auto l = getListener();
    if (l) {
        l->panelChanged(EvAreaMask, areaMask->getEnabled() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}


void LabMasksPanel::curveChanged(CurveEditor* ce)
{
    auto l = getListener();
    if (l) {
        if (ce == areaMaskContrast) {
            l->panelChanged(areaMaskEvent(), M("GENERAL_CHANGED"));
        } else if (ce == hueMask) {
            l->panelChanged(EvHMask, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == chromaticityMask) {
            l->panelChanged(EvCMask, M("HISTORY_CUSTOMCURVE"));
        } else if (ce == lightnessMask) {
            l->panelChanged(EvLMask, M("HISTORY_CUSTOMCURVE"));
        }
    }
}


void LabMasksPanel::adjusterChanged(Adjuster *a, double newval)
{
    auto l = getListener();

    if (a == maskBlur) {
        if (l) {
            l->panelChanged(EvMaskBlur, a->getTextValue());
        }
    } else if (std::find(areaMaskAdjusters.begin(), areaMaskAdjusters.end(), a) != areaMaskAdjusters.end()) {
        onAreaShapeSelectionChanged();
        populateShapeList(selected_, area_shape_index_);
        if (l) {
            l->panelChanged(areaMaskEvent(), M("GENERAL_CHANGED"));
        }
    } else if (a == deltaERange || a == deltaEDecay) {
        if (l) {
            l->panelChanged(deltaEMaskEvent(), M("GENERAL_CHANGED"));
        }
    } else if (a == contrastThreshold) {
        if (l) {
            l->panelChanged(EvContrastThresholdMask, a->getTextValue());
        }
    }
    maskShow(selected_, true);
}


void LabMasksPanel::adjusterChanged(ThresholdAdjuster *a, double newBottom, double newTop)
{
    if (a == deltaEL || a == deltaEC || a == deltaEH) {
        auto l = getListener();
        double b, L, C, H;
        deltaEL->getValue(b, L);
        deltaEC->getValue(b, C);
        deltaEH->getValue(b, H);
        static_cast<DeltaEArea *>(deltaEColor)->setColor(L, C, H);
        if (l) {
            l->panelChanged(deltaEMaskEvent(), M("GENERAL_CHANGED"));
        }
    }
}


void LabMasksPanel::adjusterAutoToggled(Adjuster *a, bool newval)
{
}


void LabMasksPanel::setMasks(const std::vector<rtengine::procparams::Mask> &masks, int show_mask_idx)
{
    disableListener();
    ConnectionBlocker b(selectionConn);
    
    masks_ = masks;
    selected_ = 0;
    if (show_mask_idx >= 0) {
        selected_ = show_mask_idx;
        showMask->set_active(true);
    } else {
        showMask->set_active(false);
    }
    static_cast<DrawnMaskPanel *>(drawnMask)->setTargetMask(nullptr);
    populateList();
    area_shape_index_ = 0;
    maskShow(selected_);
    enableListener();
}
        

void LabMasksPanel::getMasks(std::vector<rtengine::procparams::Mask> &masks, int &show_mask_idx)
{
    maskGet(selected_);
    masks = masks_;
    if (showMask->get_active()) {
        show_mask_idx = selected_;
    } else {
        show_mask_idx = -1;
    }
}


int LabMasksPanel::getSelected()
{
    return selected_;
}


void LabMasksPanel::colorForValue(double valX, double valY, enum ColorCaller::ElemType elemType, int callerId, ColorCaller *caller)
{
    float R = 0.f, G = 0.f, B = 0.f;
    double alpha = 0.f;

    const auto hue =
        [](float x) -> float
        {
            if (x >= 0.5f) {
                x -= 1.f;
            }
            x *= 2.f * rtengine::RT_PI_F;
            return rtengine::Color::huelab_to_huehsv2(x);
        };

    if (callerId == ID_HUE_MASK) {
        float x = valX - 1.f/6.f;
        if (x < 0.f) {
            x += 1.f;
        }
        x = rtengine::log2lin(x, 3.f);
        rtengine::Color::hsv2rgb01(x, 0.5f, 0.65f, R, G, B);
    } else if (callerId == ID_HUE_MASK+1) {
        rtengine::Color::hsv2rgb01(float(valY), float(valX), 0.8f, R, G, B);
    } else if (callerId == ID_HUE_MASK+2) {
        double dummy, w, h;
        deltaEH->getValue(dummy, h);
        rtengine::Color::hsv2rgb01(hue(float(h / 360.0)), 0.5f, float(valX), R, G, B);
        deltaEL->getValue(w, dummy);
        alpha = rtengine::LIM01(1.0 - w/100.0);
    } else if (callerId == ID_HUE_MASK+3) {
        double dummy, w, h;
        deltaEH->getValue(dummy, h);
        rtengine::Color::hsv2rgb01(hue(float(h / 360.0)), float(valX), 0.65f, R, G, B);
        deltaEC->getValue(w, dummy);
        alpha = rtengine::LIM01(1.0 - w/100.0);        
    } else if (callerId == ID_HUE_MASK+4) {
        double dummy, w;
        rtengine::Color::hsv2rgb01(hue(float(valX)), 0.5f, 0.65f, R, G, B);
        deltaEH->getValue(w, dummy);
        alpha = rtengine::LIM01(1.0 - w/100.0);        
    }

    caller->ccRed = double(R);
    caller->ccGreen = double(G);
    caller->ccBlue = double(B);
    caller->ccAlpha = alpha;
}


float LabMasksPanel::blendPipetteValues(CurveEditor *ce, float chan1, float chan2, float chan3)
{
    if (ce == chromaticityMask && chan1 > 0.f) {
        return rtengine::lin2log(chan1, 10.f);
    } else if (ce == hueMask && chan1 > 0.f) {
        float x = chan1 + 1.f/6.f;
        if (x > 1.f) {
            x -= 1.f;
        }
        return rtengine::lin2log(x, 3.f);
    }
    return CurveListener::blendPipetteValues(ce, chan1, chan2, chan3);
}


void LabMasksPanel::setEdited(bool yes)
{
    listEdited = yes;
    hueMask->setUnChanged(!yes);
    chromaticityMask->setUnChanged(!yes);
    lightnessMask->setUnChanged(!yes);
    maskBlur->setEditedState(yes ? Edited : UnEdited);
    maskInverted->set_inconsistent(!yes);
    showMask->set_inconsistent(!yes);
    areaMask->set_inconsistent(!yes);
    for (auto a : areaMaskAdjusters) {
        a->setEditedState(yes ? Edited : UnEdited);
    }
    areaMaskContrast->setUnChanged(!yes);
}


bool LabMasksPanel::getEdited()
{
    for (auto a : areaMaskAdjusters) {
        if (a->getEditedState() == Edited) {
            return true;
        }
    }
    return listEdited
        || !hueMask->isUnChanged()
        || !chromaticityMask->isUnChanged()
        || !lightnessMask->isUnChanged()
        || maskBlur->getEditedState() == Edited
        || !maskInverted->get_inconsistent()
        || !showMask->get_inconsistent()
        || !areaMask->get_inconsistent()
        || !areaMaskContrast->isUnChanged();
}


void LabMasksPanel::updateSelected()
{
    maskShow(selected_, true, false);
}


void LabMasksPanel::onAreaShapeSelectionChanged()
{
    if (selected_ < masks_.size() && area_shape_index_ < masks_[selected_].areaMask.shapes.size()) {
        disableListener();
        
        auto &s = masks_[selected_].areaMask.shapes[area_shape_index_];
        updateAreaMask(false);
        s.x = center_x_;
        s.y = center_y_;
        s.width = width_;
        s.height = height_;
        s.angle = angle_;
        s.roundness = areaMaskRoundness->getValue();
        s.mode = Shape::Mode(getAreaShapeMode());

        auto sel = areaMaskShapes->get_selected();
        unsigned int newidx = sel.empty() ? area_shape_index_ : sel[0];

        if (newidx != area_shape_index_) {
            areaShapeSelect(newidx, false);
        }
        enableListener();
    }
}


void LabMasksPanel::onAreaShapeResetPressed()
{
    if (selected_ < masks_.size()) {
        listEdited = true;
        disableListener();
        masks_[selected_].areaMask.shapes = {defaultAreaShape};
        area_shape_index_ = 0;
        center_x_ = defaultAreaShape.x;
        center_y_ = defaultAreaShape.y;
        width_ = defaultAreaShape.width;
        height_ = defaultAreaShape.height;
        angle_ = defaultAreaShape.angle;
        updateGeometry();
        updateAreaMask(true);
        populateShapeList(selected_, area_shape_index_);
        maskShow(selected_, true);        
        enableListener();
        auto l = getListener();
        if (l) {
            l->panelChanged(areaMaskEvent(), M("GENERAL_CHANGED"));
        }
    }
}


void LabMasksPanel::onAreaShapeAddPressed()
{
    if (selected_ < masks_.size()) {
        listEdited = true;
        auto &am = masks_[selected_].areaMask;
        auto s = am.shapes.empty() ? defaultAreaShape : am.shapes.back();
        am.shapes.push_back(s);
        populateShapeList(selected_, -1);
        areaShapeSelect(am.shapes.size()-1, true);
        maskShow(selected_, true);        
        auto l = getListener();
        if (l) {
            l->panelChanged(areaMaskEvent(), M("GENERAL_CHANGED"));
        }
    }
}


void LabMasksPanel::onAreaShapeRemovePressed()
{
    if (selected_ < masks_.size() && area_shape_index_ < masks_[selected_].areaMask.shapes.size() && masks_[selected_].areaMask.shapes.size() > 1) {
        listEdited = true;
        masks_[selected_].areaMask.shapes.erase(masks_[selected_].areaMask.shapes.begin() + area_shape_index_);
        populateShapeList(selected_, -1);
        areaShapeSelect(area_shape_index_ > 0 ? area_shape_index_ - 1 : 0, true);
        maskShow(selected_, true);
        auto l = getListener();
        if (l) {
            l->panelChanged(areaMaskEvent(), M("GENERAL_CHANGED"));
        }
    }
}


void LabMasksPanel::toggleAreaShapeMode(int i)
{
    for (int j = 0; j < 3; ++j) {
        ConnectionBlocker blocker(areaMaskModeConn[j]);
        areaMaskMode[j]->set_active(i == j);
    }
}


int LabMasksPanel::getAreaShapeMode()
{
    for (int j = 0; j < 3; ++j) {
        if (areaMaskMode[j]->get_active()) {
            return j;
        }
    }
    return 0;
}


void LabMasksPanel::onAreaShapeModeChanged(int i)
{
    if (!areaMaskMode[i]->get_active()) {
        ConnectionBlocker blocker(areaMaskModeConn[i]);
        areaMaskMode[i]->set_active(true);
        return;
    }
    toggleAreaShapeMode(i);
    onAreaShapeSelectionChanged();
    populateShapeList(selected_, area_shape_index_);
    auto l = getListener();
    if (l) {
        l->panelChanged(areaMaskEvent(), M("GENERAL_CHANGED"));
    }
    maskShow(selected_, true);
}


void LabMasksPanel::populateShapeList(int idx, int sel)
{
    ConnectionBlocker b(shapeSelectionConn);
    areaMaskShapes->clear_items();
    auto &r = masks_[idx];
    const auto rd =
        [](double v) -> double
        {
            return int(v * 10) / 10.0;
        };

    const auto m =
        [](Shape::Mode mode) -> const char *
        {
            switch (mode) {
            case Shape::SUBTRACT: return " (-)";
            case Shape::INTERSECT: return " (*)";
            default: return "";
            }
        };
    
    for (size_t i = 0; i < r.areaMask.shapes.size(); ++i) {
        auto &a = r.areaMask.shapes[i];
        auto j = areaMaskShapes->append(std::to_string(i+1));
        areaMaskShapes->set_text(j, 1, Glib::ustring::compose("%1 %2 %3 %4 %5 %6%7", rd(a.x), rd(a.y), rd(a.width), rd(a.height), rd(a.angle), rd(a.roundness), m(a.mode)));
    }
    if (sel >= 0) {
        Gtk::TreePath pth;
        pth.push_back(sel);
        areaMaskShapes->get_selection()->select(pth);
    }
}


void LabMasksPanel::onAreaMaskCopyPressed()
{
    if (selected_ < masks_.size()) {
        clipboard.setAreaMask(masks_[selected_].areaMask);
    }
}


void LabMasksPanel::onAreaMaskPastePressed()
{
    if (selected_ < masks_.size() && clipboard.hasAreaMask()) {
        listEdited = true;
        disableListener();
        masks_[selected_].areaMask = clipboard.getAreaMask();
        area_shape_index_ = 0;
        if (area_shape_index_ < masks_[selected_].areaMask.shapes.size()) {
            auto &a = masks_[selected_].areaMask;
            auto &s = a.shapes[0];
            center_x_ = s.x;
            center_y_ = s.y;
            width_ = s.width;
            height_ = s.height;
            angle_ = s.angle;
            updateGeometry();
            updateAreaMask(true);
            areaMaskFeather->setValue(a.feather);
            areaMaskBlur->setValue(a.blur);
            areaMaskContrast->setCurve(a.contrast);
            areaMaskRoundness->setValue(s.roundness);
            toggleAreaShapeMode(int(s.mode));
        }
        populateShapeList(selected_, area_shape_index_);
        maskShow(selected_, true);        
        enableListener();
        auto l = getListener();
        if (l) {
            l->panelChanged(areaMaskEvent(), M("GENERAL_CHANGED"));
        }
    }
}


void LabMasksPanel::areaShapeSelect(int sel, bool update_list)
{
    area_shape_index_ = sel;
    auto &ns = masks_[selected_].areaMask.shapes[sel];
    center_x_ = ns.x;
    center_y_ = ns.y;
    width_ = ns.width;
    height_ = ns.height;
    angle_ = ns.angle;
    updateGeometry();
    updateAreaMask(true);
    areaMaskRoundness->setValue(ns.roundness);
    toggleAreaShapeMode(int(ns.mode));
    if (areaMaskToggle->get_active()) {
        areaMaskToggle->set_active(false);
        areaMaskToggle->set_active(true);
    }

    if (update_list) {
        ConnectionBlocker b(shapeSelectionConn);
        Gtk::TreePath pth;
        pth.push_back(sel);
        areaMaskShapes->get_selection()->select(pth);
    }
}


void LabMasksPanel::setAreaDrawListener(AreaDrawListener *l)
{
    adl_ = l;
}


void LabMasksPanel::updateArea(AreaDrawUpdater::Phase phase, int x1, int y1, int x2, int y2)
{
    EditDataProvider *provider = getEditProvider();

    if (!provider) {
        return;
    }

    int imW = 0, imH = 0;
    provider->getImageSize(imW, imH);
    if (!imW || !imH) {
        return;
    }

    if (x1 > x2) {
        std::swap(x1, x2);
    }
    if (y1 > y2) {
        std::swap(y1, y2);
    }
    width_ = float(x2 - x1) * 100.f / imW;
    height_ = float(y2 - y1) * 100.f / imH;
    center_x_ = (x1 + (x2 - x1)/2 - imW/2) * 200.f / imW;
    center_y_ = (y1 + (y2 - y1)/2 - imH/2) * 200.f / imH;
    angle_ = 0;
    if (phase == AreaDrawUpdater::BEGIN) {
        subscribe();
    }
    updateGeometry();
    if (phase == AreaDrawUpdater::END) {
        updateAreaMask(true);
        onAreaShapeSelectionChanged();
        populateShapeList(selected_, area_shape_index_);
        ConnectionBlocker blocker(areaMaskDrawConn);
        areaMaskDraw->set_active(false);
        areaMaskToggle->set_active(true);

        auto l = getListener();
        if (l) {
            l->panelChanged(areaMaskEvent(), M("GENERAL_CHANGED"));
        }
    }
}


void LabMasksPanel::cancelUpdateArea()
{
    updateAreaMask(false);
    areaMaskToggle->set_active(true);
}


void LabMasksPanel::onAreaMaskDrawChanged()
{
    if (adl_) {
        if (areaMaskDraw->get_active()) {
            areaMaskToggle->set_active(false);
            adl_->startDrawingArea(this);
        } else {
            adl_->stopDrawingArea();
        }
    }
}


void LabMasksPanel::onAreaMaskDrawAddPressed()
{
    onAreaShapeAddPressed();
    areaMaskDraw->set_active(true);
}


void LabMasksPanel::on_map()
{
    Gtk::VBox::on_map();
    if (first_mask_exp_) {
        parametricMask->set_expanded(false);
        areaMask->set_expanded(false);
        deltaEMask->set_expanded(false);
        drawnMask->set_expanded(false);
        mask_exp_->set_expanded(false);
        first_mask_exp_ = false;
    }
}


void LabMasksPanel::onMaskFold(GdkEventButton *evt)
{
    if (mask_exp_->get_expanded()) {
        if (showMask->get_active()) {
            showMask->set_active(false);
        }
    }
}


void LabMasksPanel::onParametricMaskEnableToggled()
{
    auto l = getListener();
    if (l) {
        l->panelChanged(EvParametricMask, parametricMask->getEnabled() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}


void LabMasksPanel::onDeltaEMaskEnableToggled()
{
    auto l = getListener();
    if (l) {
        l->panelChanged(EvDeltaEMask, deltaEMask->getEnabled() ? M("GENERAL_ENABLED") : M("GENERAL_DISABLED"));
    }
}


inline const rtengine::ProcEvent &LabMasksPanel::deltaEMaskEvent() const
{
    return deltaEMask->getEnabled() ? EvDeltaEMask : EvDeltaEMaskVoid;
}


void LabMasksPanel::onDeltaEPickClicked()
{
    static_cast<DeltaEArea *>(deltaEColor)->activateSpot();
}


void LabMasksPanel::onDeltaESpotRequested(rtengine::Coord pos)
{
    if (deltaE_provider_) {
        deltaEMask->set_sensitive(false);
        float L, C, H;
        if (deltaE_provider_->getDeltaELCH(static_cast<DeltaEArea *>(deltaEColor)->getEditID(), pos, L, C, H)) {
            double b, t;
            deltaEL->getValue(b, t);
            deltaEL->setValue(b, L);
            deltaEC->getValue(b, t);
            deltaEC->setValue(b, C);
            deltaEH->getValue(b, t);
            deltaEH->setValue(b, H);
            static_cast<DeltaEArea *>(deltaEColor)->setColor(L, C, H);
            auto l = getListener();
            if (l) {
                l->panelChanged(EvDeltaEMask, M("GENERAL_CHANGED"));
            }
        }
        deltaEMask->set_sensitive(true);
    }
    static_cast<DeltaEArea *>(deltaEColor)->switchOffEditMode();
}


void LabMasksPanel::setDeltaEColorProvider(DeltaEColorProvider *p)
{
    deltaE_provider_ = p;
}


void LabMasksPanel::onListEnabledToggled(const Glib::ustring &path)
{
    auto it = list_model_->get_iter(path);
    if (it) {
        auto row = *it;
        row[list_model_columns_->enabled] = !row[list_model_columns_->enabled];
        auto idx = row[list_model_columns_->id] - 1;
        auto &r = masks_[idx];
        r.enabled = row[list_model_columns_->enabled];
        auto l = getListener();
        if (l) {
            l->panelChanged(EvMaskList, M("HISTORY_CHANGED"));
        }
    }
}


void LabMasksPanel::setListEnabled(Gtk::CellRenderer *renderer, const Gtk::TreeModel::iterator &it)
{
    auto row = *it;
    static_cast<Gtk::CellRendererToggle *>(renderer)->set_active(row[list_model_columns_->enabled]);
}


void LabMasksPanel::onDrawnMaskUpdated()
{
    auto l = getListener();
    if (l) {
        l->panelChanged(EvDrawnMask, M("GENERAL_CHANGED"));
    }
}


bool LabMasksPanel::onMaskNameFocusOut(GdkEventFocus *e)
{
    auto l = getListener();
    if (l) {
        l->panelChanged(EvMaskName, M("GENERAL_CHANGED"));
    }
    maskShow(selected_, true);
    return false;
}
