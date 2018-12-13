/*
 *  This file is part of RawTherapee.
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
 *
 *  (C) 2010 Emil Martinec <ejmartin@uchicago.edu>
 */

#include "dirpyrequalizer.h"
#include "eventmapper.h"

using namespace rtengine;
using namespace rtengine::procparams;


//-----------------------------------------------------------------------------
// DirPyrEqMasksContentProvider
//-----------------------------------------------------------------------------

class DirPyrEqMasksContentProvider: public LabMasksContentProvider {
public:
    DirPyrEqMasksContentProvider(DirPyrEqualizer *parent):
        parent_(parent)
    {
    }

    Gtk::Widget *getWidget() override
    {
        return parent_->box;
    }

    void getEvents(rtengine::ProcEvent &mask_list, rtengine::ProcEvent &h_mask, rtengine::ProcEvent &c_mask, rtengine::ProcEvent &l_mask, rtengine::ProcEvent &blur, rtengine::ProcEvent &show, rtengine::ProcEvent &area_mask) override
    {
        mask_list = parent_->EvList;
        h_mask = parent_->EvHueMask;
        c_mask = parent_->EvChromaticityMask;
        l_mask = parent_->EvLightnessMask;
        blur = parent_->EvMaskBlur;
        show = parent_->EvShowMask;
        area_mask = parent_->EvAreaMask;
    }

    ToolPanelListener *listener() override
    {
        if (parent_->getEnabled()) {
            return parent_->listener;
        }
        return nullptr;
    }

    void selectionChanging(int idx) override
    {
        parent_->levelsGet(idx);
    }

    void selectionChanged(int idx) override
    {
        parent_->levelsShow(idx);
    }

    bool addPressed() override
    {
        parent_->levelsData.push_back(DirPyrEqualizerParams::Levels());
        return true;
    }

    bool removePressed(int idx) override
    {
        parent_->levelsData.erase(parent_->levelsData.begin() + idx);
        return true;
    }
    
    bool copyPressed(int idx) override
    {
        parent_->levelsData.push_back(parent_->levelsData[idx]);
        return true;
    }
    
    bool moveUpPressed(int idx) override
    {
        auto r = parent_->levelsData[idx];
        parent_->levelsData.erase(parent_->levelsData.begin() + idx);
        --idx;
        parent_->levelsData.insert(parent_->levelsData.begin() + idx, r);
        return true;
    }
    
    bool moveDownPressed(int idx) override
    {
        auto r = parent_->levelsData[idx];
        parent_->levelsData.erase(parent_->levelsData.begin() + idx);
        ++idx;
        parent_->levelsData.insert(parent_->levelsData.begin() + idx, r);
        return true;
    }

    int getColumnCount() override
    {
        return 1;
    }
    
    Glib::ustring getColumnHeader(int col) override
    {
        return M("TP_CBDL_LIST_TITLE");
    }
    
    Glib::ustring getColumnContent(int col, int row) override
    {
        auto &r = parent_->levelsData[row];

        return Glib::ustring::compose(
            "%1 %2 %3\n%4 %5 %6 [%7]",
            Glib::ustring::format(std::fixed, std::setprecision(2), r.mult[0]),
            Glib::ustring::format(std::fixed, std::setprecision(2), r.mult[1]),
            Glib::ustring::format(std::fixed, std::setprecision(2), r.mult[2]),
            Glib::ustring::format(std::fixed, std::setprecision(2), r.mult[3]),
            Glib::ustring::format(std::fixed, std::setprecision(2), r.mult[4]),
            Glib::ustring::format(std::fixed, std::setprecision(2), r.mult[5]),
            Glib::ustring::format(std::fixed, std::setprecision(2), r.threshold));
    }

private:
    DirPyrEqualizer *parent_;
};


//-----------------------------------------------------------------------------
// DirPyrEqualizer
//-----------------------------------------------------------------------------

DirPyrEqualizer::DirPyrEqualizer(): FoldableToolPanel(this, "dirpyrequalizer", M("TP_DIRPYREQUALIZER_LABEL"), true, true)
{
    auto m = ProcEventMapper::getInstance();
    EvList = m->newEvent(SHARPENING, "HISTORY_MSG_CBDL_LIST");
    EvHueMask = m->newEvent(SHARPENING, "HISTORY_MSG_CBDL_HUEMASK");
    EvChromaticityMask = m->newEvent(SHARPENING, "HISTORY_MSG_CBDL_CHROMATICITYMASK");
    EvLightnessMask = m->newEvent(SHARPENING, "HISTORY_MSG_CBDL_LIGHTNESSMASK");
    EvMaskBlur = m->newEvent(SHARPENING, "HISTORY_MSG_CBDL_MASKBLUR");
    EvShowMask = m->newEvent(SHARPENING, "HISTORY_MSG_CBDL_SHOWMASK");
    EvAreaMask = m->newEvent(SHARPENING, "HISTORY_MSG_CBDL_AREAMASK");
    
    box = Gtk::manage(new Gtk::VBox());

    Gtk::HBox * buttonBox1 = Gtk::manage (new Gtk::HBox(true, 10));
    box->pack_start(*buttonBox1);

    Gtk::Button * lumacontrastMinusButton = Gtk::manage (new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMACONTRAST_MINUS")));
    buttonBox1->pack_start(*lumacontrastMinusButton);
    lumacontrastMinusButton->signal_pressed().connect( sigc::mem_fun(*this, &DirPyrEqualizer::lumacontrastMinusPressed));

    Gtk::Button * lumaneutralButton = Gtk::manage (new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMANEUTRAL")));
    buttonBox1->pack_start(*lumaneutralButton);
    lumaneutralButton->signal_pressed().connect( sigc::mem_fun(*this, &DirPyrEqualizer::lumaneutralPressed));

    Gtk::Button * lumacontrastPlusButton = Gtk::manage (new Gtk::Button(M("TP_DIRPYREQUALIZER_LUMACONTRAST_PLUS")));
    buttonBox1->pack_start(*lumacontrastPlusButton);
    lumacontrastPlusButton->signal_pressed().connect( sigc::mem_fun(*this, &DirPyrEqualizer::lumacontrastPlusPressed));

    buttonBox1->show_all_children();

    Gtk::HSeparator *separator2 = Gtk::manage (new  Gtk::HSeparator());
    box->pack_start(*separator2, Gtk::PACK_SHRINK, 2);

    for(int i = 0; i < 6; i++) {
        Glib::ustring ss;
        ss = Glib::ustring::format(i);

        if     (i == 0) {
            ss += Glib::ustring::compose(" (%1)", M("TP_DIRPYREQUALIZER_LUMAFINEST"));
        } else if(i == 5) {
            ss += Glib::ustring::compose(" (%1)", M("TP_DIRPYREQUALIZER_LUMACOARSEST"));
        }

        multiplier[i] = Gtk::manage ( new Adjuster (ss, 0, 4, 0.01, 1.0) );
        multiplier[i]->setAdjusterListener(this);
        box->pack_start(*multiplier[i]);
    }

    Gtk::HSeparator *separator3 = Gtk::manage (new  Gtk::HSeparator());
    box->pack_start(*separator3, Gtk::PACK_SHRINK, 2);

    threshold = Gtk::manage ( new Adjuster (M("TP_DIRPYREQUALIZER_THRESHOLD"), 0, 1, 0.01, 0.2) );
    threshold->setAdjusterListener(this);
    box->pack_start(*threshold);

    labMasksContentProvider.reset(new DirPyrEqMasksContentProvider(this));
    labMasks = Gtk::manage(new LabMasksPanel(labMasksContentProvider.get()));
    pack_start(*labMasks, Gtk::PACK_EXPAND_WIDGET, 4);   

    show_all_children ();
    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
}


DirPyrEqualizer::~DirPyrEqualizer ()
{

}

void DirPyrEqualizer::read (const ProcParams* pp, const ParamsEdited* pedited)
{
    disableListener();

    setEnabled(pp->dirpyrequalizer.enabled);

    levelsData = pp->dirpyrequalizer.levels;
    auto m = pp->dirpyrequalizer.labmasks;
    if (levelsData.empty()) {
        levelsData.emplace_back(rtengine::DirPyrEqualizerParams::Levels());
        m.emplace_back(rtengine::LabCorrectionMask());
    }
    labMasks->updateAreaMaskDefaults(pp);
    labMasks->setMasks(m, pp->dirpyrequalizer.showMask);
        
    
    if (pedited) {
        set_inconsistent(multiImage && !pedited->dirpyrequalizer.enabled);
        if (pedited->dirpyrequalizer.levels) {
            labMasks->setEdited(true);
            for (int i = 0; i < 6; i++) {
                multiplier[i]->setEditedState(Edited);
            }
            threshold->setEditedState(Edited);
        } else {
            labMasks->setEdited(false);
            for (int i = 0; i < 6; i++) {
                multiplier[i]->setEditedState(UnEdited);
            }
            threshold->setEditedState(UnEdited);
        }
    }

    enableListener();
}


void DirPyrEqualizer::write (ProcParams* pp, ParamsEdited* pedited)
{
    pp->dirpyrequalizer.enabled = getEnabled();

    levelsGet(labMasks->getSelected());
    pp->dirpyrequalizer.levels = levelsData;
    labMasks->getMasks(pp->dirpyrequalizer.labmasks, pp->dirpyrequalizer.showMask);
    assert(pp->dirpyrequalizer.levels.size() == pp->dirpyrequalizer.labmasks.size());

    labMasks->updateSelected();
        
    if (pedited) {
        pedited->dirpyrequalizer.enabled = !get_inconsistent();
        pedited->dirpyrequalizer.levels = labMasks->getEdited() || threshold->getEditedState();
        for (int i = 0; i < 6; i++) {
            if (multiplier[i]->getEditedState()) {
                pedited->dirpyrequalizer.levels = true;
            }
        }
    }
}


void DirPyrEqualizer::setDefaults (const ProcParams* defParams, const ParamsEdited* pedited)
{

    if (defParams->dirpyrequalizer.levels.size() == 1) {
        for (int i = 0; i < 6; i++) {
            multiplier[i]->setDefault(defParams->dirpyrequalizer.levels[0].mult[i]);
        }

        threshold->setDefault(defParams->dirpyrequalizer.levels[0].threshold);
    }

    if (pedited) {
        for (int i = 0; i < 6; i++) {
            multiplier[i]->setDefaultEditedState(pedited->dirpyrequalizer.levels ? Edited : UnEdited);
        }

        threshold->setDefaultEditedState(pedited->dirpyrequalizer.levels ? Edited : UnEdited);
    } else {
        for (int i = 0; i < 6; i++) {
            multiplier[i]->setDefaultEditedState(Irrelevant);
        }

        threshold->setDefaultEditedState(Irrelevant);
    }
}


void DirPyrEqualizer::setBatchMode (bool batchMode)
{
    // TODO
    ToolPanel::setBatchMode(batchMode);

    for (int i = 0; i < 6; i++) {
        multiplier[i]->showEditedCB();
    }

    threshold->showEditedCB();
    labMasks->setBatchMode();
}


void DirPyrEqualizer::adjusterChanged(Adjuster* a, double newval)
{
    if (listener && getEnabled()) {
        if (a == threshold) {
            listener->panelChanged (EvDirPyrEqualizerThreshold,
                                    Glib::ustring::compose("%1",
                                            Glib::ustring::format(std::fixed, std::setprecision(2), threshold->getValue()))
                                   );
        } else {
            listener->panelChanged (EvDirPyrEqualizer,
                                    Glib::ustring::compose("%1, %2, %3, %4, %5, %6",
                                            Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[0]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[1]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[2]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[3]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[4]->getValue()),
                                            Glib::ustring::format(std::fixed, std::setprecision(2), multiplier[5]->getValue()))
                                   );
        }
    }
}


void DirPyrEqualizer::adjusterAutoToggled(Adjuster* a, bool newval)
{
}


void DirPyrEqualizer::enabledChanged ()
{

    if (listener) {
        if (get_inconsistent()) {
            listener->panelChanged (EvDirPyrEqlEnabled, M("GENERAL_UNCHANGED"));
        } else if (getEnabled()) {
            listener->panelChanged (EvDirPyrEqlEnabled, M("GENERAL_ENABLED"));
        } else {
            listener->panelChanged (EvDirPyrEqlEnabled, M("GENERAL_DISABLED"));
        }
    }
}


void DirPyrEqualizer::lumaneutralPressed ()
{

    for (int i = 0; i < 6; i++) {
        multiplier[i]->setValue(1.0);
        adjusterChanged(multiplier[i], 1.0);
    }
}


void DirPyrEqualizer::lumacontrastPlusPressed ()
{

    for (int i = 0; i < 6; i++) {
        float inc = 0.05 * (6 - i);
        multiplier[i]->setValue(multiplier[i]->getValue() + inc);
        adjusterChanged(multiplier[i], multiplier[i]->getValue());
    }
}


void DirPyrEqualizer::lumacontrastMinusPressed ()
{

    for (int i = 0; i < 6; i++) {
        float inc = -0.05 * (6 - i);
        multiplier[i]->setValue(multiplier[i]->getValue() + inc);
        adjusterChanged(multiplier[i], multiplier[i]->getValue());
    }
}


void DirPyrEqualizer::setAdjusterBehavior (bool multiplieradd, bool thresholdadd, bool skinadd)
{
}

void DirPyrEqualizer::trimValues (rtengine::procparams::ProcParams* pp)
{
}


void DirPyrEqualizer::setEditProvider(EditDataProvider *provider)
{
    labMasks->setEditProvider(provider);
}


void DirPyrEqualizer::procParamsChanged(
    const rtengine::procparams::ProcParams* params,
    const rtengine::ProcEvent& ev,
    const Glib::ustring& descr,
    const ParamsEdited* paramsEdited)
{
    labMasks->updateAreaMaskDefaults(params);
}


void DirPyrEqualizer::updateGeometry(int fw, int fh)
{
    labMasks->updateGeometry(fw, fh);
}


void DirPyrEqualizer::levelsGet(int idx)
{
    if (idx < 0 || size_t(idx) >= levelsData.size()) {
        return;
    }
    
    auto &r = levelsData[idx];
    for (int i = 0; i < 6; ++i) {
        r.mult[i] = multiplier[i]->getValue();
    }
    r.threshold = threshold->getValue();
}


void DirPyrEqualizer::levelsShow(int idx)
{
    const bool disable = listener;
    if (disable) {
        disableListener();
    }

    auto &r = levelsData[idx];
    for (int i = 0; i < 6; ++i) {
        multiplier[i]->setValue(r.mult[i]);
    }
    threshold->setValue(r.threshold);
    
    if (disable) {
        enableListener();
    }
}
