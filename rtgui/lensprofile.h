/*
*  This file is part of RawTherapee.
*
*  Copyright (c) 2012 Oliver Duis <oduis@oliverduis.de>
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
#ifndef _LENSPROFILE_H_
#define _LENSPROFILE_H_

#include <gtkmm.h>
#include "toolpanel.h"
#include "guiutils.h"
#include "lensgeom.h"

class LensProfilePanel : public ToolParamBlock, public FoldableToolPanel
{

protected:

    MyFileChooserButton *fcbLCPFile;
    Gtk::CheckButton *ckbUseDist, *ckbUseVign, *ckbUseCA;
    Gtk::HBox *hbLCPFile;
    Gtk::Label *lLCPFileHead;
    bool lcModeChanged, lcpFileChanged, useDistChanged, useVignChanged, useCAChanged;
    sigc::connection conLCPFile, conUseDist, conUseVign, conUseCA;
    void updateDisabled(bool enable);
    bool allowFocusDep;
    bool isRaw;
    const rtengine::FramesMetaData* metadata;

    Gtk::RadioButton::Group corrGroup;
    Gtk::RadioButton *corrOff;
    Gtk::RadioButton *corrLensfunAuto;
    Gtk::RadioButton *corrLensfunManual;
    Gtk::RadioButton *corrLcpFile;
    Gtk::RadioButton *corrUnchanged;
    MyComboBox *lensfunCameras;
    MyComboBox *lensfunLenses;
    Gtk::Image *warning;

    class LFDbHelper {
    public:
        class LFModelCam: public Gtk::TreeModel::ColumnRecord {
        public:
            LFModelCam() { add(make); add(model); }
            Gtk::TreeModelColumn<Glib::ustring> make;
            Gtk::TreeModelColumn<Glib::ustring> model;
        };

        class LFModelLens: public Gtk::TreeModel::ColumnRecord {
        public:
            LFModelLens() { add(lens); add(prettylens); }
            Gtk::TreeModelColumn<Glib::ustring> lens;
            Gtk::TreeModelColumn<Glib::ustring> prettylens;
        };

        LFModelCam lensfunModelCam;
        LFModelLens lensfunModelLens;
    
        Glib::RefPtr<Gtk::TreeStore> lensfunCameraModel;
        Glib::RefPtr<Gtk::TreeStore> lensfunLensModel;

        LFDbHelper();
        void fillLensfunCameras();
        void fillLensfunLenses();
    };
    static LFDbHelper *lf;

    bool useLensfunChanged;
    bool lensfunAutoChanged;
    bool lensfunCameraChanged;
    bool lensfunLensChanged;

    bool setLensfunCamera(const Glib::ustring &make, const Glib::ustring &model);
    bool setLensfunLens(const Glib::ustring &lens);
    bool checkLensfunCanCorrect(bool automatch);
    void updateLensfunWarning();
    
public:

    LensProfilePanel ();

    void read           (const rtengine::procparams::ProcParams* pp, const ParamsEdited* pedited = nullptr);
    void write          (rtengine::procparams::ProcParams* pp, ParamsEdited* pedited = nullptr);
    void setRawMeta     (bool raw, const rtengine::FramesMetaData* pMeta);

    void onLCPFileChanged ();
    void onUseDistChanged();
    void onUseVignChanged();
    void onUseCAChanged();

    void setBatchMode(bool yes);

    void onLensfunCameraChanged();
    void onLensfunLensChanged();
    void onCorrModeChanged(const Gtk::RadioButton *rbChanged);
};

#endif
