/* -*- C++ -*-
 *  
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *  Copyright (c) 2010 Oliver Duis <www.oliverduis.de>
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
#include <memory>
#include "imageareapanel.h"
#include "toolpanelcoord.h"
#include "profilepanel.h"
#include "../rtengine/rtengine.h"
#include "history.h"
#include "histogrampanel.h"
#include "thumbnail.h"
#include "saveasdlg.h"
#include "batchqueueentry.h"
#include "thumbnaillistener.h"
#include "navigator.h"
#include "progressconnector.h"
#include "filepanel.h"

class EditorPanel;
class MyProgressBar;

struct EditorPanelIdleHelper {
    EditorPanel* epanel;
    bool destroyed;
    int pending;
};

class RTWindow;
using rtengine::array2D;

class EditorPanel final :
    public Gtk::VBox,
    public PParamsChangeListener,
    public rtengine::ProgressListener,
    public ThumbnailListener,
    public HistoryBeforeAfterListener,
    public rtengine::HistogramListener,
    public HistogramPanelListener,
    public rtengine::SizeListener
{
public:
    explicit EditorPanel (FilePanel* filePanel = nullptr);
    ~EditorPanel () override;

    void open (Thumbnail* tmb, rtengine::InitialImage* isrc);
    void setAspect ();
    void on_realize () override;
    void leftPaneButtonReleased (GdkEventButton *event);
    void rightPaneButtonReleased (GdkEventButton *event);

    void setParent(RTWindow* p);

    void setParentWindow (Gtk::Window* p)
    {
        parentWindow = p;
    }

    void writeOptions();
    void writeToolExpandedStatus (std::vector<int> &tpOpen);

    void showTopPanel (bool show);
    bool isRealized()
    {
        return realized;
    }
    // ProgressListener interface
    void setProgress(double p) override;
    void setProgressStr(const Glib::ustring& str) override;
    void setProgressState(bool inProcessing) override;
    void error(const Glib::ustring& descr) override;

    void error(const Glib::ustring& title, const Glib::ustring& descr);
    void displayError(const Glib::ustring& title, const Glib::ustring& descr);  // this is called by error in the gtk thread
    void refreshProcessingState (bool inProcessing); // this is called by setProcessingState in the gtk thread

    // PParamsChangeListener interface
    void procParamsChanged(
        const rtengine::procparams::ProcParams* params,
        const rtengine::ProcEvent& ev,
        const Glib::ustring& descr,
        const ParamsEdited* paramsEdited = nullptr
    ) override;
    void clearParamChanges() override;

    // thumbnaillistener interface
    void procParamsChanged (Thumbnail* thm, int whoChangedIt) override;

    // HistoryBeforeAfterListener
    void historyBeforeAfterChanged(const rtengine::procparams::ProcParams& params) override;

    // HistogramListener
    void histogramChanged(
        const LUTu& histRed,
        const LUTu& histGreen,
        const LUTu& histBlue,
        const LUTu& histLuma,
        const LUTu& histToneCurve,
        const LUTu& histLCurve,
        const LUTu& histCCurve,
        const LUTu& histLCAM,
        const LUTu& histCCAM,
        const LUTu& histRedRaw,
        const LUTu& histGreenRaw,
        const LUTu& histBlueRaw,
        const LUTu& histChroma,
        const LUTu& histLRETI,
        int vectorscopeScale,
        const array2D<int>& vectorscopeHC,
        const array2D<int>& vectorscopeHS,
        int waveformScale,
        const array2D<int>& waveformRed,
        const array2D<int>& waveformGreen,
        const array2D<int>& waveformBlue,
        const array2D<int>& waveformLuma
    ) override;
    void setObservable(rtengine::HistogramObservable* observable) override;
    bool updateHistogram(void) const override;
    bool updateHistogramRaw(void) const override;
    bool updateVectorscopeHC(void) const override;
    bool updateVectorscopeHS(void) const override;
    bool updateWaveform(void) const override;

    // HistogramPanelListener
    void scopeTypeChanged(Options::ScopeType new_type) override;

    // SizeListener
    void sizeChanged(int w, int h, int ow, int oh) override;
    
    // event handlers
    void info_toggled ();
    void hideHistoryActivated ();
    void tbRightPanel_1_toggled ();
    void tbTopPanel_1_toggled ();
    void beforeAfterToggled ();
    void tbBeforeLock_toggled();
    void saveAsPressed(GdkEventButton *event);
    void queueImgPressed(GdkEventButton *event);
    void sendToGimpPressed ();
    void openNextEditorImage ();
    void openPreviousEditorImage ();
    void syncFileBrowser ();

    void tbTopPanel_1_visible (bool visible);
    bool CheckSidePanelsVisibility();
    void tbShowHideSidePanels_managestate();
    void toggleSidePanels();
    void toggleSidePanelsZoomFit();

    void saveProfile ();
    Glib::ustring getShortName ();
    Glib::ustring getFileName ();
    bool handleShortcutKey (GdkEventKey* event);
    bool keyPressedBefore(GdkEventKey *event);
    bool keyReleased(GdkEventKey *event);
    bool scrollPressed(GdkEventScroll *event);

    bool getIsProcessing() const
    {
        return isProcessing || !can_open_now();
    }
    void setIsProcessing() { isProcessing = true; }
    
    void updateProfiles (const Glib::ustring &printerProfile, rtengine::RenderingIntent printerIntent, bool printerBPC);
    void updateTPVScrollbar (bool hide);
    void updateHistogramPosition (int oldPosition, int newPosition);

    void defaultMonitorProfileChanged (const Glib::ustring &profile_name, bool auto_monitor_profile);

    bool saveImmediately (const Glib::ustring &filename, const SaveFormat &sf);

    Gtk::Paned* catalogPane;

    void cleanup();
    void close ();
    
private:

    BatchQueueEntry *createBatchQueueEntry(bool fast_export, bool use_batch_queue_profile, const rtengine::procparams::PartialProfile *export_profile);
    bool idle_imageSaved(ProgressConnector<int> *pc, rtengine::IImagefloat* img, Glib::ustring fname, SaveFormat sf, rtengine::procparams::ProcParams &pparams);
    bool idle_saveImage(ProgressConnector<rtengine::IImagefloat*> *pc, Glib::ustring fname, SaveFormat sf, rtengine::procparams::ProcParams &pparams);
    bool idle_sendToGimp( ProgressConnector<rtengine::IImagefloat*> *pc, Glib::ustring fname);
    bool idle_sentToGimp(ProgressConnector<int> *pc, rtengine::IImagefloat* img, Glib::ustring filename);
    void histogramProfile_toggled();

    void do_save_image(bool fast_export);
    void do_queue_image(bool fast_export);
    bool autosave();

    bool can_open_now() const;
    
    Glib::ustring lastSaveAsFileName;
    bool realized;

    MyProgressBar  *progressLabel;
    Gtk::ToggleButton* info;
    Gtk::ToggleButton* hidehp;
    Gtk::ToggleButton* tbShowHideSidePanels;
    Gtk::ToggleButton* tbTopPanel_1;
    Gtk::ToggleButton* tbRightPanel_1;
    Gtk::ToggleButton* tbBeforeLock;
    //bool bAllSidePanelsVisible;
    Gtk::ToggleButton* beforeAfter;
    Gtk::Paned* hpanedl;
    Gtk::Paned* hpanedr;
    Gtk::Image *iHistoryShow, *iHistoryHide;
    Gtk::Image *iTopPanel_1_Show, *iTopPanel_1_Hide;
    Gtk::Image *iRightPanel_1_Show, *iRightPanel_1_Hide;
    Gtk::Image *iShowHideSidePanels;
    Gtk::Image *iShowHideSidePanels_exit;
    Gtk::Image *iBeforeLockON, *iBeforeLockOFF;
    Gtk::Paned *leftbox;
    Gtk::Box *leftsubbox;
    Gtk::Paned *vboxright;
    Gtk::Box *vsubboxright;

    Gtk::Button* queueimg;
    Gtk::Button* saveimgas;
    Gtk::Button* sendtogimp;
    Gtk::Button* navSync;
    Gtk::Button* navNext;
    Gtk::Button* navPrev;

    class ColorManagementToolbar;
    std::unique_ptr<ColorManagementToolbar> colorMgmtToolBar;

    ImageAreaPanel* iareapanel;
    PreviewHandler* previewHandler;
    PreviewHandler* beforePreviewHandler;   // for the before-after view
    Navigator* navigator;
    ImageAreaPanel* beforeIarea;    // for the before-after view
    Gtk::Box* beforeBox;
    Gtk::Box* afterBox;
    Gtk::Label* beforeLabel;
    Gtk::Label* afterLabel;
    Gtk::Box* beforeAfterBox;
    Gtk::Box* beforeHeaderBox;
    Gtk::Box* afterHeaderBox;
    Gtk::ToggleButton* toggleHistogramProfile;

    Gtk::Frame* ppframe;
    ProfilePanel* profilep;
    History* history;
    HistogramPanel* histogramPanel;
    ToolPanelCoordinator* tpc;
    RTWindow* parent;
    Gtk::Window* parentWindow;
    //SaveAsDialog* saveAsDialog;
    FilePanel* fPanel;

    bool firstProcessingDone;

    Thumbnail* openThm;  // may get invalid on external delete event
    Glib::ustring fname;  // must be saved separately

    int selectedFrame;

    rtengine::InitialImage* isrc;
    std::shared_ptr<rtengine::StagedImageProcessor> ipc;
    std::shared_ptr<rtengine::StagedImageProcessor> beforeIpc;    // for the before-after view

    EditorPanelIdleHelper* epih;

    int err;

    time_t processingStartedTime;

    sigc::connection ShowHideSidePanelsconn;

    bool isProcessing;

    IdleRegister idle_register;

    rtengine::HistogramObservable* histogram_observable;
    Options::ScopeType histogram_scope_type;

    sigc::connection autosave_conn_;

    std::unique_ptr<ToolShortcutManager> shortcut_mgr_;
};

