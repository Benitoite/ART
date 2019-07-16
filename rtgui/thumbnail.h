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
#ifndef _THUMBNAIL_
#define _THUMBNAIL_

#include <string>
#include <glibmm.h>
#include "cachemanager.h"
#include "options.h"
#include "../rtengine/rtengine.h"
#include "../rtengine/rtthumbnail.h"
#include "cacheimagedata.h"
#include "thumbnaillistener.h"
#include "threadutils.h"
#include "pparamschangelistener.h"

class CacheManager;
class Thumbnail: public PParamsSnapshotListener {
    MyMutex mutex;
    Glib::ustring fname;              // file name corresponding to the thumbnail
    CacheImageData cfs;                // cache entry corresponding to the thumbnail
    CacheManager *cachemgr;           // parent
    int ref;                // variable for reference counting
    int enqueueNumber;      // the number of instances in the batch queue corresponding to this thumbnail

    // if the thumbnail is in processed mode, this class holds its data:
    rtengine::Thumbnail *tpp;
    int tw, th;             // dimensions of timgdata (it stores tpp->width and tpp->height in processed mode for simplicity)
    float imgRatio;           // hack to avoid rounding error
//        double          scale;            // portion of the sizes of the processed thumbnail image and the full scale image

    //rtengine::procparams::ProcParams pparams;
    rtengine::procparams::ProcParamsCollection pparams;
    
    bool pparamsValid;
    bool needsReProcessing;
    bool imageLoading;

    // these are the data of the result image of the last getthumbnailimage  call (for caching purposes)
    unsigned char *lastImg;
    int lastW;
    int lastH;
    double lastScale;

    // exif & date/time strings
    Glib::ustring exifString;
    Glib::ustring dateTimeString;

    bool initial_;

    // rating info
    struct Rating {
        template <class T> struct Param {
            T value;
            bool edited;
            Param(T v): value(v), edited(false) {}
            Param &operator=(T v)
            {
                value = v;
                edited = true;
                return *this;
            }
            operator T() const { return value; }
        };
        Param<int> rank;
        Param<int> color;
        Param<bool> trash;

        explicit Rating(int r=0, int c=0, bool t=false):
            rank(r), color(c), trash(t) {}
        bool edited() const { return rank.edited || color.edited || trash.edited; }
    };
    Rating rating_;

    // vector of listeners
    std::vector<ThumbnailListener*> listeners;

    void _loadThumbnail(bool firstTrial = true);
    void _saveThumbnail();
    void _generateThumbnailImage();
    int infoFromImage(const Glib::ustring& fname);
    void loadThumbnail(bool firstTrial = true);
    void generateExifDateTimeStrings();

    Glib::ustring getCacheFileName(const Glib::ustring& subdir, const Glib::ustring& fext) const;

    void saveRating();
    void loadRating();
    void saveMetadata();

public:
    Thumbnail(CacheManager *cm, const Glib::ustring &fname, CacheImageData *cf);
    Thumbnail(CacheManager *cm, const Glib::ustring &fname, const std::string &md5);
    ~Thumbnail();

    bool hasProcParams();
    const rtengine::procparams::ProcParams &getProcParams();
    const rtengine::procparams::ProcParams &getProcParamsU();  // Unprotected version
    const std::vector<std::pair<Glib::ustring, rtengine::procparams::ProcParams>> &getProcParamsSnapshots();

    // Use this to create params on demand for update ; if flaggingMode=true, the procparams is created for a file being flagged (inTrash, rank, colorLabel)
    rtengine::procparams::ProcParams *createProcParamsForUpdate(bool returnParams, bool force, bool flaggingMode = false);

    void setProcParams(const rtengine::procparams::PartialProfile &pp, int whoChangedIt=-1, bool updateCacheNow=true, bool resetToDefault=false);
    void setProcParams(const rtengine::procparams::ProcParams &pp, int whoChangedIt=-1, bool updateCacheNow=true, bool resetToDefault=false);
    void clearProcParams(int whoClearedIt=-1);
    void loadProcParams(bool load_rating=true);

    void notifylisterners_procParamsChanged(int whoChangedIt);

    bool isQuick() { return cfs.thumbImgType == CacheImageData::QUICK_THUMBNAIL;}
    bool isPParamsValid() { return pparamsValid; }
    bool isRecentlySaved();
    void imageDeveloped();
    void imageEnqueued();
    void imageRemovedFromQueue();
    bool isEnqueued();
    bool isPixelShift();
    bool isHDR();

//        unsigned char*  getThumbnailImage (int &w, int &h, int fixwh=1); // fixwh = 0: fix w and calculate h, =1: fix h and calculate w
    rtengine::IImage8 *processThumbImage(const rtengine::procparams::ProcParams& pparams, int h, double& scale);
    rtengine::IImage8 *upgradeThumbImage(const rtengine::procparams::ProcParams& pparams, int h, double& scale);
    void getThumbnailSize(int &w, int &h, const rtengine::procparams::ProcParams *pparams = nullptr);
    void getFinalSize(const rtengine::procparams::ProcParams& pparams, int& w, int& h);
    void getOriginalSize(int& w, int& h);

    const Glib::ustring &getExifString();
    const Glib::ustring &getDateTimeString();
    void getCamWB(double& temp, double& green)
    {
        if (tpp) {
            tpp->getCamWB(temp, green);
        } else {
            temp = green = -1.0;
        }
    }
    void getAutoWB(double& temp, double& green, double equal, double tempBias);
    void getSpotWB(int x, int y, int rect, double& temp, double& green)
    {
        if (tpp) {
            tpp->getSpotWB (getProcParams(), x, y, rect, temp, green);
        } else {
            temp = green = -1.0;
        }
    }
    void applyAutoExp(rtengine::procparams::ProcParams& pparams)
    {
        if (tpp) {
            tpp->applyAutoExp (pparams);
        }
    }

    ThFileType getType();
    Glib::ustring getFileName() { return fname; }
    void setFileName(const Glib::ustring &fn);

    bool isSupported();

    const CacheImageData *getCacheImageData() { return &cfs; }
    std::string getMD5() { return cfs.md5; }

    int getRank()
    {
        return rating_.rank;
    }
    void setRank(int rank)
    {
        rating_.rank = rank;
    }

    int getColorLabel()
    {
        return rating_.color;
    }
    void setColorLabel(int colorlabel)
    {
        rating_.color = colorlabel;
    }

    bool getInTrash()
    {
        return rating_.trash;
    }
    void setInTrash(bool stage)
    {
        rating_.trash = stage;
    }

    void addThumbnailListener(ThumbnailListener *tnl);
    void removeThumbnailListener(ThumbnailListener *tnl);

    void increaseRef();
    void decreaseRef();

    void updateCache(bool updatePParams = true, bool updateCacheImageData = true);
    void saveThumbnail();

    bool openDefaultViewer(int destination);
    bool imageLoad(bool loading);

    std::shared_ptr<rtengine::FramesMetaData> getMetaData();

    static Glib::ustring getXmpSidecarPath(const Glib::ustring &path);

    void snapshotsChanged(const std::vector<std::pair<Glib::ustring, rtengine::procparams::ProcParams>> &snapshots) override;
};


#endif

