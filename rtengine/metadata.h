/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright (c) 2019 Alberto Griggio <alberto.griggio@gmail.com>
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

#include <glibmm.h>
#include <exiv2/exiv2.hpp>
#include <memory>
#include "procparams.h"

namespace rtengine {

class Exiv2Metadata {
public:
    Exiv2Metadata();
    explicit Exiv2Metadata(const Glib::ustring &path);
    Exiv2Metadata(const Glib::ustring &path, bool merge_xmp_sidecar);

    void load() const;
    
    Exiv2::ExifData &exifData() { return image_.get() ? image_->exifData() : exif_data_; }
    const Exiv2::ExifData &exifData() const { return const_cast<Exiv2Metadata *>(this)->exifData(); }
    
    Exiv2::IptcData &iptcData() { return image_.get() ? image_->iptcData() : iptc_data_; }
    const Exiv2::IptcData &iptcData() const { return const_cast<Exiv2Metadata *>(this)->iptcData(); }
    
    Exiv2::XmpData &xmpData() { return image_.get() ? image_->xmpData() : xmp_data_; }
    const Exiv2::XmpData &xmpData() const { return const_cast<Exiv2Metadata *>(this)->xmpData(); }

    const Glib::ustring &filename() const { return src_; }
    const rtengine::procparams::ExifPairs &exif() const { return exif_; }
    const rtengine::procparams::IPTCPairs &iptc() const { return iptc_; }
    void setExif(const rtengine::procparams::ExifPairs &exif) { exif_ = exif; }
    void setIptc(const rtengine::procparams::IPTCPairs &iptc) { iptc_ = iptc; }
    
    void saveToImage(const Glib::ustring &path) const;
    void saveToXmp(const Glib::ustring &path) const;

    static Glib::ustring xmpSidecarPath(const Glib::ustring &path);
    static Exiv2::XmpData getXmpSidecar(const Glib::ustring &path);

    static void init(const Glib::ustring &base_dir);
    static void cleanup();
   
private:
    void do_merge_xmp(Exiv2::Image *dst) const;
    void import_exif_pairs(Exiv2::ExifData &out) const;
    void import_iptc_pairs(Exiv2::IptcData &out) const;
    void remove_unwanted(Exiv2::Image *dst) const;
    
    Glib::ustring src_;
    bool merge_xmp_;
    mutable std::shared_ptr<Exiv2::Image> image_;
    rtengine::procparams::ExifPairs exif_;
    rtengine::procparams::IPTCPairs iptc_;
    Exiv2::ExifData exif_data_;
    Exiv2::IptcData iptc_data_;
    Exiv2::XmpData xmp_data_;
};

// Glib::ustring get_xmp_sidecar_path(const Glib::ustring &path);
// Exiv2::Image::AutoPtr open_exiv2(const Glib::ustring &fname,
//                                  bool merge_xmp_sidecar);
// Exiv2::XmpData read_exiv2_xmp(const Glib::ustring &fname);

} // namespace rtengine
