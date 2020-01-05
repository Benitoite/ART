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

#include <stdio.h>
#include <glib/gstdio.h>
#include <iostream>
#include <unistd.h>
#include <giomm.h>

#ifdef WIN32
#  include <windows.h>
#  include <io.h>
#  include <fcntl.h>
#  include <set>
#endif

#include "metadata.h"
#include "settings.h"
#include "../rtgui/version.h"
#include "../rtgui/pathutils.h"


namespace rtengine {

extern const Settings *settings;

std::unique_ptr<Exiv2Metadata::ImageCache> Exiv2Metadata::cache_(nullptr);

namespace {

constexpr size_t IMAGE_CACHE_SIZE = 200;

#ifdef WIN32
std::wstring to_wstr(const Glib::ustring &s)
{
    auto *ws = g_utf8_to_utf16(s.c_str(), -1, NULL, NULL, NULL);
    std::wstring ret(reinterpret_cast<wchar_t *>(ws));
    g_free(ws);
    return ret;
}
#endif // WIN32


Exiv2::Image::AutoPtr open_exiv2(const Glib::ustring &fname)
{
#if defined WIN32 && defined EXV_UNICODE_PATH
    std::wstring wfname = to_wstr(fname);
    auto image = Exiv2::ImageFactory::open(wfname);
#else
    auto image = Exiv2::ImageFactory::open(Glib::filename_from_utf8(fname));
#endif
    image->readMetadata();
    if (!image->good()) {
        throw Exiv2::Error(Exiv2::kerErrorMessage, "exiv2: invalid image");
    }
    return image;
}


Glib::ustring exiftool_base_dir;
#ifdef WIN32
  const Glib::ustring exiftool_default = "exiftool.exe";
#else
  const Glib::ustring exiftool_default = "exiftool";
#endif


#ifdef WIN32

// Glib::spawn_sync opens a console window for command-line apps, I wasn't
// able to find out how not to do that (see also:
// http://gtk.10911.n7.nabble.com/g-spawn-on-windows-td84743.html).
// Therefore, we roll our own
bool exec_subprocess(const std::vector<Glib::ustring> &argv, std::string &out, std::string &err)
{
    const auto add_quoted =
        [](std::wostream &out, const std::wstring &ws) -> void
        {
            out << '"';
            for (size_t j = 0; j < ws.size(); ) {
                int backslashes = 0;
                while (j < ws.size() && ws[j] == '\\') {
                    ++backslashes;
                    ++j;
                }
                if (j == ws.size()) {
                    backslashes = backslashes * 2;
                } else if (ws[j] == '"') {
                    backslashes = backslashes * 2 + 1;
                }
                for (int i = 0; i < backslashes; ++i) {
                    out << '\\';
                }
                if (j < ws.size()) {
                    out << ws[j];
                    ++j;
                } else {
                    break;
                }
            }
            out << '"';
        };

    struct HandleCloser {
        ~HandleCloser()
        {
            for (auto h : toclose) {
                CloseHandle(h);
            }
        }
        std::set<HANDLE> toclose;
    };

    HANDLE fds_from[2];
    HANDLE fds_from_e[2];
    SECURITY_ATTRIBUTES sa;
    HandleCloser hc;

    sa.nLength = sizeof(SECURITY_ATTRIBUTES);
    sa.lpSecurityDescriptor = NULL;
    sa.bInheritHandle = TRUE;

    const auto mkpipe =
        [&](HANDLE *fd) -> bool
        {
            if (!CreatePipe(&(fd[0]), &(fd[1]), &sa, 0)) {
                return false;
            }
            hc.toclose.insert(fd[0]);
            hc.toclose.insert(fd[1]);
            if (!SetHandleInformation(fd[0], HANDLE_FLAG_INHERIT, 0)) {
                return false;
            }
            return true;
        };

    if (!mkpipe(fds_from) || !mkpipe(fds_from_e)) {
        return false;
    }

    PROCESS_INFORMATION pi;
    STARTUPINFOW si;

    ZeroMemory(&si, sizeof(STARTUPINFOW));
    si.cb = sizeof(STARTUPINFOW);
    si.dwFlags = STARTF_USESTDHANDLES;
    si.wShowWindow = SW_HIDE;
    si.hStdInput = GetStdHandle(STD_INPUT_HANDLE);
    si.hStdOutput = fds_from[1];
    si.hStdError = fds_from_e[1];

    std::wstring pth = to_wstr(argv[0]);
    if (Glib::path_get_basename(argv[0]) == argv[0]) {
        wchar_t pathbuf[MAX_PATH+1];
        int n = SearchPathW(nullptr, pth.c_str(), nullptr, MAX_PATH+1, pathbuf, nullptr);
        if (n > 0) {
            pth = pathbuf;
        }
    }

    wchar_t *cmdline = nullptr;
    {
        std::wostringstream cmdlinebuf;
        add_quoted(cmdlinebuf, pth);
        for (size_t i = 1; i < argv.size(); ++i) {
            cmdlinebuf << ' ';
            add_quoted(cmdlinebuf, to_wstr(argv[i]));
        }
        std::wstring s = cmdlinebuf.str();
        cmdline = new wchar_t[s.size()+1];
        memcpy(cmdline, s.c_str(), s.size() * sizeof(wchar_t));
        cmdline[s.size()] = 0;
    }

    if (!CreateProcessW(pth.c_str(), cmdline, nullptr, nullptr, TRUE,
                        CREATE_NO_WINDOW,
                        (LPVOID)nullptr, nullptr, &si, &pi)) {
        delete[] cmdline;
        return false;
    } else {
        hc.toclose.insert(pi.hProcess);
        hc.toclose.insert(pi.hThread);
    }
    delete[] cmdline;

    const auto read_pipe =
        [&](HANDLE *fd) -> std::string
        {
            constexpr size_t bufsize = 4096;
            unsigned char buf[bufsize];
            std::ostringstream sbuf;
            DWORD n;

            hc.toclose.erase(fd[1]);
            CloseHandle(fd[1]);

            while (ReadFile(fd[0], buf, bufsize, &n, nullptr)) {
                buf[n] = 0;
                sbuf << buf;
                if (n < bufsize) {
                    break;
                }
            }
            return sbuf.str();
        };

    out = read_pipe(fds_from);
    err = read_pipe(fds_from_e);

    unsigned int status = 255;
    const DWORD wait_timeout_ms = 1000; // INFINITE
    if (WaitForSingleObject(pi.hProcess, wait_timeout_ms) != WAIT_OBJECT_0) {
        TerminateProcess(pi.hProcess, status);
    }
    if (!GetExitCodeProcess(pi.hProcess, (LPDWORD)&status)) {
        status = 255;
    }

    return status == 0;
}

#else // WIN32

bool exec_subprocess(const std::vector<Glib::ustring> &argv, std::string &out, std::string &err)
{
    std::vector<std::string> args;
    args.reserve(argv.size());
    for (auto &s : argv) {
        args.push_back(Glib::filename_from_utf8(s));
    }
    int exit_status = -1;
    Glib::spawn_sync("", args, Glib::SPAWN_DEFAULT|Glib::SPAWN_SEARCH_PATH, Glib::SlotSpawnChildSetup(), &out, &err, &exit_status);
    return WIFEXITED(exit_status) && WEXITSTATUS(exit_status) == 0;
}

#endif // WIN32


Exiv2::Image::AutoPtr exiftool_import(const Glib::ustring &fname, const std::exception &exc)
{
    Glib::ustring exiftool = settings->exiftool_path;
    if (exiftool == exiftool_default) {
        Glib::ustring e = Glib::build_filename(exiftool_base_dir, exiftool);
        if (Glib::file_test(e, Glib::FILE_TEST_EXISTS)) {
            exiftool = e;
        }
    }
        
    std::string templ = Glib::build_filename(Glib::get_tmp_dir(), Glib::ustring::compose("ART-exiftool-%1-XXXXXX", Glib::path_get_basename(fname)));
    int fd = Glib::mkstemp(templ);
    if (fd < 0) {
        throw exc;
    }
    Glib::ustring outname = Glib::filename_to_utf8(templ) + ".xmp";
    // int exit_status = -1;
    std::vector<Glib::ustring> argv = {
        exiftool,
        "-TagsFromFile",
        fname,
        "-xmp:all<all",       
        outname
    };
    if (settings->verbose) {
        std::cout << "importing metadata for " << fname << " with exiftool"
                  << std::endl;
    }
    std::string out, err;
    bool ok = exec_subprocess(argv, out, err);
    close(fd);
    g_remove(templ.c_str());
    if (settings->verbose) {
        if (!out.empty()) {
            std::cout << "  exiftool stdout: " << out << std::flush;
        }
        if (!err.empty()) {
            std::cout << "  exiftool stderr: " << err << std::flush;
        }
    }
    if (!ok) {
        if (Glib::file_test(outname, Glib::FILE_TEST_EXISTS)) {
            g_remove(outname.c_str());
        }
        throw exc;
    }
    try {
        auto image = Exiv2::ImageFactory::open(outname);
        image->readMetadata();
        auto &exif = image->exifData();
        auto &xmp = image->xmpData();
        const auto set_from =
            [&](const char *src, const char *dst) -> void
            {
                auto dk = Exiv2::ExifKey(dst);
                auto pos = exif.findKey(dk);
                if (pos == exif.end() || !pos->size()) {
                    auto sk = Exiv2::XmpKey(src);
                    auto it = xmp.findKey(sk);
                    if (it != xmp.end() && it->size()) {
                        exif[dst] = it->toString();
                    }
                }
            };
        set_from("Xmp.exifEX.LensModel", "Exif.Photo.LensModel");
        xmp.clear();
        g_remove(outname.c_str());
        return image;
    } catch (std::exception &) {
        if (Glib::file_test(outname, Glib::FILE_TEST_EXISTS)) {
            g_remove(outname.c_str());
        }
        throw exc;
    }
    return Exiv2::Image::AutoPtr();
}

} // namespace


Exiv2Metadata::Exiv2Metadata():
    src_(""),
    merge_xmp_(false),
    image_(nullptr)
{
}


Exiv2Metadata::Exiv2Metadata(const Glib::ustring &path):
    src_(path),
    merge_xmp_(settings->metadata_xmp_sync != Settings::MetadataXmpSync::NONE),
    image_(nullptr)
{
}


Exiv2Metadata::Exiv2Metadata(const Glib::ustring &path, bool merge_xmp_sidecar):
    src_(path),
    merge_xmp_(merge_xmp_sidecar),
    image_(nullptr)
{
}


void Exiv2Metadata::load() const
{
    if (!src_.empty() && !image_.get() && Glib::file_test(src_.c_str(), Glib::FILE_TEST_EXISTS)) {
        CacheVal val;
        auto finfo = Gio::File::create_for_path(src_)->query_info(G_FILE_ATTRIBUTE_TIME_MODIFIED);
        if (cache_ && cache_->get(src_, val) && val.second >= finfo->modification_time()) {
            if (settings->verbose) {
                std::cout << "Metadata for " << src_ << " found in cache" << std::endl;
            }
            image_ = val.first;
        } else {
            try {
                auto img = open_exiv2(src_);
                image_.reset(img.release());
            } catch (std::exception &exc) {
                auto img = exiftool_import(src_, exc);
                image_.reset(img.release());
            }
            if (cache_) {
                cache_->set(src_, CacheVal(image_, finfo->modification_time()));
            }
        }

        if (merge_xmp_) {
            do_merge_xmp(image_.get());
        }
    }
}


void Exiv2Metadata::do_merge_xmp(Exiv2::Image *dst) const
{
    try { 
        auto xmp = getXmpSidecar(src_);
        Exiv2::ExifData exif;
        Exiv2::IptcData iptc;
        Exiv2::moveXmpToIptc(xmp, iptc);
        Exiv2::moveXmpToExif(xmp, exif);

        for (auto &datum : exif) {
            dst->exifData()[datum.key()] = datum;
        }
        for (auto &datum : iptc) {
            dst->iptcData()[datum.key()] = datum;
        }
        for (auto &datum : xmp) {
            dst->xmpData()[datum.key()] = datum;
        }
    } catch (std::exception &exc) {
        if (settings->verbose) {
            std::cerr << "Error loading metadata from XMP sidecar: "
                      << exc.what() << std::endl;
        }
    }
}


void Exiv2Metadata::saveToImage(const Glib::ustring &path) const
{
    auto dst = open_exiv2(path);
    if (image_.get()) {
        dst->setMetadata(*image_);
        if (merge_xmp_) {
            do_merge_xmp(dst.get());
        }
        remove_unwanted(dst.get());
    } else {
        dst->setExifData(exif_data_);
        dst->setIptcData(iptc_data_);
        dst->setXmpData(xmp_data_);
    }

    dst->exifData()["Exif.Image.Software"] = RTNAME " " RTVERSION;
    import_exif_pairs(dst->exifData());
    import_iptc_pairs(dst->iptcData());
    dst->writeMetadata();    
}


void Exiv2Metadata::remove_unwanted(Exiv2::Image *dst) const
{
    static const std::vector<std::string> keys = {
        "Exif.Image.Orientation",
        "Exif.Image2.JPEGInterchangeFormat",
        "Exif.Image2.JPEGInterchangeFormatLength"
    };
    for (auto &k : keys) {
        auto it = dst->exifData().findKey(Exiv2::ExifKey(k));
        if (it != dst->exifData().end()) {
            dst->exifData().erase(it);
        }
    }
    Exiv2::ExifThumb thumb(dst->exifData());
    thumb.erase();
}


void Exiv2Metadata::import_exif_pairs(Exiv2::ExifData &out) const
{
    for (auto &p : exif_) {
        try {
            out[p.first] = p.second;
        } catch (std::exception &exc) {}
    }
}


void Exiv2Metadata::import_iptc_pairs(Exiv2::IptcData &out) const
{
    for (auto &p : iptc_) {
        try {
            auto &v = p.second;
            if (v.size() >= 1) {
                out[p.first] = v[0];
                for (size_t j = 1; j < v.size(); ++j) {
                    Exiv2::Iptcdatum d(Exiv2::IptcKey(p.first));
                    d.setValue(v[j]);
                    out.add(d);
                }
            }
        } catch (std::exception &exc) {}
    }
}


void Exiv2Metadata::saveToXmp(const Glib::ustring &path) const
{
    Exiv2::XmpData xmp;
    Exiv2::copyExifToXmp(exifData(), xmp);
    Exiv2::copyIptcToXmp(iptcData(), xmp);
    for (auto &datum : xmpData()) {
        xmp[datum.key()] = datum;
    }
    Exiv2::ExifData exif;
    Exiv2::IptcData iptc;
    import_exif_pairs(exif);
    import_iptc_pairs(iptc);
    Exiv2::copyExifToXmp(exif, xmp);
    Exiv2::copyIptcToXmp(iptc, xmp);

    std::string data;
    bool err = false;
    if (Exiv2::XmpParser::encode(data, xmp, Exiv2::XmpParser::omitPacketWrapper|Exiv2::XmpParser::useCompactFormat) != 0) {
        err = true;
    } else {
        FILE *out = g_fopen(path.c_str(), "wb");
        if (!out || fputs(data.c_str(), out) == EOF) {
            err = true;
        }
        if (out) {
            fclose(out);
        }
    }

    class Error: public Exiv2::AnyError {
    public:
        Error(const std::string &msg): msg_(msg) {}
        const char *what() const throw() { return msg_.c_str(); }
        int code() const throw() { return 0; }

    private:
        std::string msg_;
    };
    if (err) {
        throw Error("error saving XMP sidecar " + path);
    }
}


Glib::ustring Exiv2Metadata::xmpSidecarPath(const Glib::ustring &path)
{
    Glib::ustring fn = path;
    if (settings->xmp_sidecar_style == Settings::XmpSidecarStyle::STD) {
        fn = removeExtension(fn);
    }
    return fn + ".xmp";
}


Exiv2::XmpData Exiv2Metadata::getXmpSidecar(const Glib::ustring &path)
{
    Exiv2::XmpData ret;
    auto fname = xmpSidecarPath(path);
    if (Glib::file_test(fname, Glib::FILE_TEST_EXISTS)) {
        auto image = open_exiv2(fname);
        ret = image->xmpData();
    }
    return ret;
}


void Exiv2Metadata::init(const Glib::ustring &base_dir)
{
    cache_.reset(new ImageCache(IMAGE_CACHE_SIZE));
    exiftool_base_dir = base_dir;
    Exiv2::XmpParser::initialize();
}


void Exiv2Metadata::cleanup()
{
    Exiv2::XmpParser::terminate();
}

} // namespace rtengine
