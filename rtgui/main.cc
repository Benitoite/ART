/*
 *  This file is part of ART.
 *
 *  Copyright (c) 2004-2010 Gabor Horvath <hgabor@rawtherapee.com>
 *
 *  ART is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  ART is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with ART.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef __GNUC__
#if defined(__FAST_MATH__)
#error Using the -ffast-math CFLAG is known to lead to problems. Disable it to compile ART.
#endif
#endif

#include "config.h"
#include <gtkmm.h>
#include <giomm.h>
#include <iostream>
#include <tiffio.h>
#include "rtwindow.h"
#include <cstring>
#include <cstdlib>
#include <locale.h>
#include <lensfun.h>
#include "options.h"
#include "soundman.h"
#include "rtimage.h"
#include "version.h"
#include "extprog.h"
#include "../rtengine/dynamicprofile.h"
#include "printhelp.h"
#include "wbprovider.h"
#include "pathutils.h"
#include "session.h"

#ifndef WIN32
#include <glibmm/fileutils.h>
#include <glib.h>
#include <glib/gstdio.h>
#include <glibmm/threads.h>
#else
#include <glibmm/thread.h>
#include "conio.h"
#endif

#ifdef WITH_MIMALLOC
#  include <mimalloc.h>
#endif

// Set this to 1 to make RT work when started with Eclipse and arguments, at least on Windows platform
#define ECLIPSE_ARGS 0

extern Options options;

// stores path to data files
Glib::ustring argv0;
Glib::ustring creditsPath;
Glib::ustring licensePath;
Glib::ustring argv1;
Glib::ustring argv2;
bool simpleEditor = false;
bool gimpPlugin = false;
bool remote = false;
bool is_session = false;
unsigned char initialGdkScale = 1;

namespace {

// This recursive mutex will be used by gdk_threads_enter/leave instead of a simple mutex
static Glib::Threads::RecMutex myGdkRecMutex;

static void myGdkLockEnter()
{
    myGdkRecMutex.lock();
}

static void myGdkLockLeave()
{
    myGdkRecMutex.unlock();
}


void process_help_params(int argc, char **argv)
{
    for (int i = 1; i < argc; ++i) {
        Glib::ustring currParam(argv[i]);
#if ECLIPSE_ARGS
        currParam = currParam.substr(1, currParam.length() - 2);
#endif
        if (currParam.length() > 1 && currParam[0] == '-') {
            switch (currParam[1]) {
            case 'v':
                std::cout << RTNAME << ", version " << RTVERSION << std::endl;
                exit(0);
            case '?':
            case 'h':
                ART_print_help(argv[0], true);
                exit(0);
            default:
                break;
            }
        }
    }
}


/* Process line command options
 * Returns
 *  0 if process in batch has executed
 *  1 to start GUI (with a dir or file option)
 *  2 to start GUI because no files found
 *  -1 if there is an error in parameters
 *  -2 if an error occurred during processing
 *  -3 if at least one required procparam file was not found */
int processLineParams(int argc, char **argv)
{
    Glib::setenv("ART_IS_SESSION", "0");
    int ret = 1;
    for (int iArg = 1; iArg < argc; iArg++) {
        Glib::ustring currParam (argv[iArg]);
        if (currParam.empty()) {
            continue;
        }
#if ECLIPSE_ARGS
        currParam = currParam.substr(1, currParam.length() - 2);
#endif

        if (currParam[0] == '-' && currParam.size() > 1) {
            switch (currParam[1]) {
#ifdef WIN32

            case 'w': // This case is handled outside this function
                break;
#endif

            case 'v':
                printf("%s, version %s\n", RTNAME, RTVERSION);
                ret = 0;
                break;

#ifndef __APPLE__ // TODO agriggio - there seems to be already some "single instance app" support for OSX in rtwindow. Disabling it here until I understand how to merge the two

            case 'R':
                if (!gimpPlugin) {
                    remote = true;
                }

                break;
#endif
            case 's':
                simpleEditor = true;
                remote = false;
                break;
                
            case '-':
                if (currParam.substr(5) == "--gtk" || currParam == "--g-fatal-warnings") {
                    break;
                }

            case 'g':
                if (currParam == "-gimp") {
                    gimpPlugin = true;
                    simpleEditor = true;
                    remote = false;
                    break;
                }

            case 'S':
                if (remote) {
                    Glib::setenv("ART_IS_SESSION", "1");
                    if (currParam == "-S") {
                        if (iArg+1 < argc && argv[iArg+1][0] != '-') {
                            ++iArg;
                            art::session::load(Glib::ustring(fname_to_utf8(argv[iArg])));
                        } else {
                            art::session::load(art::session::filename() + ".last");
                        }
                        break;
                    } else if (currParam == "-Sc") {
                        art::session::clear();
                        break;
                    } else if (currParam == "-Sa" || currParam == "-Sr") {
                        std::vector<Glib::ustring> fnames;
                        ++iArg;
                        while (iArg < argc && argv[iArg][0] != '-') {
                            fnames.emplace_back(fname_to_utf8(argv[iArg]));
                            ++iArg;
                        }
                        if (currParam == "-Sa") {
                            art::session::add(fnames);
                        } else {
                            art::session::remove(fnames);
                        }
                        break;
                    }
                }
                // no break here on purpose

            case 'h':
            case '?':
            default:
                ART_print_help(argv[0], true);
                ret = -1;
                break;
            }
        } else {
            if (argv1.empty()) {
                argv1 = Glib::ustring(fname_to_utf8(argv[iArg]));
#if ECLIPSE_ARGS
                argv1 = argv1.substr(1, argv1.length() - 2);
#endif
            } else if (gimpPlugin) {
                argv2 = Glib::ustring(fname_to_utf8(argv[iArg]));
                break;
            }

            if (!gimpPlugin) {
                break;
            }
        }
    }

    return ret;
}


bool init_rt()
{
    UserCommandStore::getInstance()->init(Glib::build_filename(options.rtdir, "usercommands"));
    SoundManager::init();
    wb_presets::init(argv0, options.rtdir);

    if ( !options.rtSettings.verbose ) {
        TIFFSetWarningHandler (nullptr);   // avoid annoying message boxes
    }

#ifndef WIN32

    // Move the old path to the new one if the new does not exist
    if (Glib::file_test (Glib::build_filename (options.rtdir, "cache"), Glib::FILE_TEST_IS_DIR) && !Glib::file_test (options.cacheBaseDir, Glib::FILE_TEST_IS_DIR)) {
        g_rename (Glib::build_filename (options.rtdir, "cache").c_str (), options.cacheBaseDir.c_str ());
    }

#endif

    return true;
}


void cleanup_rt()
{
    rtengine::cleanup();
}


RTWindow *create_rt_window()
{
    Glib::ustring icon_path = Glib::build_filename (argv0, "images");
    Glib::RefPtr<Gtk::IconTheme> defaultIconTheme = Gtk::IconTheme::get_default();
    defaultIconTheme->append_search_path (icon_path);

    //gdk_threads_enter ();
    RTWindow *rtWindow = new RTWindow();
    return rtWindow;
}


class RTApplication: public Gtk::Application
{
public:
    RTApplication():
        Gtk::Application("us.pixls.art.application",
                         Gio::APPLICATION_SEND_ENVIRONMENT|Gio::APPLICATION_HANDLES_COMMAND_LINE),
        rtWindow (nullptr)
    {
    }

    ~RTApplication() override
    {
        if (rtWindow) {
//            delete rtWindow;
            art::session::save(art::session::filename() + ".last");
            art::session::clear();
            cleanup_rt();
        }
    }
    
private:
    bool create_window()
    {
        if (rtWindow) {
            return true;
        }

        if (!init_rt()) {
            Gtk::MessageDialog msgd ("Fatal error!\nThe ART_SETTINGS environment variable is set, but uses a relative path. The path must be absolute!", true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
            add_window (msgd);
            msgd.run ();
            return false;
        } else {
            rtWindow = create_rt_window();
            add_window(*rtWindow);
            return true;
        }
    }

    int on_command_line(const Glib::RefPtr<Gio::ApplicationCommandLine> &command_line) override
    {
        auto s = command_line->getenv("ART_IS_SESSION");
        if (!s.empty() && atoi(s.c_str())) {
            is_session = true;
        }
        int argc = 0;
        auto argv = command_line->get_arguments(argc);
        if (is_session) {
            bool raise = !rtWindow;
            if (create_window()) {
                const auto doit =
                    [] (gpointer data) -> gboolean {
                        FileCatalog *filecatalog = static_cast<FileCatalog *>(data);
                        if (!art::session::check(filecatalog->lastSelectedDir())) {
                            filecatalog->dirSelected(art::session::path(), "");
                        }
                        return FALSE;
                    };
                gdk_threads_add_idle(doit, rtWindow->fpanel->fileCatalog);
                if (raise) {
                    rtWindow->present();
                }
                return 0;
            }
        } else if (argc == 2) {
            if (create_window()) {
                struct Data {
                    Glib::ustring dirname;
                    Glib::ustring fname;
                    FileCatalog *filecatalog;
                };
                Data *d = new Data;
                d->fname = argv[1];
                if (Glib::file_test(d->fname, Glib::FILE_TEST_IS_DIR)) {
                    d->dirname = d->fname;
                    d->fname = "";
                } else {
                    d->dirname = Glib::path_get_dirname(d->fname);
                }
                d->filecatalog = rtWindow->fpanel->fileCatalog;

                const auto doit =
                    [] (gpointer data) -> gboolean {
                        Data *d = static_cast<Data *>(data);
                        d->filecatalog->dirSelected(d->dirname, d->fname);
                        delete d;
                        return FALSE;
                    };
                gdk_threads_add_idle(doit, d);
                rtWindow->present();
                return 0;
            }
        } else {
            if (create_window()) {
                rtWindow->present();
                return 0;
            }
        }
        return 1;
    }

private:
    RTWindow *rtWindow;
};

void show_gimp_plugin_info_dialog(Gtk::Window *parent)
{
    if (options.gimpPluginShowInfoDialog) {
        Gtk::MessageDialog info(*parent, M("GIMP_PLUGIN_INFO"), false, Gtk::MESSAGE_INFO, Gtk::BUTTONS_OK, true);
        Gtk::Box *box = info.get_message_area();
        Gtk::CheckButton dontshowagain(M("DONT_SHOW_AGAIN"));
        dontshowagain.show();
        box->pack_start(dontshowagain);
        info.run();
        options.gimpPluginShowInfoDialog = !dontshowagain.get_active();
    }
}


} // namespace


int main (int argc, char **argv)
{
#ifdef WITH_MIMALLOC
    mi_version();
#endif
    
    setlocale (LC_ALL, "");
    setlocale (LC_NUMERIC, "C"); // to set decimal point to "."

    simpleEditor = false;
    gimpPlugin = false;
    remote = true;
    argv0 = "";
    argv1 = "";
    argv2 = "";

    process_help_params(argc, argv);

    Glib::init();  // called by Gtk::Main, but this may be important for thread handling, so we call it ourselves now
    Gio::init();

#ifdef WIN32
    if (GetFileType (GetStdHandle (STD_OUTPUT_HANDLE)) == 0x0003) {
        // started from msys2 console => do not buffer stdout
        setbuf(stdout, NULL);
    }
#endif

#ifdef BUILD_BUNDLE
    char exname[512] = {0};
    Glib::ustring exePath;
    // get the path where the rawtherapee executable is stored
#ifdef WIN32
    WCHAR exnameU[512] = {0};
    GetModuleFileNameW (NULL, exnameU, 511);
    WideCharToMultiByte (CP_UTF8, 0, exnameU, -1, exname, 511, 0, 0 );
#else

    if (readlink ("/proc/self/exe", exname, 511) < 0) {
        strncpy (exname, argv[0], 511);
    }

#endif
    exePath = Glib::path_get_dirname (exname);

    // set paths
    if (Glib::path_is_absolute (DATA_SEARCH_PATH)) {
        argv0 = DATA_SEARCH_PATH;
    } else if (strcmp(DATA_SEARCH_PATH, ".") == 0) {
        argv0 = exePath;
    } else {
        argv0 = Glib::build_filename (exePath, DATA_SEARCH_PATH);
    }

    if (Glib::path_is_absolute (CREDITS_SEARCH_PATH)) {
        creditsPath = CREDITS_SEARCH_PATH;
    } else {
        creditsPath = Glib::build_filename (exePath, CREDITS_SEARCH_PATH);
    }

    if (Glib::path_is_absolute (LICENCE_SEARCH_PATH)) {
        licensePath = LICENCE_SEARCH_PATH;
    } else {
        licensePath = Glib::build_filename (exePath, LICENCE_SEARCH_PATH);
    }

    options.rtSettings.lensfunDbDirectory = LENSFUN_DB_PATH;

#else
    argv0 = DATA_SEARCH_PATH;
    creditsPath = CREDITS_SEARCH_PATH;
    licensePath = LICENCE_SEARCH_PATH;
    options.rtSettings.lensfunDbDirectory = LENSFUN_DB_PATH;
#endif

    Glib::ustring fatalError;
    try {
        Options::load();
    } catch (Options::Error &e) {
        fatalError = e.get_msg();
    }

#ifdef WIN32
    bool consoleOpened = false;

    // suppression of annoying error boxes
    SetErrorMode (SEM_FAILCRITICALERRORS | SEM_NOGPFAULTERRORBOX | SEM_NOOPENFILEERRORBOX);

    if (argc > 1) {
        if (!remote && !Glib::file_test (argv1, Glib::FILE_TEST_EXISTS ) && !Glib::file_test (argv1, Glib::FILE_TEST_IS_DIR)) {
            bool stdoutRedirecttoConsole = (GetFileType (GetStdHandle (STD_OUTPUT_HANDLE)) == 0x0000);
            // open console, if stdout is invalid
            if (stdoutRedirecttoConsole) {
                // check if parameter -w was passed.
                // We have to do that in this step, because it controls whether to open a console to show the output of following steps
                bool Console = true;

                for (int i = 1; i < argc; i++)
                    if (!strcmp (argv[i], "-w") || !strcmp (argv[i], "-R") || !strcmp (argv[i], "-gimp")) {
                        Console = false;
                        break;
                    }

                if (Console && AllocConsole()) {
                    AttachConsole ( GetCurrentProcessId() ) ;
                    // Don't allow CTRL-C in console to terminate RT
                    SetConsoleCtrlHandler ( NULL, true );
                    // Set title of console
                    char consoletitle[128];
                    sprintf (consoletitle, "%s %s Console", RTNAME, RTVERSION);
                    SetConsoleTitle (consoletitle);
                    // increase size of screen buffer
                    COORD c;
                    c.X = 200;
                    c.Y = 1000;
                    SetConsoleScreenBufferSize ( GetStdHandle ( STD_OUTPUT_HANDLE ), c );
                    // Disable console-Cursor
                    CONSOLE_CURSOR_INFO cursorInfo;
                    cursorInfo.dwSize = 100;
                    cursorInfo.bVisible = false;
                    SetConsoleCursorInfo ( GetStdHandle ( STD_OUTPUT_HANDLE ), &cursorInfo );

                    if (stdoutRedirecttoConsole) { // if stdout is Redirect to console, we also redirect stderr to console
                        freopen ( "CON", "w", stdout ) ;
                        freopen ( "CON", "w", stderr ) ;
                    }

                    freopen ( "CON", "r", stdin ) ;

                    consoleOpened = true;
                }
            }
        }
        int ret = processLineParams ( argc, argv);

        if ( ret <= 0 ) {
            fflush(stdout);
            if (consoleOpened) {
                printf ("Press any key to exit " RTNAME "\n");
                FlushConsoleInputBuffer (GetStdHandle (STD_INPUT_HANDLE));
                getch();
            }

            return ret;
        }
    }

#else

    if (argc > 1) {
        int ret = processLineParams ( argc, argv);

        if ( ret <= 0 ) {
            return ret;
        }
    }

#endif

    if (gimpPlugin) {
        if (!Glib::file_test (argv1, Glib::FILE_TEST_EXISTS) || Glib::file_test (argv1, Glib::FILE_TEST_IS_DIR)) {
            printf ("Error: argv1 doesn't exist\n");
            return 1;
        }

        if (argv2.empty()) {
            printf ("Error: -gimp requires two arguments\n");
            return 1;
        }
    } else if (!remote && !(Glib::file_test(argv1, Glib::FILE_TEST_EXISTS) && !Glib::file_test(argv1, Glib::FILE_TEST_IS_DIR))) {
        simpleEditor = false;
    }

    int ret = 0;

    if (options.pseudoHiDPISupport) {
        // Reading/updating GDK_SCALE early if it exists
        const gchar *gscale = g_getenv("GDK_SCALE");
        if (gscale && gscale[0] == '2') {
            initialGdkScale = 2;
        }
        // HOMBRE: On Windows, if resolution is set to 200%, Gtk internal variables are SCALE=2 and DPI=96
        g_setenv("GDK_SCALE", "1", true);
    }

    gdk_threads_set_lock_functions (G_CALLBACK (myGdkLockEnter), (G_CALLBACK (myGdkLockLeave)));
    gdk_threads_init();
    gtk_init (&argc, &argv);  // use the "--g-fatal-warnings" command line flag to make warnings fatal

    if (fatalError.empty() && remote) {
        char *app_argv[2] = { const_cast<char *> (argv0.c_str()) };
        int app_argc = 1;

        if (!argv1.empty()) {
            app_argc = 2;
            app_argv[1] = const_cast<char *> (argv1.c_str());
        }

        RTApplication app;
        ret = app.run(app_argc, app_argv);
    } else {
        if (fatalError.empty() && init_rt()) {
            Gtk::Main m (&argc, &argv);
            gdk_threads_enter();
            const std::unique_ptr<RTWindow> rtWindow (create_rt_window());
            if (gimpPlugin) {
                show_gimp_plugin_info_dialog(rtWindow.get());
            }
            m.run (*rtWindow);
            gdk_threads_leave();

            if (gimpPlugin && rtWindow->epanel && rtWindow->epanel->isRealized()) {
                if (!rtWindow->epanel->saveImmediately(argv2, SaveFormat())) {
                    ret = -2;
                }
            }

            cleanup_rt();
        } else {
            Gtk::Main m (&argc, &argv);
            Gtk::MessageDialog msgd (Glib::ustring::compose("FATAL ERROR!\n\n%1", fatalError), true, Gtk::MESSAGE_ERROR, Gtk::BUTTONS_OK, true);
            msgd.run ();
            ret = -2;
        }
    }

#ifdef WIN32

    if (consoleOpened) {
        printf ("Press any key to exit " RTNAME "\n");
        fflush(stdout);
        FlushConsoleInputBuffer (GetStdHandle (STD_INPUT_HANDLE));
        getch();
    }

#endif

    return ret;
}
