#!/usr/bin/env python3

import os, sys
import shutil
import subprocess
import argparse

def getdlls(opts):
    res = []
    d = os.getcwd()
    p = subprocess.Popen(['ldd', os.path.join(d, 'ART.exe')],
                         stdout=subprocess.PIPE)
    out, _ = p.communicate()
    for line in out.decode('utf-8').splitlines():
        if ' => ' in line:
            bits = line.split()
            lib = bits[2]
            if not lib.lower().startswith('/c/windows/'):
                res.append(opts.msys + lib)
    return res

def getopts():
    p = argparse.ArgumentParser()
    p.add_argument('-o', '--outdir', required=True,
                   help='output directory for the bundle')
    p.add_argument('-m', '--msys', default='c:/msys64',
                   help='msys installation directory (default: c:/msys64)')
    p.add_argument('-e', '--exiftool',
                   help='path to exiftool.exe (default: search in PATH)')
    p.add_argument('-v', '--verbose', action='store_true')
    ret = p.parse_args()
    if ret.exiftool is None:
        ret.exiftool = shutil.which('exiftool')
    if ret.msys.endswith('/'):
        ret.msys = ret.msys[:-1]
    return ret

def extra_files(opts):
    def D(s): return opts.msys + '/' + s
    return [
        ('.', [
            D('mingw64/bin/gdbus.exe'),
            D('mingw64/bin/gspawn-win64-helper.exe'),
            D('mingw64/bin/gspawn-win64-helper-console.exe')
            ] + ([(opts.exiftool, 'exiftool.exe')] if opts.exiftool else []),
        ),
        ('share/icons/Adwaita', [
            D('mingw64/share/icons/Adwaita/scalable'),
            D('mingw64/share/icons/Adwaita/index.theme'), 
            D('mingw64/share/icons/Adwaita/cursors'),
        ]),
        ('lib', [
            D('mingw64/lib/gdk-pixbuf-2.0'),
        ]),
        ('share/glib-2.0/schemas', [
            D('mingw64/share/glib-2.0/schemas/gschemas.compiled'),
        ]),
        ('share', [
            (D('mingw64/var/lib/lensfun-updates/version_1'), 'lensfun'),
        ]),
    ]


def main():
    opts = getopts()
    d = os.getcwd()
    if not os.path.exists('ART.exe'):
        sys.stderr.write('ERROR: ART.exe not found! Please run this script '
                         'from the build directory of ART\n')
        sys.exit(1)
    if opts.verbose:
        print('copying %s to %s' % (os.getcwd(), opts.outdir))
    shutil.copytree(d, opts.outdir)
    for lib in getdlls(opts):
        if opts.verbose:
            print('copying: %s' % lib)
        shutil.copy2(lib,
                     os.path.join(opts.outdir, os.path.basename(lib)))
    for key, elems in extra_files(opts):
        for elem in elems:
            name = None
            if isinstance(elem, tuple):
                elem, name = elem
            else:
                name = os.path.basename(elem)
            if opts.verbose:
                print('copying: %s' % elem)
            if os.path.isdir(elem):
                shutil.copytree(elem, os.path.join(opts.outdir, key, name))
            else:
                dest = os.path.join(opts.outdir, key, name)
                destdir = os.path.dirname(dest)
                if not os.path.exists(destdir):
                    os.makedirs(destdir)
                shutil.copy2(elem, dest)
    os.makedirs(os.path.join(opts.outdir, 'share/gtk-3.0'))
    with open(os.path.join(opts.outdir, 'share/gtk-3.0/settings.ini'), 'w') \
         as out:
        out.write('[Settings]\ngtk-button-images=1\n')


if __name__ == '__main__':
    main()
