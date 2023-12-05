#!/usr/bin/env python

import os, sys, argparse
import PyInstaller.__main__


def getopts():
    p = argparse.ArgumentParser()
    p.add_argument('imageio_dir')
    p.add_argument('-o', '--outdir', default='.')
    ret = p.parse_args()
    ret.imageio_dir = os.path.normpath(os.path.abspath(ret.imageio_dir))
    return ret


def main():
    mydir = os.path.abspath(os.path.dirname(__file__))
    opts = getopts()
    hidden_imports = []
    for (dirpath, dirnames, filenames) in os.walk(opts.imageio_dir):
        for name in filenames:
            if name.endswith('.py'):
                sys.path.append(dirpath)
                m = os.path.splitext(name)[0]
                hidden_imports.append('--hidden-import=' + m)
                
    if opts.outdir:
        os.chdir(opts.outdir)
        print(';; directory is: %s' % os.getcwd())

    sep = os.pathsep
    tool = os.path.join(mydir, 'imageio_driver.py')

    args = ['--name=python', '--clean'] + hidden_imports + [tool]
    PyInstaller.__main__.run(args)    


if __name__ == '__main__':
    main()
