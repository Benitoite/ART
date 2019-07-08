from __future__ import print_function
import os
import subprocess

def get_keys():
    keys = {}
    with open('rtdata/languages/default') as f:
        for line in f:
            if line.strip().startswith('#'):
                continue
            bits = line.split(';', 1)
            assert len(bits) == 2, bits
            if bits:
                keys[bits[0]] = bits[1].strip()
    return keys


def get_files():
    p = subprocess.Popen(['hg', 'manifest'], stdout=subprocess.PIPE)
    out, _ = p.communicate()
    ret = {}
    for name in out.splitlines():
        if name.startswith('rtgui/'):
            with open(name) as f:
                ret[name] = f.read()
    return ret


def find_match(key, files):
    for name in files:
        if key in files[name]:
            return True
    return False


def main():
    os.chdir(os.path.join(os.path.dirname(__file__), '..'))
    keys = get_keys()
    files = get_files()
    ret = []
    for key in keys:
        if find_match(key, files):
            ret.append(key)
    for key in sorted(ret):
        print('%s;%s' % (key, keys[key]))


if __name__ == '__main__':
    main()
