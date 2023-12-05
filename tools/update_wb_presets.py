#!/usr/bin/env python

import argparse
import json
import time
from collections import OrderedDict as odict


def getopts():
    p = argparse.ArgumentParser()
    p.add_argument('-u', '--update', required=True,
                   help='wb presets to update')
    p.add_argument('filename')
    p.add_argument('-f', '--force', action='store_true',
                   help='overwrite existing presets')
    return p.parse_args()


def main():
    opts = getopts()
    comments = []
    db = odict()
    with open(opts.update) as f:
        data = []
        for line in f:
            if not line.startswith('//'):
                data.append(line)
            else:
                comments.append(line)
        data = json.loads("".join(data))
    for entry in data:
        db[entry["make_model"]] = entry

    with open(opts.filename) as f:
        data = []
        for line in f:
            if not line.startswith('//'):
                data.append(line)
        data = json.loads("".join(data))

    updated = []
    for entry in data:
        key = entry['make_model']
        if opts.force or key not in db:
            db[key] = entry
            updated.append(key)
    
    print("".join(comments))
    if updated:
        print('// updated on %s:' % time.asctime())
        for make_model in sorted(updated):
            print('// %s' % make_model)
    print("\n[")
    sep = ""
    for key in db:
        entry = db[key]
        print(sep, end='')
        sep = ',\n'
        print('  {\n    "make_model" : "%s",' % entry['make_model'])
        print('    "presets" : [')
        print('      ', end='')
        print(',\n      '.join('{ "name" : "%s", "multipliers" : [ %s ] }' %
                                (p['name'],
                                 ', '.join(map(str, p['multipliers'])))
                               for p in entry['presets']))
        print('    ]\n  }', end='')
    print("\n]")


if __name__ == '__main__':
    main()
