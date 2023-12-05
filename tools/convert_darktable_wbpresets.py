import json
import sys


def main():
    db = json.load(sys.stdin)
    res = []
    assert 'version' in db and db['version'] == 1
    presets = db['wb_presets']
    for entry in presets:
        make = entry['maker']
        for data in entry['models']:
            camname = make + ' ' + data['model']
            out = { 'make_model' : camname.upper() }
            presets = []
            for p in data['presets']:
                if 'tuning' not in p:
                    name = p['name'].lower()
                    mults = p['channels'][:-1]
                    presets.append({'name' : name, 'multipliers' : mults})
            out['presets'] = presets
            res.append(out)
    json.dump(res, sys.stdout, indent=2)
    sys.stdout.write('\n')


if __name__ == '__main__':
    main()

