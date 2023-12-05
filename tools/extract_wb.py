#!/usr/bin/python3
# direct translation of extract_wb in python using as little external deps as possible

#############################################################################
# ART
"""
Adapted to ART from tools/extract_wb.py of darktable

    Copyright (C) 2009-2020 darktable developers.
    
    darktable is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    darktable is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    You should have received a copy of the GNU General Public License
    along with darktable.  If not, see <http://www.gnu.org/licenses/>.
"""
#############################################################################


from __future__ import print_function
import sys
from sys import argv
import os
import xml.etree.ElementTree as ET
import subprocess
from subprocess import PIPE
import shlex

#############################################################################
# ART
import argparse
import json
import time
from collections import OrderedDict as odict

def getopts():
    p = argparse.ArgumentParser()
    p.add_argument('-u', '--update',
                   help='update presets already found in file')
    p.add_argument('-f', '--force', action='store_true',
                   help='overwrite existing presets')
    p.add_argument('filenames', nargs='+')
    return p.parse_args()

opts = getopts()
argv = [None] + opts.filenames
#############################################################################

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

if len(argv) < 2 :
    sys.exit("Usage: extract_wb <file1> [file2] ...")

IGNORED_PRESETS = {"Auto", "Kelvin", "Measured", "AsShot", "As Shot", "Preset",
                 "Natural Auto", "Multi Auto", "Color Temperature Enhancement",
                 "One Touch WB 1", "One Touch WB 2", "One Touch WB 3",
                 "One Touch WB 4", "Custom WB 1", "Auto0", "Auto1", "Auto2",
                 "Custom", "CWB1", "CWB2", "CWB3", "CWB4", "Black", 
                 "Illuminator1", "Illuminator2", "Uncorrected"}

FL_PRESET_REPLACE = {
  "Fluorescent" : "CoolWhiteFluorescent",
  "FluorescentP1" : "DayWhiteFluorescent",
  "FluorescentP2" : "DaylightFluorescent",
  "FluorescentM1" : "WarmWhiteFluorescent",
  "FluorescentD"  : "DaylightFluorescent",
  "FluorescentN"  : "NeutralFluorescent",
  "FluorescentW"  : "WhiteFluorescent",
  "Daylight Fluorescent" : "DaylightFluorescent",
  "Day White Fluorescent" : "DayWhiteFluorescent",
  "White Fluorescent" : "WhiteFluorescent",
  "Unknown (0x600)" : "Underwater",
  "Sunny" : "DirectSunlight",
  "Fine Weather" : "DirectSunlight",
  "Tungsten (Incandescent)" : "Tungsten",
  "ISO Studio Tungsten" : "Tungsten",
  "Cool WHT FL" : "CoolWhiteFluorescent",
  "Daylight FL" : "DaylightFluorescent",
  "Warm WHT FL" : "WarmWhiteFluorescent",
  "Warm White Fluorescent" : "WarmWhiteFluorescent",
  "White FL" : "WhiteFluorescent",
  "Mercury Lamp" : "HighTempMercuryVaporFluorescent",
  "Day White FL" : "DayWhiteFluorescent",
  "Sodium Lamp" : "SodiumVaporFluorescent",
  "3000K (Tungsten light)" : "Tungsten",
  "4000K (Cool white fluorescent)" : "CoolWhiteFluorescent",
  "5300K (Fine Weather)" : "Daylight",
  "5500K (Flash)" : "Flash",
  "6000K (Cloudy)" : "Cloudy",
  "7500K (Fine Weather with Shade)" : "Shade",
  }

PRESET_ORDER = ["DirectSunlight", "Daylight", "D55", "Shade","Cloudy",
              "Tungsten", "Incandescent","Fluorescent", 
              "WarmWhiteFluorescent", "CoolWhiteFluorescent",
              "DayWhiteFluorescent","DaylightFluorescent",
              "DaylightFluorescent", "NeutralFluorescent", "WhiteFluorescent",
              "HighTempMercuryVaporFluorescent", "HTMercury", 
              "SodiumVaporFluorescent", "Underwater", "Flash", "Unknown"]

PRESET_SORT_MAPPING = {}

for index,name in enumerate(PRESET_ORDER):
    PRESET_SORT_MAPPING[name] = index + 1

##############################################################################
# ART
if False:
    cams_from_source = os.path.dirname(os.path.abspath(__file__)) + "/../src/external/rawspeed/data/cameras.xml"
    cams_from_dist = os.path.dirname(os.path.abspath(__file__)) + "/../rawspeed/cameras.xml"

    CAMERAS = os.path.abspath(cams_from_source) if os.path.exists(os.path.abspath(cams_from_source)) else os.path.abspath(cams_from_dist)

    if not os.path.exists(CAMERAS):
        sys.exit("Can't find cameras mapping file, should be in {0}".format(CAMERAS))

    exif_name_map = {}
    xml_doc = ET.parse(CAMERAS)
    for camera in xml_doc.getroot().findall('Camera'):
        maker = exif_maker = camera.get('make')
        model = exif_model = camera.get('model')
        exif_id = maker,model
        if camera.find('ID') is not None:
            cid = camera.find('ID')
            maker = cid.get('make')
            model = cid.get('model')
        exif_name_map[exif_id] = maker,model
        for alias in camera.findall('Alias'):
            exif_model = alias.text
            exif_id = exif_maker, exif_model
            exif_name_map[exif_id] = maker,model
else:
    class DummyDict(object):
        def __getitem__(self, key):
            make, model = key
            make = make.upper()
            model = model.upper()
            if model.startswith(make + ' '):
                model = model[len(make)+1:]
            return make, model
    exif_name_map = DummyDict()
############################################################################## 
    

found_presets = []

for filename in argv[1:]:
    red = green = blue = maker = model = preset = None
    finetune = fl_count = rlevel = blevel = glevel = 0
    listed_presets = []
    preset_names = {}
    gm_skew = False
    command = "exiftool -Make -Model \"-WBType*\" \"-WB_*\" \"-ColorTemp*\"    "\
    "-WhiteBalance -WhiteBalance2 -WhitePoint -ColorCompensationFilter       "\
    "-WBShiftAB -WBShiftAB_GM -WBShiftAB_GM_Precise -WBShiftGM -WBScale      "\
    "-WhiteBalanceFineTune -WhiteBalanceComp -WhiteBalanceSetting            "\
    "-WhiteBalanceBracket -WhiteBalanceBias -WBMode -WhiteBalanceMode        "\
    "-WhiteBalanceTemperature -WhiteBalanceDetected -ColorTemperature        "\
    "-WBShiftIntelligentAuto -WBShiftCreativeControl -WhiteBalanceSetup      "\
    "-WBRedLevel -WBBlueLevel -WBGreenLevel -RedBalance -BlueBalance         "\
    "\"{0}\"".format(filename)
    if filename.endswith(('.txt','.TXT')):
        command = 'cat "{0}"'.format(filename)
    command = shlex.split(command)
    proc = subprocess.check_output(command, universal_newlines=True)
    for io in proc.splitlines():
        lineparts = io.split(':')
        tag = lineparts[0].strip()
        values = lineparts[1].strip().split(' ')
        if 'Make' in tag.split():
            maker = lineparts[1].strip()
        elif 'Model' in tag.split():
            model = lineparts[1].strip()
        elif tag == "WB RGGB Levels":
            green = (float(values[1])+float(values[2]))/2.0
            red = float(values[0])/green
            blue = float(values[3])/green
            green = 1
        elif tag == "WB RB Levels":
            red = float(values[0])
            blue = float(values[1])
            if len(values) == 4 and values[2] == "256" and values[3] == "256":
                red /= 256.0
                blue /= 256.0
            green = 1
        elif tag == "WB GRB Levels":
            green = float(values[0])
            red = float(values[1])/green
            blue = float(values[2])/green
            green = 1
        # elif tag == "WB GRB Levels Auto" and maker == "FUJIFILM" # fuji seems to use "WB GRB Levels Auto to describe manual finetuning
        #  green = float(values[0])
        #  red = float(values[1])/green
        #  blue = float(values[2])/green
        #  green = 1
        elif tag == "White Point" and len(values) > 3:
            green = (float(values[1])+float(values[2]))/2.0
            red = float(values[0])/green
            blue = float(values[3])/green
            green = 1
        elif tag == "White Balance" or tag == "White Balance 2":
            preset = ' '.join(values)
            if preset in FL_PRESET_REPLACE:
                preset = FL_PRESET_REPLACE[preset]
        elif ' '.join(tag.split()[:2]) == "WB Type":
            preset_names[' '.join(tag.split()[2:])] = ' '.join(values)
        elif ' '.join(tag.split()[:3]) in ['WB RGB Levels', 'WB RGGB Levels', 'WB RB Levels']:
            # todo - this codepath is weird
            p = ''.join(tag.split()[3:])
            if( p in preset_names):
                p = preset_names[p]

            r=g=b=0

            if len(values) == 4 and ' '.join(tag.split()[:3]) in ['WB RB Levels']:
                g = (float(values[2])+float(values[3]))/2.0
                r = float(values[0])/g
                b = float(values[1])/g
                g = 1
            elif len(values) == 4:
                g = (float(values[1])+float(values[2]))/2.0
                r = float(values[0])/g
                b = float(values[3])/g
                g = 1
            elif len(values) == 3:
                g = float(values[1])
                r = float(values[0])/g
                b = float(values[2])/g
                g = 1
            elif len(values) == 2 and ' '.join(tag.split()[:3]) in ['WB RB Levels']:
                r = float(values[0])
                b = float(values[2])
                g = 1
            else:
                eprint("Found RGB tag '{0}' with {1} values instead of 2, 3 or 4".format(p, len(values)))
            
            if 'Fluorescent' in p:
                fl_count += 1
            
            if not p:
                p= 'Unknown'
            if p not in IGNORED_PRESETS:
                listed_presets.append(tuple([p,r,g,b]))
        elif tag == "WB Red Level":
            rlevel = float(values[0])
        elif tag == "WB Blue Level":
            blevel = float(values[0])
        elif tag == "WB Green Level":
            glevel = float(values[0])
        elif tag == "WB Shift AB": # canon - positive is towards amber, panasonic/leica/pentax - positive is towards blue?
            finetune = values[0]
        elif tag == "WB Shift GM": # detect GM shift and warn about it
            gm_skew = gm_skew or (int(values[0]) != 0)
        elif tag == "WB Shift AB GM": # sony
            finetune = values[0]
            gm_skew = gm_skew or (int(values[1]) != 0)
        elif tag == "White Balance Fine Tune" and maker.startswith("NIKON"): # nikon
            finetune = 0-(int(values[0]) * 2) # nikon lies about half-steps (eg 6->6->5 instead of 6->5.5->5, need to address this later on, so rescalling this now)
            gm_skew = gm_skew or (int(values[1]) != 0)
        elif tag == "White Balance Fine Tune" and maker == "FUJIFILM" and int(values[3]) != 0: # fuji
            eprint("Warning: Fuji does not seem to produce any sensible data for finetuning! If all finetuned values are identical, use one with no finetuning (0)")
            finetune = int(values[3]) / 20 # Fuji has -180..180 but steps are every 20
            gm_skew = gm_skew or (int(values[1].replace(',','')) != 0)
        elif tag == "White Balance Fine Tune" and maker == "SONY" and preset == "CoolWhiteFluorescent":
            # Sony's Fluorescent Fun
            if values[0] == "-1":
                preset = "WarmWhiteFluorescent"
            elif values[0] == "0":
                preset = "CoolWhiteFluorescent"
            elif values[0] == "1":
                preset = "DayWhiteFluorescent"
            elif values[0] == "2":
                preset = "DaylightFluorescent"
            else:
                eprint("Warning: Unknown Sony Fluorescent WB Preset!")
        elif tag == "White Balance Bracket": # olympus
            finetune = values[0]
            gm_skew = gm_skew or (int(values[1]) != 0)
        elif tag == "Color Compensation Filter": # minolta?
            gm_skew = gm_skew or (int(values[0]) != 0)

        if rlevel > 0 and glevel > 0 and blevel > 0:
            red = rlevel/glevel
            blue = blevel/glevel
            green = 1

    if gm_skew:
        eprint('WARNING: {0} has finetuning over GM axis! Data is skewed!'.format(filename))

    # Adjust the maker/model we found with the map we generated before
    if exif_name_map[maker,model]:
        enm = exif_name_map[maker,model]
        maker = enm[0]
        model = enm[1]
    else:
        eprint("WARNING: Couldn't find model in cameras.xml ('{0}', '{1}')".format(maker, model))

    for preset_arr in listed_presets:
        # ugly hack. Canon's Fluorescent is listed as WhiteFluorescent in usermanual
        preset_arrv = list(preset_arr)
        if maker and maker == "Canon" and preset_arrv[0] == "Fluorescent":
            preset_arrv[0] = "WhiteFluorescent"
        if preset_arrv[0] in FL_PRESET_REPLACE:
            preset_arrv[0] = FL_PRESET_REPLACE[preset_arrv[0]]
        if preset_arrv[0] not in IGNORED_PRESETS:
            found_presets.append(tuple([maker,model,preset_arrv[0], 0, preset_arrv[1], preset_arrv[2], preset_arrv[3]]))
    
    # Print out the WB value that was used in the file
    if not preset:
        preset = filename
    if red and green and blue and preset not in IGNORED_PRESETS:
        found_presets.append(tuple([maker, model, preset, int(finetune), red, green, blue]))

# get rid of duplicate presets

found_presets = list(set(found_presets))

def preset_to_sort(preset):
    sort_for_preset = 0
    if preset[2] in IGNORED_PRESETS:
        sort_for_preset = 0
    elif preset[2] in PRESET_SORT_MAPPING:
        sort_for_preset = PRESET_SORT_MAPPING[preset[2]]
    elif preset[2].endswith('K'):
        sort_for_preset = int(preset[2][:-1])
    else:
        eprint("WARNING: no defined sort order for '{0}'".format(preset[2]))
    return tuple([preset[0], preset[1], sort_for_preset, preset[3], preset[4], preset[5], preset[6]])

found_presets.sort(key=preset_to_sort)

min_padding = 0
for preset in found_presets:
    if len(preset[2]) > min_padding:
        min_padding = len(preset[2])

#dealing with Nikon half-steps
for index in range(len(found_presets)-1):
    if (found_presets[index][0] == 'Nikon' and #case now translated
        found_presets[index+1][0] == found_presets[index][0] and 
        found_presets[index+1][1] == found_presets[index][1] and
        found_presets[index+1][2] == found_presets[index][2] and
        found_presets[index+1][3] == found_presets[index][3]) :
       
        curr_finetune = int(found_presets[index][3])

        if curr_finetune < 0:
            found_presets[index+1] = list(found_presets[index+1])
            found_presets[index+1][3] = (int(found_presets[index+1][3]) + 1)
            found_presets[index+1] = tuple(found_presets[index+1])
        elif curr_finetune > 0:
            found_presets[index] = list(found_presets[index])
            found_presets[index][3] = (curr_finetune) - 1
            found_presets[index] = tuple(found_presets[index])

for index in range(len(found_presets)-1):
    if (found_presets[index][0] == "Nikon" and #case now translated
        found_presets[index+1][0] == found_presets[index][0] and ##
        found_presets[index+1][1] == found_presets[index][1] and
        found_presets[index+1][2] == found_presets[index][2]) :

        found_presets[index] = list(found_presets[index])
        found_presets[index+1] = list(found_presets[index+1])

        if (found_presets[index+1][3] % 2 == 0 and
            found_presets[index][3] % 2 == 0 and
            found_presets[index+1][3] == found_presets[index][3] + 2):
    
            #detected gap eg -12 -> -10. slicing in half to undo multiplication done earlier
            found_presets[index][3] = int(found_presets[index][3] / 2)
            found_presets[index+1][3] = int(found_presets[index+1][3] / 2)
        elif (found_presets[index+1][3] % 2 == 0 and
              found_presets[index][3] % 2 == 1 and
              found_presets[index+1][3] == (found_presets[index][3] + 1)*2 and
              (index + 2 == len(found_presets) or
               found_presets[index+2][2] != found_presets[index+1][2] ) ):

            #dealing with corner case of last-halfstep not being dealth with earlier
            found_presets[index+1][3] = int(found_presets[index+1][3] / 2)
        
        found_presets[index] = tuple(found_presets[index])
        found_presets[index+1] = tuple(found_presets[index+1])

#detect lazy finetuning (will not complain if there's no finetuning)
lazy_finetuning = []
for index in range(len(found_presets)-1):
    if (found_presets[index+1][0] == found_presets[index][0] and ##
        found_presets[index+1][1] == found_presets[index][1] and
        found_presets[index+1][2] == found_presets[index][2] and
        found_presets[index+1][3] != ((found_presets[index][3])+1) ):
        
        # found gap. complain about needing to interpolate
        lazy_finetuning.append(tuple([found_presets[index][0], found_presets[index][1], found_presets[index][2]]))

# Get rid of duplicate lazy finetuning reports
lazy_finetuning = list(set(lazy_finetuning))

# $stderr.puts lazy_finetuning.inspect.gsub("], ", "],\n") # debug content

for lazy in lazy_finetuning:
  eprint("Gaps detected in finetuning for {0} {1} preset {2}, dt will need to interpolate!".format(lazy[0], lazy[1], lazy[2]))

#############################################################################
# ART
if False:
    for preset in found_presets:
        if preset[2] in IGNORED_PRESETS:
            eprint("Ignoring preset '{0}'".format(preset[2]))
        else:
            preset_name = ''
            if preset[2].endswith('K'):
                preset_name = '"'+preset[2]+'"'
            else:
                preset_name = preset[2]
            print('  {{ "{0}", "{1}", {2:<{min_pad}}, {3}, {{ {4}, {5}, {6}, 0 }} }},'.format(preset[0], preset[1], preset_name, preset[3], preset[4], preset[5], preset[6], min_pad=min_padding))
else:
    IGNORED_PRESETS.add("Unknown")
    names = {
        "Daylight" : "daylight",
        "DirectSunlight" : "direct sunlight",
        "Cloudy" : "cloudy",
        "Shade" : "shade",
        "Incandescent" : "incandescent",
        "IncandescentWarm" : "incandescent warm",
        "Tungsten" : "tungsten",
        "Fluorescent" : "fluorescent",
        "FluorescentHigh" : "fluorescent high",
        "CoolWhiteFluorescent" : "cool white fluorescent",
        "WarmWhiteFluorescent" : "warm white fluorescent",
        "DaylightFluorescent" : "daylight fluorescent",
        "NeutralFluorescent" : "neutral fluorescent",
        "WhiteFluorescent" : "white fluorescent",
        "SodiumVaporFluorescent" : "sodium-vapor fluorescent",
        "DayWhiteFluorescent" : "day white fluorescent",
        "HighTempMercuryVaporFluorescent" : "high temp. mercury-vapor fluorescent",
        "HTMercury" : "high temp. mercury-vapor",
        "D55" : "D55",
        "Flash" : "flash",
        "FlashAuto" : "flash (auto mode)",
        "EveningSun" : "evening sun",
        "Underwater" : "underwater",
        "BlackNWhite" : "black & white",
        "uf_spot_wb" : "spot WB",
        "uf_manual_wb" : "manual WB",
        "uf_camera_wb" : "camera WB",
        "uf_auto_wb" : "auto WB",
        }

    comments = []
    db = odict()
    if opts.update:
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

    updated = set()
    for preset in found_presets:
        if preset[2] in IGNORED_PRESETS:
            eprint("Ignoring preset '{0}'".format(preset[2]))
        else:
            preset_name = ''
            ##if not preset[2].endswith('K'):
            ##     preset_name = '"'+preset[2]+'"'
            ## else:
            preset_name = names.get(preset[2], preset[2])
            make_model = "{0} {1}".format(preset[0], preset[1])
            data = db.get(make_model)
            entry = {"name" : preset_name,
                     "multipliers" : list(preset[4:7])}
            if not data:
                db[make_model] = {
                    "make_model" : make_model,
                    "presets" : [ entry ]
                    }
                updated.add(make_model)
            else:
                found = False
                for p in data["presets"]:
                    if p["name"] == preset_name:
                        if opts.force:
                            p["multipliers"] = entry["multipliers"]
                            updated.add(make_model)
                        found = True
                        break
                if not found:
                    data["presets"].append(entry)
                    updated.add(make_model)

    print("".join(comments))
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
#############################################################################
