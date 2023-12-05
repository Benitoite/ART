/** -*- C++ -*-
 *  
 *  This file is part of ART.
 *
 *  Copyright (c) 2020 Alberto Griggio <alberto.griggio@gmail.com>
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

/*
 * HOW TO USE:
 *
 * Download the White Balance presets C file from darktable (in
 * src/external/wb_presets.c), put it in the same directory as this file,
 * compile and run. The program will generate the JSON data for WB presets on
 * stdout
 *
 */

#include <iostream>
#include <string.h>
#include <map>
#include <vector>
#include <sstream>
#include <time.h>
#include "wb_presets.c"


std::string key(const char *make, const char *model)
{
    std::ostringstream buf;
    for (const char *c = make; *c; ++c) {
        buf << char(toupper(*c));
    }
    buf << ' ';
    for (const char *c = model; *c; ++c) {
        buf << char(toupper(*c));
    }
    return buf.str();
}


int main()
{
    struct WBData {
        std::string label;
        std::array<double, 3> mult;
    };
    std::map<std::string, std::vector<WBData>> data;

    for (int i = 0; i < wb_preset_count; ++i) {
        const auto &p = wb_preset[i];
        if (p.tuning == 0 && p.channel[3] == 0) {
            data[key(p.make, p.model)].emplace_back(
            WBData{p.name, {p.channel[0], p.channel[1], p.channel[2]}});
        }
    }

    time_t t = time(nullptr);
    std::cout << "// White balance preset values for various cameras\n"
              << "// generated from src/external/wb_preset.c of Darktable\n"
              << "// see http://darktable.org for more information\n"
              << "// last updated on " << asctime(localtime(&t)) << std::endl;

    std::cout << "[\n";
    const char *sep = "";
    for (auto &p : data) {
        std::cout << sep;
        sep = ",\n";
        std::cout << "  {\n";
        std::cout << "    \"make_model\" : \"" << p.first << "\",\n";
        std::cout << "    \"presets\" : [\n";
        const char *sep2 = "";
        for (auto &w : p.second) {
            std::cout << sep2;
            sep2 = ",\n";
            std::cout << "      { \"name\" : \"" << w.label << "\","
                      << " \"multipliers\" : [ "
                      << w.mult[0] << ", " << w.mult[1] << ", "
                      << w.mult[2] << " ] }";
        }
        std::cout << "\n    ]\n  }";
    }
    std::cout << "\n]" << std::endl;
    
    return 0;
}
