/* -*- C++ -*-
 *
 *  This file is part of ART.
 *
 *  Copyright 2023 Alberto Griggio <alberto.griggio@gmail.com>
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
#pragma once

#include <string>
#include <glibmm/ustring.h>
#include <vector>

namespace rtengine {

enum class CLUTParamType {
    PT_INT,
    PT_FLOAT,
    PT_BOOL,
    PT_CHOICE
};

struct CLUTParamDescriptor {
    std::string name;
    CLUTParamType type;
    double value_min;
    double value_max;
    double value_default;
    std::vector<Glib::ustring> choices;
    Glib::ustring gui_name;
    Glib::ustring gui_group;
    double gui_step;
};

} // namespace rtengine
