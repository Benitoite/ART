/* -*- C++ -*-
 *
 *  This file is part of ART.
 *
 *  Copyright (c) 2019 Alberto Griggio <alberto.griggio@gmail.com>
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

#include <vector>
#include <memory>
#include "array2D.h"


namespace rtengine {

class Pyramid {
public:
    explicit Pyramid(int nlevels):
        size_(nlevels)
    {
        data_ = new array2D<float>*[size_];
        for (size_t i = 0; i < size_; ++i) {
            data_[i] = new array2D<float>();
        }
    }

    ~Pyramid()
    {
        for (size_t i = 0; i < size_; ++i) {
            delete data_[i];
        }
        delete[] data_;
    }
    
    size_t size() const { return size_; }

    array2D<float> &operator[](size_t k)
    {
        return *(data_[k]);
    }

private:
    array2D<float> **data_;
    size_t size_;
};

void gaussianPyramid(int nlevels, array2D<float> &src, Pyramid &out, bool multithread);
void gaussianPyramid(float threshold, int nlevels, array2D<float> &src, Pyramid &out, bool multithread);

void laplacianPyramid(int nlevels, array2D<float> &src, Pyramid &out, bool multithread);
void laplacianPyramid(float threshold, int nlevels, array2D<float> &src, Pyramid &out, bool multithread);

void collapseLaplacianPyramid(int nlevels, Pyramid &pyramid, bool multithread);

} // namespace rtengine
