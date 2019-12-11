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

// code taken from:
// 1. tmo_fattal02, and
//
// 2. the Matlab source accompanying the paper:
//
//   Paris, Hasinoff, and Kautz, "Local Laplacian Filters:
//   Edge-aware Image Processing with a Laplacian Pyramid", ACM 
//   Transactions on Graphics (Proc. SIGGRAPH 2011), 30(4), 2011.
//
// LICENSE:
//
// 1.
/* @author Grzegorz Krawczyk, <krawczyk@mpi-sb.mpg.de>
 *
 *
 * This file is a part of LuminanceHDR package, based on pfstmo.
 * ----------------------------------------------------------------------
 * Copyright (C) 2003,2004 Grzegorz Krawczyk
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 * ----------------------------------------------------------------------
 */
//
// 2. 
// Copyright (c) 2011 Sam Hasinoff

// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use,
// copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following
// conditions:

// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.

#ifdef _OPENMP
#include <omp.h>
#endif
#include <iostream>

#include "pyramids.h"
#include "rt_math.h"
#include "rescale.h"


namespace rtengine {

namespace {

void downSample (const array2D<float> &A, array2D<float> &B)
{
    const int width = B.width();
    const int height = B.height();

    // Note, I've uncommented all omp directives. They are all ok but are
    // applied to too small problems and in total don't lead to noticeable
    // speed improvements. The main issue is the pde solver and in case of the
    // fft solver uses optimised threaded fftw routines.
    //#pragma omp parallel for
    for ( int y = 0 ; y < height ; y++ ) {
        for ( int x = 0 ; x < width ; x++ ) {
            float p = A[2 * y][2 * x];
            p += A[2 * y][2 * x + 1];
            p += A[2 * y + 1][2 * x];
            p += A[2 * y + 1][2 * x + 1];
            B[y][x] = p * 0.25f; // p / 4.0f;
        }
    }
}


void upSample(const array2D<float> &A, array2D<float> &B)
{
    const int width = B.width();
    const int height = B.height();
    const int awidth = A.width();
    const int aheight = A.height();

    //#pragma omp parallel for shared(A, B)
    for ( int y = 0 ; y < height ; y++ ) {
        for ( int x = 0 ; x < width ; x++ ) {
            int ax = static_cast<int> (x * 0.5f); //x / 2.f;
            int ay = static_cast<int> (y * 0.5f); //y / 2.f;
            ax = (ax < awidth) ? ax : awidth - 1;
            ay = (ay < aheight) ? ay : aheight - 1;

            B[y][x] = A[ay][ax];
        }
    }
}


void gaussianBlur(float threshold, const array2D<float> &I, array2D<float> &L, bool multithread)
{
    const int width = I.width();
    const int height = I.height();
    const float tval = threshold;// / 2.f;

    if (width < 3 || height < 3) {
        if (&I != &L) {
            for (int y = 0; y < width; ++y) {
                for (int x = 0; x < height; ++x) {
                    L[y][x] = std::min(I[y][x], tval);
                }
            }
        }
        return;
    }

    array2D<float> T(width, height);

    //--- X blur
#ifdef _OPENMP
    #pragma omp parallel for shared(I, T) if(multithread)
#endif
    for (int y = 0; y < height; ++y) {
        int count = 0;
        float t = 0.f;
        
        for (int x = 1; x < width - 1; ++x) {
            t = 0.f;
            count = 0;
            
            if (I[y][x] < threshold) { 
                t = 2.f * I[y][x];
                count = 2;
            }
            if (I[y][x - 1] < threshold) {
                t += I[y][x - 1];
                ++count;
            }
            if (I[y][x + 1] < threshold) {
                t += I[y][x + 1];
                ++count;
            }
            if (!count) {
                T[y][x] = tval;//threshold;
            } else {
                T[y][x] = t / float(count);
            }
            // float t = 2.f * I[y][x];
            // t += I[y][x - 1];
            // t += I[y][x + 1];
            // T[y][x] = t * 0.25f; // t / 4.f;
        }

        t = 0.f;
        count = 0;
        if (I[y][0] < threshold) {
            t = 3.f * I[y][0];
            count = 3;
        }
        if (I[y][1] < threshold) {
            t += I[y][1];
            ++count;
        }
        if (count) {
            T[y][0] = t / float(count);
        } else {
            T[y][0] = tval;//threshold;
        }

        t = 0.f;
        count = 0;
        if (T[y][width - 1] < threshold) {
            t = 3.f * T[y][width - 1];
            count = 3;
        }
        if (I[y][width - 1] < threshold) {
            t += I[y][width - 1];
            ++count;
        }
        if (count) {
            T[y][width - 1] = t / float(count);
        } else {
            T[y][width - 1] = tval;//threshold;
        }
        // T[y][0] = (3.f * I[y][0] + I[y][1]) * 0.25f; // / 4.f;
        // T[y][width - 1] = (3.f * I[y][width - 1] + I[y][width - 2]) * 0.25f; // / 4.f;
    }

    //--- Y blur
#ifdef _OPENMP
    #pragma omp parallel for if(multithread)
#endif
    for (int x = 0 ;x < width - 7; x += 8) {
        float t = 0.f;
        int count = 0;
        
        for (int y = 1; y < height - 1; ++y) {
            for (int xx = 0; xx < 8; ++xx) {
                t = 0.f;
                count = 0;

                if (T[y][x + xx] < threshold) {
                    t = 2.f * T[y][x + xx];
                    count = 2;
                }
                if (T[y - 1][x + xx] < threshold) {
                    t += T[y - 1][x + xx];
                    ++count;
                }
                if (T[y + 1][x + xx] < threshold) {
                    t += T[y + 1][x + xx];
                    ++count;
                }
                if (count) {
                    L[y][x + xx] = t / float(count);
                } else {
                    L[y][x + xx] = tval;//threshold;
                }
                
                // float t = 2.f * T[y][x + xx];
                // t += T[y - 1][x + xx];
                // t += T[y + 1][x + xx];
                // L[y][x + xx] = t * 0.25f; // t/4.0f;
            }
        }

        for (int xx = 0; xx < 8; ++xx) {
            count = 0;
            t = 0.f;

            if (T[0][x + xx] < threshold) {
                t = 3.f * T[0][x + xx];
                count = 3;
            }
            if (T[1][x + xx] < threshold) {
                t += T[1][x + xx];
                ++count;
            }
            if (count) {
                L[0][x + xx] = t / float(count);
            } else {
                L[0][x + xx] = tval;//threshold;
            }

            count = 0;
            t = 0.f;
            if (T[height - 1][x + xx] < threshold) {
                t = 3.f * T[height - 1][x + xx];
                count = 3;
            }
            if (T[height - 2][x + xx] < threshold) {
                t += T[height - 2][x + xx];
                ++count;
            }
            if (count) {
                L[height - 1][x + xx] = t / float(count);
            } else {
                L[height - 1][x + xx] = tval;//threshold;
            }
            // L[0][x + xx] = (3.f * T[0][x + xx] + T[1][x + xx]) * 0.25f; // / 4.0f;
            // L[height - 1][x + xx] = (3.f * T[height - 1][x + xx] + T[height - 2][x + xx]) * 0.25f; // / 4.0f;
        }
    }

    for (int x = width - (width % 8); x < width; ++x) {
        int count = 0;
        float t = 0.f;
        
        for (int y = 1; y < height - 1; ++y) {
            count = 0;
            t = 0.f;

            if (T[y][x] < threshold) {
                t = 2.f * T[y][x];
                count = 2;
            }
            if (T[y - 1][x] < threshold) {
                t += T[y - 1][x];
                ++count;
            }
            if (T[y + 1][x] < threshold) {
                t += T[y + 1][x];
                ++count;
            }
            if (count) {
                L[y][x] = t / float(count);
            } else {
                L[y][x] = tval;//threshold;
            }
            // float t = 2.f * T[y][x];
            // t += T[y - 1][x];
            // t += T[y + 1][x];
            // L[y][x] = t * 0.25f; // t/4.0f;
        }

        count = 0;
        t = 0.f;
        if (T[0][x] < threshold) {
            t = 3.f * T[0][x];
            count = 3;
        }
        if (T[1][x] < threshold) {
            t += T[1][x];
            ++count;
        }
        if (count) {
            L[0][x] = t / float(count);
        } else {
            L[0][x] = tval;//threshold;
        }

        count = 0;
        t = 0.f;
        if (T[height - 1][x] < threshold) {
            t = 3.f * T[height - 1][x];
            count = 3;
        }
        if (T[height - 2][x] < threshold) {
            t += T[height - 2][x];
            ++count;
        }
        if (count) {
            L[height - 1][x] = t / float(count);
        } else {
            L[height - 1][x] = tval;//threshold;
        }

        // L[0][x] = (3.f * T[0][x] + T[1][x]) * 0.25f; // / 4.0f;
        // L[height - 1][x] = (3.f * T[height - 1][x] + T[height - 2][x]) * 0.25f; // / 4.0f;
    }
}

void copy(const array2D<float> &src, array2D<float> &dst, bool multithread)
{
    dst(src.width(), src.height());
    const int H = src.height();
    const int W = src.width();
    
#ifdef _OPENMP
    #pragma omp parallel for if(multithread)
#endif
    for (int y = 0; y < H; ++y) {
        for (int x = 0; x < W; ++x) {
            dst[y][x] = src[y][x];
        }
    }
}

} // namespace


void gaussianPyramid(float threshold, int nlevels, array2D<float> &src, Pyramid &out, bool multithread)
{
    int width = src.width();
    int height = src.height();

    out[0](src.width(), src.height());
    copy(src, out[0], multithread);

    array2D<float> *L = new array2D<float>(width, height);
    gaussianBlur(threshold, out[0], *L, multithread);

    for (int k = 1; k < nlevels; ++k) {
        if (width > 2 && height > 2) {
            width /= 2;
            height /= 2;
            out[k](width, height);
            downSample(*L, out[k]);
        } else {
            out[k](L->width(), L->height());
            copy(*L, out[k], multithread);
        }

        if (k < nlevels - 1) {
            delete L;
            L = new array2D<float>(width, height);
            gaussianBlur(threshold, out[k], *L, multithread);
        }
    }

    delete L;
}


void gaussianPyramid(int nlevels, array2D<float> &src, Pyramid &out, bool multithread)
{
    gaussianPyramid(RT_INFINITY, nlevels, src, out, multithread);
}


void laplacianPyramid(float threshold, int nlevels, array2D<float> &src, Pyramid &out, bool multithread)
{
    gaussianPyramid(threshold, nlevels, src, out, multithread);
    for (int k = 0; k < nlevels-1; ++k) {
        array2D<float> &cur = out[k];
        array2D<float> tmp(cur.width(), cur.height());
        upSample(out[k+1], tmp);

        const int h = cur.height();
        const int w = cur.width();
#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                cur[y][x] -= tmp[y][x];
            }
        }
    }
}


void laplacianPyramid(int nlevels, array2D<float> &src, Pyramid &out, bool multithread)
{
    laplacianPyramid(RT_INFINITY, nlevels, src, out, multithread);
}


void collapseLaplacianPyramid(int nlevels, Pyramid &pyramid, bool multithread)
{
    for (int k = nlevels - 2; k >= 0; --k) {
        array2D<float> tmp(pyramid[k].width(), pyramid[k].height());
        upSample(pyramid[k+1], tmp);

        const int h = tmp.height();
        const int w = tmp.width();
#ifdef _OPENMP
#       pragma omp parallel for if (multithread)
#endif
        for (int y = 0; y < h; ++y) {
            for (int x = 0; x < w; ++x) {
                pyramid[k][y][x] += tmp[y][x];
            }
        }
    }
}

} // namespace rtengine
