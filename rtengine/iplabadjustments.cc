/* -*- C++ -*-
 *
 *  This file is part of RawTherapee.
 *
 *  Copyright 2019 Alberto Griggio <alberto.griggio@gmail.com>
 *
 *  RawTherapee is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  RawTherapee is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with RawTherapee.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef _OPENMP
#include <omp.h>
#endif

#include "improcfun.h"
#include "curves.h"

namespace rtengine {

void ImProcFunctions::labAdjustments(LabImage *lab, int pW, LUTu *histCCurve, LUTu *histLurve)
{
    LUTu dummy;
    LUTu hist16;
    LUTf lumacurve;
    LUTf clcurve;
    LUTf satcurve;
    LUTf lhskcurve;
    LUTf curve1;
    LUTf curve2;

    hist16(65536);
    lumacurve(32770, 0); // lumacurve[32768] and lumacurve[32769] will be set to 32768 and 32769 later to allow linear interpolation
    clcurve(65536, 0);
    satcurve(65536, 0);
    lhskcurve(65536, 0);
    curve1(65536);
    curve2(65536);

    if (params->labCurve.contrast != 0) { //only use hist16 for contrast
        hist16.clear();
        int fh = lab->H;
        int fw = lab->W;
        
#ifdef _OPENMP
#       pragma omp parallel if (multiThread)
#endif
        {
            LUTu hist16thr (hist16.getSize());  // one temporary lookup table per thread
            hist16thr.clear();
#ifdef _OPENMP
#           pragma omp for schedule(static) nowait
#endif
            for (int i = 0; i < fh; i++) {
                for (int j = 0; j < fw; j++) {
                    hist16thr[(int)((lab->L[i][j]))]++;
                }
            }
            
#ifdef _OPENMP
#           pragma omp critical
#endif
            {
                hist16 += hist16thr;
            }
        }
    }  

    bool utili;
    CurveFactory::complexLCurve(params->labCurve.brightness, params->labCurve.contrast, params->labCurve.lcurve, hist16, lumacurve, dummy, scale == 1 ? 1 : 16, utili);

    bool clcutili;
    CurveFactory::curveCL(clcutili, params->labCurve.clcurve, clcurve, scale == 1 ? 1 : 16);

    bool autili, butili;
    bool ccutili, cclutili;
    CurveFactory::complexsgnCurve(autili, butili, ccutili, cclutili, params->labCurve.acurve, params->labCurve.bcurve, params->labCurve.cccurve,
                                  params->labCurve.lccurve, curve1, curve2, satcurve, lhskcurve, scale == 1 ? 1 : 16);

    chromiLuminanceCurve(pW, lab, lab, curve1, curve2, satcurve, lhskcurve, clcurve, lumacurve, utili, autili, butili, ccutili, cclutili, clcutili, histCCurve ? *histCCurve : dummy, histLurve ? *histLurve : dummy);
}

} // namespace rtengine
