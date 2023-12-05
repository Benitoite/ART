/* -*- C -*-
 * 
 * Simplified porting of darktable's sigmoid module to ART CTL
 * Copyright of the original code follows
 */
/*
    This file is part of darktable,
    Copyright (C) 2020-2023 darktable developers.

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
*/

// @ART-label: "$CTL_SIGMOID"
// @ART-colorspace: "rec2020"

//-----------------------------------------------------------------------------
// constants for color space conversions and display mapping
//-----------------------------------------------------------------------------

const float xyz_rec2020[3][3] = {
    {0.6734241,  0.1656411,  0.1251286},
    {0.2790177,  0.6753402,  0.0456377},
    { -0.0019300,  0.0299784, 0.7973330}
};

const float rec709_xyz[3][3] = {
    {3.1338561, -1.6168667, -0.4906146},
    { -0.9787684,  1.9161415,  0.0334540},
    {0.0719453, -0.2289914,  1.4052427}
};

const float primaries_list[2][3][2] = {
    { // Rec.2020
        {0.708486, 0.293545},
        {0.190182, 0.775398},
        {0.129252, 0.047142}
    },
    { // Rec.709
        {0.64, 0.33},
        {0.30, 0.60},
        {0.15, 0.06}
    }
};

const float D50_whitepoint[2] = {0.34567, 0.35850};

const float rec2020_xyz[3][3] = invert_f33(xyz_rec2020);

const float MIDDLE_GREY = 0.1845;

const float display_black_target = 0.0152;
const float display_white_target = 100;


//-----------------------------------------------------------------------------
// functions for the sigmoid curve
//-----------------------------------------------------------------------------

float max(float a, float b)
{
    if (a > b) {
        return a;
    } else {
        return b;
    }
}


float generalized_loglogistic_sigmoid(float value,
                                      float magnitude,
                                      float paper_exp,
                                      float film_fog,
                                      float film_power,
                                      float paper_power)
{
    const float clamped_value = max(value, 0.0);
    // The following equation can be derived as a model for film + paper but it has a pole at 0
    // magnitude * powf(1.0f + paper_exp * powf(film_fog + value, -film_power), -paper_power);
    // Rewritten on a stable around zero form:
    const float film_response = pow(film_fog + clamped_value, film_power);
    const float paper_response = magnitude * pow(film_response / (paper_exp + film_response), paper_power);

    // Safety check for very large floats that cause numerical errors
    if (isnan_f(paper_response)) {
        return magnitude;
    } else {
        return paper_response;
    }
}


void calculate_params(float middle_grey_contrast,
                      float contrast_skewness,
                      float display_black_target,
                      float display_white_target,
                      output float film_power,
                      output float white_target,
                      output float black_target,
                      output float film_fog,
                      output float paper_exposure,
                      output float paper_power)
{
    /* Calculate actual skew log logistic parameters to fulfill the following:
     * f(scene_zero) = display_black_target
     * f(scene_grey) = MIDDLE_GREY
     * f(scene_inf)  = display_white_target
     * Slope at scene_grey independent of skewness i.e. only changed by the contrast parameter.
     */

    // Calculate a reference slope for no skew and a normalized display
    const float ref_film_power = middle_grey_contrast;
    const float ref_paper_power = 1.0;
    const float ref_magnitude = 1.0;
    const float ref_film_fog = 0.0;
    const float ref_paper_exposure
        = pow(ref_film_fog + MIDDLE_GREY, ref_film_power) * ((ref_magnitude / MIDDLE_GREY) - 1.0);
    const float delta = 1e-6;
    const float ref_slope
        = (generalized_loglogistic_sigmoid(MIDDLE_GREY + delta, ref_magnitude, ref_paper_exposure, ref_film_fog,
                                           ref_film_power, ref_paper_power)
           - generalized_loglogistic_sigmoid(MIDDLE_GREY - delta, ref_magnitude, ref_paper_exposure, ref_film_fog,
                                             ref_film_power, ref_paper_power))
        / 2.0 / delta;

    // Add skew
    paper_power = pow(5.0, -contrast_skewness);

    // Slope at low film power
    const float temp_film_power = 1.0;
    const float temp_white_target = 0.01 * display_white_target;
    const float temp_white_grey_relation
        = pow(temp_white_target / MIDDLE_GREY, 1.0 / paper_power) - 1.0;
    const float temp_paper_exposure = pow(MIDDLE_GREY, temp_film_power) * temp_white_grey_relation;
    const float temp_slope
        = (generalized_loglogistic_sigmoid(MIDDLE_GREY + delta, temp_white_target, temp_paper_exposure,
                                           ref_film_fog, temp_film_power, paper_power)
           - generalized_loglogistic_sigmoid(MIDDLE_GREY - delta, temp_white_target, temp_paper_exposure,
                                             ref_film_fog, temp_film_power, paper_power))
        / 2.0 / delta;

    // Figure out what film power fulfills the target slope
    // (linear when assuming display_black = 0.0)
    film_power = ref_slope / temp_slope;

    // Calculate the other parameters now that both film and paper power is known
    white_target = 0.01 * display_white_target;
    black_target = 0.01 * display_black_target;
    const float white_grey_relation
        = pow(white_target / MIDDLE_GREY, 1.0 / paper_power) - 1.0;
    const float white_black_relation
        = pow(black_target / white_target, -1.0 / paper_power) - 1.0;

    film_fog = MIDDLE_GREY * pow(white_grey_relation, 1.0 / film_power)
        / (pow(white_black_relation, 1.0 / film_power)
           - pow(white_grey_relation, 1.0 / film_power));
    paper_exposure
        = pow(film_fog + MIDDLE_GREY, film_power) * white_grey_relation;    
}


//-----------------------------------------------------------------------------
// functions for calculating custom primaries
//-----------------------------------------------------------------------------

float determinant(float a, float b, float c, float d)
{
    return a * d - b * c;
}


float intersect_line_segments(float x1, float y1,
                              float x2, float y2,
                              float x3, float y3,
                              float x4, float y4)
{
    const float denominator = determinant(x1 - x2, x3 - x4, y1 - y2, y3 - y4);
    if (denominator == 0.0) {
        return FLT_MAX; // lines don't intersect
    }

    const float t = determinant(x1 - x3, x3 - x4, y1 - y3, y3 - y4) / denominator;
    if (t >= 0.0) {
        return t;
    }
    return FLT_MAX; // intersection is in the wrong direction
}


float find_distance_to_edge(float primaries[3][2],
                            float cos_angle, float sin_angle)
{
    const float x1 = D50_whitepoint[0];
    const float y1 = D50_whitepoint[1];
    const float x2 = x1 + cos_angle;
    const float y2 = y1 + sin_angle;

    float distance_to_edge = FLT_MAX;
    for (int i = 0; i < 3; i = i+1) {
        const int next_i = (i + 1) % 3;
        const float x3 = primaries[i][0];
        const float y3 = primaries[i][1];
        const float x4 = primaries[next_i][0];
        const float y4 = primaries[next_i][1];
        const float distance = intersect_line_segments(x1, y1, x2, y2, x3, y3, x4, y4);
        if (distance < distance_to_edge) {
            distance_to_edge = distance;
        }
    }

    return distance_to_edge;
}


float[3][3] get_matrix_from_xy_coords(float primaries[3][2])
{
    // http://www.brucelindbloom.com/index.html?Eqn_RGB_XYZ_Matrix.html
    float primaries_matrix[3][3] = {{0, 0, 0},
                                    {0, 0, 0},
                                    {0, 0, 0}};
    for (int i = 0; i < 3; i = i+1) {
        // N.B. compared to linked equations, our matrix is transposed
        primaries_matrix[i][0] = primaries[i][0] / primaries[i][1];
        primaries_matrix[i][1] = 1.0;
        primaries_matrix[i][2] = (1.0 - primaries[i][0] - primaries[i][1]) / primaries[i][1];
    }

    float primaries_inverse[3][3] = invert_f33(primaries_matrix);
                                               
    const float XYZ_white[3] = { D50_whitepoint[0] / D50_whitepoint[1], 1.0, (1.0 - D50_whitepoint[0] - D50_whitepoint[1]) / D50_whitepoint[1] };
    float scale[3] = mult_f3_f33(XYZ_white, primaries_inverse);

    float res[3][3];
    for (int i = 0; i < 3; i = i+1) {
        for (int j = 0; j < 3; j = j+1) {
            res[i][j] = scale[i] * primaries_matrix[i][j];
        }
    }
    return res;
}


float[2] rotate_and_scale_primary(float primaries[3][2],
                                  float scaling, float rotation,
                                  int primary_index)
{
    // Generate a custom set of tone mapping primaries by scaling
    // and rotating the primaries of the given profile.
    const float px = primaries[primary_index][0];
    const float py = primaries[primary_index][1];
    const float dx = px - D50_whitepoint[0];
    const float dy = py - D50_whitepoint[1];
    const float angle = atan2(dy, dx) + rotation;
    const float cos_angle = cos(angle);
    const float sin_angle = sin(angle);
    const float distance_to_edge = find_distance_to_edge(primaries,
                                                         cos_angle, sin_angle);
    const float dx_new = scaling * distance_to_edge * cos_angle;
    const float dy_new = scaling * distance_to_edge * sin_angle;
    float new_xy[2] = { dx_new + D50_whitepoint[0],
                        dy_new + D50_whitepoint[1] };
    return new_xy;
}


void calculate_adjusted_primaries(float inset[3], float rotation[3],
                                  int primaries_base_profile, float purity,
                                  output float pipe_to_base[3][3],
                                  output float base_to_rendering[3][3],
                                  output float rendering_to_pipe[3][3])
{
    // Make adjusted primaries for generating the inset matrix
    //
    // References:
    // AgX by Troy Sobotka - https://github.com/sobotka/AgX-S2O3
    // Related discussions on Blender Artists forums -
    // https://blenderartists.org/t/feedback-development-filmic-baby-step-to-a-v2/1361663
    //
    // The idea is to "inset" the work RGB data toward achromatic
    // along spectral lines before per-channel curves. This makes
    // handling of bright, saturated colors much better as the
    // per-channel process desaturates them.
    // The primaries are also rotated to compensate for Abney etc.
    // and achieve a favourable shift towards yellow.

    if (primaries_base_profile == 1) {
        // convert from Rec.2020 to Rec.709
        pipe_to_base = transpose_f33(mult_f33_f33(rec709_xyz, xyz_rec2020));
    } else {
        for (int i = 0; i < 3; i = i+1) {
            for (int j = 0; j < 3; j = j+1) {
                pipe_to_base[i][j] = 0;
            }
            pipe_to_base[i][i] = 1;
        }
    }
    float base_to_pipe[3][3] = invert_f33(pipe_to_base);

    float base_primaries[3][2] = primaries_list[primaries_base_profile];

    // Rotated, scaled primaries are calculated based on the "base profile"
    float custom_primaries[3][2];
    for (int i = 0; i < 3; i = i+1) {
        custom_primaries[i] = rotate_and_scale_primary(base_primaries,
                                                       1.0 - inset[i],
                                                       rotation[i], i);
    }
  
    float custom_to_XYZ[3][3] = get_matrix_from_xy_coords(custom_primaries);

    float out_profile[3][3];
    if (primaries_base_profile == 1) {
        out_profile = transpose_f33(rec709_xyz);
    } else {
        out_profile = transpose_f33(rec2020_xyz);
    }
    base_to_rendering = mult_f33_f33(custom_to_XYZ, out_profile);

    for (int i = 0; i < 3; i = i+1) {
        const float scaling = 1 - purity / 100 * inset[i];
        custom_primaries[i] = rotate_and_scale_primary(base_primaries,
                                                       scaling, rotation[i], i);
    }
    custom_to_XYZ = get_matrix_from_xy_coords(custom_primaries);
    float tmp[3][3] = mult_f33_f33(custom_to_XYZ, out_profile);
    float rendering_to_base[3][3] = invert_f33(tmp);
    rendering_to_pipe = mult_f33_f33(rendering_to_base, base_to_pipe);
}


//-----------------------------------------------------------------------------
// script entry point and parameters definition
//-----------------------------------------------------------------------------

// @ART-param: ["middle_grey_contrast", "$TP_LABCURVE_CONTRAST", 0.7, 3, 1.5, 0.1]
// @ART-param: ["contrast_skewness", "$CTL_SIGMOID_SKEW", -1, 1, -0.2, 0.01]
// @ART-param: ["custom_primaries", "$CTL_SIGMOID_USE_PRIMARIES", true]
// @ART-param: ["primaries_base_profile", "$CTL_SIGMOID_BASE_PRIMARIES", ["Rec. 2020", "Rec. 709 / sRGB"]]
// @ART-param: ["r_inset", "$CTL_SIGMOID_INSET", 0, 50, 10, 0.1, "$TP_CHMIXER_PRIMARY_R"]
// @ART-param: ["r_rotation", "$CTL_SIGMOID_ROTATION", -22.9, 22.9, 2, 0.1, "$TP_CHMIXER_PRIMARY_R"]
// @ART-param: ["g_inset", "$CTL_SIGMOID_INSET", 0, 50, 10, 0.1, "$TP_CHMIXER_PRIMARY_G"]
// @ART-param: ["g_rotation", "$CTL_SIGMOID_ROTATION", -22.9, 22.9, -1, 0.1, "$TP_CHMIXER_PRIMARY_G"]
// @ART-param: ["b_inset", "$CTL_SIGMOID_INSET", 0, 50, 15, 0.1, "$TP_CHMIXER_PRIMARY_B"]
// @ART-param: ["b_rotation", "$CTL_SIGMOID_ROTATION", -22.9, 22.9, -3, 0.1, "$TP_CHMIXER_PRIMARY_B"]
// @ART-param: ["purity", "$CTL_SIGMOID_PURITY", 0, 100, 0, 1]

void ART_main(varying float r,
              varying float g,
              varying float b,
              output varying float rout,
              output varying float gout,
              output varying float bout,
              float middle_grey_contrast,
              float contrast_skewness,
              bool custom_primaries,
              int primaries_base_profile,
              float r_inset, float r_rotation,
              float g_inset, float g_rotation,
              float b_inset, float b_rotation,
              float purity)
{
    float film_power;
    float white_target;
    float black_target;
    float film_fog;
    float paper_exposure;
    float paper_power;

    // compute the sigmoid parameters from the UI controls
    calculate_params(middle_grey_contrast, contrast_skewness,
                     display_black_target, display_white_target, film_power,
                     white_target, black_target, film_fog,
                     paper_exposure, paper_power);

    const float inset[3] = { r_inset / 100, g_inset / 100, b_inset / 100 };
    const float rad = M_PI / 180.0;
    const float rotation[3] =
        { r_rotation * rad, g_rotation * rad, b_rotation * rad };
    float pipe_to_base[3][3];
    float base_to_rendering[3][3];
    float rendering_to_pipe[3][3];
    calculate_adjusted_primaries(inset, rotation, primaries_base_profile, purity,
                                 pipe_to_base, base_to_rendering,
                                 rendering_to_pipe);

    float rgb[3] = {r, g, b};

    if (custom_primaries && primaries_base_profile > 0) {
        rgb = mult_f3_f33(rgb, pipe_to_base);
    }

    for (int i = 0; i < 3; i = i+1) {
        rgb[i] = max(rgb[i], 0);
    }

    if (custom_primaries) {
        rgb = mult_f3_f33(rgb, base_to_rendering);
    }

    for (int i = 0; i < 3; i = i+1) {
        rgb[i] = generalized_loglogistic_sigmoid(rgb[i], white_target,
                                                 paper_exposure, film_fog,
                                                 film_power, paper_power);
    }

    if (custom_primaries) {
        rgb = mult_f3_f33(rgb, rendering_to_pipe);
    }

    rout = rgb[0];
    gout = rgb[1];
    bout = rgb[2];
}

