#ifndef MKE2_INCLUDE_CONSTANTS_HPP_
#define MKE2_INCLUDE_CONSTANTS_HPP_

#include <common_types.hpp>

// #define MKE2_DENSE_SOLVE

namespace cnst {

static constexpr Index scale = 10;
static constexpr Index max_iter = 100000;

static constexpr Value xdim = 200.;
static constexpr Value ydim = 40.;
static constexpr Value zdim = 40.;
static constexpr Value init = 0.;

static constexpr Value pi = 3.141592653589; // 1
static constexpr Value p = 1e5;             // Па
static constexpr Value mu = 1e3;            // Па/м
static constexpr Value E = 212e9;           // Па
static constexpr Value omega = 2e4 * pi;    // Гц
static constexpr Value nu = 0.29;           // 1
static constexpr Value rho = 10210;         // кг / м^3

static constexpr Value L[2]
    = {E * nu / ((1 + nu) * (1 - 2 * nu)), E / (2 * (1 + nu))};
static constexpr Value C[36] = {
    L[0] + 2 * L[1], L[0], L[0], 0, 0, 0, //
    L[0], L[0] + 2 * L[1], L[0], 0, 0, 0, //
    L[0], L[0], L[0] + 2 * L[1], 0, 0, 0, //
    0, 0, 0, L[1], 0, 0,                  //
    0, 0, 0, 0, L[1], 0,                  //
    0, 0, 0, 0, 0, L[1]                   //
};
}

#endif // MKE2_INCLUDE_CONSTANTS_HPP_