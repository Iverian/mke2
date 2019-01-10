#ifndef MKE2_SRC_MKE2_CONSTANTS_HPP_
#define MKE2_SRC_MKE2_CONSTANTS_HPP_

namespace cnst {

static constexpr double pi = 3.141592653589;
static constexpr double p = 1e5;          // Па
static constexpr double mu = 1e3;         // Па/м
static constexpr double E = 212e9;        // Па
static constexpr double omega = 2e4 * pi; // Гц
static constexpr double nu = 0.29;        // 1
static constexpr double rho = 10210;      // кг / м^3

static constexpr double L[2]
    = {E * nu / ((1 + nu) * (1 - 2 * nu)), E / (2 * (1 + nu))};
static constexpr double C[36] = {L[0] + 2 * L[1],
                                 L[0],
                                 L[0],
                                 0,
                                 0,
                                 0, //
                                 L[0],
                                 L[0] + 2 * L[1],
                                 L[0],
                                 0,
                                 0,
                                 0, //
                                 L[0],
                                 L[0],
                                 L[0] + 2 * L[1],
                                 0,
                                 0,
                                 0, //
                                 0,
                                 0,
                                 0,
                                 L[1],
                                 0,
                                 0, //
                                 0,
                                 0,
                                 0,
                                 0,
                                 L[1],
                                 0, //
                                 0,
                                 0,
                                 0,
                                 0,
                                 0,
                                 L[1]};
}

#endif // MKE2_SRC_MKE2_CONSTANTS_HPP_