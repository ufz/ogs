/**
 *  \file
 *  \copyright
 *   Copyright (c) 2012-2024, OpenGeoSys Communty (http://www.opengeosys.org)
 *              Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "DimensionlessGibbsFreeEnergyRegion2.h"

#include <array>
#include <cmath>

#include "BaseLib/Error.h"

namespace MaterialLib
{
namespace Fluid
{
namespace DimensionlessGibbsFreeEnergyRegion2
{
static constexpr std::array<int, 9> J0 = {0, 1, -5, -4, -3, -2, -1, 2, 3};
static constexpr std::array<double, 9> n0 = {
    -0.96927686500217e1, 0.10086655968018e2, -0.56087911283020e-2,
    0.71452738081455e-1, -0.40710498223928,  0.14240819171444e1,
    -0.43839511319450e1, -0.28408632460772,  0.21268463753307e-1};

static constexpr std::array<double, 43> n = {
    -0.17731742473213e-2,  -0.17834862292358e-1,  -0.45996013696365e-1,
    -0.57581259083432e-1,  -0.50325278727930e-1,  -0.33032641670203e-4,
    -0.18948987516315e-3,  -0.39392777243355e-2,  -0.43797295650573e-1,
    -0.26674547914087e-4,  0.20481737692309e-7,   0.43870667284435e-6,
    -0.32277677238570e-4,  -0.15033924542148e-2,  -0.40668253562649e-1,
    -0.78847309559367e-9,  0.12790717852285e-7,   0.48225372718507e-6,
    0.22922076337661e-5,   -0.16714766451061e-10, -0.21171472321355e-2,
    -0.23895741934104e2,   -0.59059564324270e-17, -0.12621808899101e-5,
    -0.38946842435739e-1,  0.11256211360459e-10,  -0.82311340897998e1,
    0.19809712802088e-7,   0.10406965210174e-18,  -0.10234747095929e-12,
    -0.10018179379511e-8,  -0.80882908646985e-10, 0.10693031879409,
    -0.33662250574171,     0.89185845355421e-24,  0.30629316876232e-12,
    -0.42002467698208e-5,  -0.59056029685639e-25, 0.37826947613457e-5,
    -0.12768608934681e-14, 0.73087610595061e-28,  0.55414715350778e-16,
    -0.94369707241210e-6};

static constexpr std::array<int, 43> I = {
    1,  1,  1,  1,  1,  2,  2,  2,  2,  2,  3,  3,  3, 3,  3,
    4,  4,  4,  5,  6,  6,  6,  7,  7,  7,  8,  8,  9, 10, 10,
    10, 16, 16, 18, 20, 20, 20, 21, 22, 23, 24, 24, 24};

static constexpr std::array<int, 43> J = {
    0,  1,  2,  3,  6,  1,  2,  4,  7,  36, 0,  1,  3,  6, 35,
    1,  2,  3,  7,  3,  16, 35, 0,  11, 25, 8,  36, 13, 4, 10,
    14, 29, 50, 57, 20, 35, 48, 21, 53, 39, 26, 40, 58};

double computeGamma0(const double tau, const double pi)
{
    if ((pi < 0.) || (pi == 0.))
    {
        OGS_FATAL(
            "The dimensionless Gibbs free energy in IAPWS-IF97 region2 can not "
            "be calculated from a non-positive pressure.");
    }

    double val = std::log(pi);
    for (int i = 0; i < 9; i++)
    {
        val += n0[i] * std::pow(tau, J0[i]);
    }

    return val;
}

double getGamma(const double tau, const double pi)
{
    double val = computeGamma0(tau, pi);
    for (int i = 0; i < 43; i++)
    {
        val += n[i] * std::pow(pi, I[i]) * std::pow(tau - 0.5, J[i]);
    }

    return val;
}

double getdGammadTau(const double tau, const double pi)
{
    double val1 = 0.;
    double val2 = 0.;
    for (int i = 0; i < 9; i++)
    {
        val1 += n0[i] * J0[i] * std::pow(tau, J0[i] - 1);
    }

    for (int i = 0; i < 43; i++)
    {
        val2 +=
            n[i] * J[i] * std::pow(pi, I[i]) * std::pow(tau - 0.5, J[i] - 1);
    }

    return val1 + val2;
}

double getdGammadPi(const double tau, const double pi)
{
    if ((pi < 0.) || (pi == 0.))
    {
        OGS_FATAL(
            "The dimensionless Gibbs free energy in IAPWS-IF97 region2 can not "
            "be calculated from a non-positive pressure.");
    }

    double val = 1 / pi;
    for (int i = 0; i < 43; i++)
    {
        val += n[i] * I[i] * std::pow(pi, I[i] - 1) * std::pow(tau - 0.5, J[i]);
    }

    return val;
}

}  // namespace DimensionlessGibbsFreeEnergyRegion2
}  // namespace Fluid
}  // namespace MaterialLib
