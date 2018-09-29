/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "MaterialLib/PhysicalConstant.h"
#include "TESAssemblyParams.h"

namespace ProcessLib
{
namespace TES
{
inline double fluid_density(const double p, const double T, const double x)
{
    // OGS-5 density model 26

    const double M0 = MaterialLib::PhysicalConstant::MolarMass::N2;
    const double M1 = MaterialLib::PhysicalConstant::MolarMass::Water;

    const double xn = M0 * x / (M0 * x + M1 * (1.0 - x));

    return p / (MaterialLib::PhysicalConstant::IdealGasConstant * T) *
           (M1 * xn + M0 * (1.0 - xn));
    ;
}

template <int i>
double mypow(const double x)
{
    if (i < 0)
    {
        return 1.0 / mypow<-i>(x);
    }

    const double p = mypow<(i >> 1)>(x);
    return (i & 1) ? p * p * x : p * p;
}

template <>
inline double mypow<0>(const double /*x*/)
{
    return 1.0;
}

struct FluidViscosityN2
{
    static double get(double rho, double T)
    {
        const double rho_c = 314;   // [kg/m3]
        const double CVF = 14.058;  // [1e-3 Pa-s]

        const double sigma = 0.36502496e-09;
        const double k = 1.38062e-23;
        const double eps = 138.08483e-23;
        const double c1 = 0.3125;
        const double c2 = 2.0442e-49;

        const double T_star = T * k / eps;
        rho = rho / rho_c;

        double Omega = loop1_term<0>(T_star);
        Omega += loop1_term<1>(T_star);
        Omega += loop1_term<2>(T_star);
        Omega += loop1_term<3>(T_star);
        Omega += loop1_term<4>(T_star);

        Omega = std::exp(Omega);

        // eta in [Pa*s]
        const double eta_0 = c1 * std::sqrt(c2 * T) / (sigma * sigma * Omega);

        double sum = loop2_term<2>(rho);
        sum += loop2_term<3>(rho);
        sum += loop2_term<4>(rho);

        //
        const double eta_r =
            CVF * 1e-6 * (C[0] / (rho - C[1]) + C[0] / C[1] + sum);

        return eta_0 + eta_r;  // [Pa*s]
    }

private:
    template <unsigned i>
    static double loop1_term(double T_star)
    {
        return A[i] * mypow<i>(log(T_star));
    }

    template <unsigned i>
    static double loop2_term(double rho)
    {
        return C[i] * mypow<i - 1>(rho);
    }

    static const double A[5];
    static const double C[5];
};

struct FluidViscosityH2O
{
    static double get(double rho, double T)
    {
        double my, my_0, my_1;
        double H[4];

        T = T / 647.096;
        rho = rho / 322.0;

        H[0] = 1.67752;
        H[1] = 2.20462;
        H[2] = 0.6366564;
        H[3] = -0.241605;

        double h[6][7] = {{0.0}};
        h[0][0] = 0.520094000;
        h[1][0] = 0.085089500;
        h[2][0] = -1.083740000;
        h[3][0] = -0.289555000;
        h[0][1] = 0.222531000;
        h[1][1] = 0.999115000;
        h[2][1] = 1.887970000;
        h[3][1] = 1.266130000;
        h[5][1] = 0.120573000;
        h[0][2] = -0.281378000;
        h[1][2] = -0.906851000;
        h[2][2] = -0.772479000;
        h[3][2] = -0.489837000;
        h[4][2] = -0.257040000;
        h[0][3] = 0.161913000;
        h[1][3] = 0.257399000;
        h[0][4] = -0.032537200;
        h[3][4] = 0.069845200;
        h[4][5] = 0.008721020;
        h[3][6] = -0.004356730;
        h[5][6] = -0.000593264;

        double sum1 = H[0] / mypow<0>(T);
        sum1 += H[1] / mypow<1>(T);
        sum1 += H[2] / mypow<2>(T);
        sum1 += H[3] / mypow<3>(T);

        my_0 = 100 * std::sqrt(T) / sum1;

        double sum2 = inner_loop<0>(rho, T, h);
        sum2 += inner_loop<1>(rho, T, h);
        sum2 += inner_loop<2>(rho, T, h);
        sum2 += inner_loop<3>(rho, T, h);
        sum2 += inner_loop<4>(rho, T, h);
        sum2 += inner_loop<5>(rho, T, h);

        my_1 = std::exp(rho * sum2);

        my = (my_0 * my_1) / 1e6;
        return my;
    }

private:
    template <int i>
    static double inner_loop(const double rho,
                             const double T,
                             const double (&h)[6][7])
    {
        const double base = rho - 1.0;

        double sum3 = h[i][0] * mypow<0>(base);
        sum3 += h[i][1] * mypow<1>(base);
        sum3 += h[i][2] * mypow<2>(base);
        sum3 += h[i][3] * mypow<3>(base);
        sum3 += h[i][4] * mypow<4>(base);
        sum3 += h[i][5] * mypow<5>(base);
        sum3 += h[i][6] * mypow<6>(base);

        return mypow<i>(1 / T - 1) * sum3;
    }
};

inline double fluid_viscosity(const double p, const double T, const double x)
{
    // OGS 5 viscosity model 26

    const double M0 = MaterialLib::PhysicalConstant::MolarMass::N2;
    const double M1 = MaterialLib::PhysicalConstant::MolarMass::Water;
    const double R = MaterialLib::PhysicalConstant::IdealGasConstant;

    // reactive component
    const double x0 =
        M0 * x / (M0 * x + M1 * (1.0 - x));  // mass in mole fraction
    const double V0 = FluidViscosityH2O::get(M1 * p / (R * T), T);
    // inert component
    const double x1 = 1.0 - x0;
    const double V1 = FluidViscosityN2::get(M0 * p / (R * T), T);

    const double M0_over_M1(M1 / M0);  // reactive over inert
    const double V0_over_V1(V0 / V1);

    const double phi_12 =
        mypow<2>(1.0 +
                 std::sqrt(V0_over_V1) * std::pow(1.0 / M0_over_M1, 0.25)) /
        std::sqrt(8.0 * (1.0 + M0_over_M1));

    const double phi_21 = phi_12 * M0_over_M1 / V0_over_V1;

    return V0 * x0 / (x0 + x1 * phi_12) + V1 * x1 / (x1 + x0 * phi_21);
}

struct FluidHeatConductivityN2
{
    static double get(double rho, double T)
    {
        const double X1 = 0.95185202;
        const double X2 = 1.0205422;

        const double rho_c = 314;  // [kg/m3]
        const double M = 28.013;
        const double k = 1.38062e-23;
        const double eps = 138.08483e-23;
        const double N_A = 6.02213E26;
        const double R = 8.31434;
        // const double R = MaterialLib::PhysicalConstant::IdealGasConstant;
        const double CCF = 4.173;  // mW/m/K

        const double c1 = 0.3125;
        const double c2 = 2.0442e-49;
        const double sigma = 0.36502496e-09;

        rho /= rho_c;

        // dilute heat conductivity
        const double sum1 = loop1_term<0>(T) + loop1_term<1>(T) +
                            loop1_term<2>(T) + loop1_term<3>(T) +
                            loop1_term<4>(T) + loop1_term<5>(T) +
                            loop1_term<6>(T);
        const double temp(std::exp((f[8] / T)) - 1);
        const double c_v0 =
            R *
            (sum1 + ((f[7] * (f[8] / T) * (f[8] / T) * (std::exp((f[8] / T)))) /
                         (temp * temp) -
                     1));

        double cvint;
        cvint = c_v0 * 1000 / N_A;

        // dilute gas viscosity
        const double log_T_star = std::log(T * k / eps);

        const double Omega =
            std::exp(loop2_term<0>(log_T_star) + loop2_term<1>(log_T_star) +
                     loop2_term<2>(log_T_star) + loop2_term<3>(log_T_star) +
                     loop2_term<4>(log_T_star));

        // eta in [Pa*s]
        const double eta_0 =
            1e6 * (c1 * std::sqrt(c2 * T) / (sigma * sigma * Omega));

        const double F = eta_0 * k * N_A / (M * 1000);

        const double lambda_tr = 2.5 * (1.5 - X1);
        const double lambda_in = X2 * (cvint / k + X1);

        const double lambda_0 = F * (lambda_tr + lambda_in);

        const double sum2 = loop3_term<0>(rho) + loop3_term<1>(rho) +
                            loop3_term<2>(rho) + loop3_term<3>(rho);
        const double lambda_r = sum2 * CCF;

        return (lambda_0 + lambda_r) / 1000;  // lambda in [W/m/K]
    }

private:
    template <int i>
    static double loop1_term(const double T)
    {
        return f[i] * mypow<i - 3>(T);
    }

    template <int i>
    static double loop2_term(const double log_T_star)
    {
        return A[i] * mypow<i>(log_T_star);
    }

    template <int i>
    static double loop3_term(const double rho)
    {
        return C[i] * mypow<i + 1>(rho);
    }

    const static double A[5];
    const static double f[9];
    const static double C[4];
};

struct FluidHeatConductivityH2O
{
    static double get(double rho, double T)
    {
        double S, Q;
        double b[3], B[2], d[4], C[6];

        T /= 647.096;
        rho /= 317.11;

        b[0] = -0.397070;
        b[1] = 0.400302;
        b[2] = 1.060000;

        B[0] = -0.171587;
        B[1] = 2.392190;

        d[0] = 0.0701309;
        d[1] = 0.0118520;
        d[2] = 0.00169937;
        d[3] = -1.0200;

        C[0] = 0.642857;
        C[1] = -4.11717;
        C[2] = -6.17937;
        C[3] = 0.00308976;
        C[4] = 0.0822994;
        C[5] = 10.0932;

        const double sum1 = loop_term<0>(T) + loop_term<1>(T) +
                            loop_term<2>(T) + loop_term<3>(T);

        const double lambda_0 = std::sqrt(T) * sum1;
        const double lambda_1 =
            b[0] + b[1] * rho +
            b[2] * std::exp(B[0] * (rho + B[1]) * (rho + B[1]));

        const double dT = fabs(T - 1) + C[3];
        const double dT_pow_3_5 = std::pow(dT, 3. / 5.);
        Q = 2 + (C[4] / dT_pow_3_5);

        if (T >= 1)
            S = 1 / dT;
        else
            S = C[5] / dT_pow_3_5;

        const double rho_pow_9_5 = std::pow(rho, 9. / 5.);
        const double rho_pow_Q = std::pow(rho, Q);
        const double T_pow_3_2 = T * std::sqrt(T);
        const double lambda_2 =
            (d[0] / mypow<10>(T) + d[1]) * rho_pow_9_5 *
                std::exp(C[0] * (1 - rho * rho_pow_9_5)) +
            d[2] * S * rho_pow_Q *
                std::exp((Q / (1. + Q)) * (1 - rho * rho_pow_Q)) +
            d[3] * std::exp(C[1] * T_pow_3_2 + C[2] / mypow<5>(rho));

        return lambda_0 + lambda_1 + lambda_2;  // lambda in [W/m/K]
    }

private:
    template <unsigned i>
    static double loop_term(const double T)
    {
        return a[i] * mypow<i>(T);
    }

    static const double a[4];
};

inline double fluid_heat_conductivity(const double p,
                                      const double T,
                                      const double x)
{
    // OGS 5 fluid heat conductivity model 11

    const double M0 = MaterialLib::PhysicalConstant::MolarMass::N2;
    const double M1 = MaterialLib::PhysicalConstant::MolarMass::Water;
    const double R = MaterialLib::PhysicalConstant::IdealGasConstant;

    // TODO [CL] max() is redundant if the fraction is guaranteed to be between
    // 0 and 1.
    // reactive component
    const double x0 = std::max(M0 * x / (M0 * x + M1 * (1.0 - x)),
                               0.);  // convert mass to mole fraction
    const double k0 = FluidHeatConductivityH2O::get(M1 * p / (R * T), T);
    // inert component
    const double x1 = 1.0 - x0;
    const double k1 = FluidHeatConductivityN2::get(M0 * p / (R * T), T);

    const double M1_over_M2 = M1 / M0;  // reactive over inert
    const double V1_over_V2 = FluidViscosityH2O::get(M1 * p / (R * T), T) /
                              FluidViscosityN2::get(M0 * p / (R * T), T);
    const double L1_over_L2 = V1_over_V2 / M1_over_M2;

    const double M12_pow_mquarter = std::pow(M1_over_M2, -0.25);
    const double phi_12 = (1.0 + std::sqrt(L1_over_L2) * M12_pow_mquarter) *
                          (1.0 + std::sqrt(V1_over_V2) * M12_pow_mquarter) /
                          std::sqrt(8.0 * (1.0 + M1_over_M2));
    const double phi_21 = phi_12 * M1_over_M2 / V1_over_V2;

    return k0 * x0 / (x0 + x1 * phi_12) + k1 * x1 / (x1 + x0 * phi_21);
}

}  // TES
}  // ProcessLib
