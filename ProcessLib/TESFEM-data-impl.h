/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * The code of this file is used to decouple the evaluation of matrix elements from the rest of OGS6,
 * not all of OGS6 has to be recompiled every time a small change is done.
 */

#ifndef PROCESS_LIB_TESFEM_DATA_IMPL_H_
#define PROCESS_LIB_TESFEM_DATA_IMPL_H_

#include <iostream>
#include <cstdio>

#include <typeinfo>

#include "logog/include/logog.hpp"

#include "TESFEM-data-fwd.h"

#include "NumLib/Function/Interpolation.h"

#include "TESFEMReactionAdaptor.h"

namespace
{
const double GAS_CONST = 8.3144621;

enum class MatOutType { OGS5, PYTHON };

const MatOutType MATRIX_OUTPUT_FORMAT = MatOutType::PYTHON;
}


static double fluid_density(const double p, const double T, const double x)
{
	// OGS-5 density model 26

	const double M0 = ProcessLib::TES::M_N2;
	const double M1 = ProcessLib::TES::M_H2O;

	const double xn = M0*x/(M0*x + M1*(1.0-x));

	return p / (GAS_CONST * T) * (M1*xn + M0*(1.0-xn));;
}



template<int i>
inline double mypow(const double x)
{
    if (i<0) {
        return 1.0/mypow<-i>(x);
    } else {
        const double p = mypow<(i>>1)>(x);
        return (i&1) ? p*p*x : p*p;
    }
}

template<>
inline double mypow<0>(const double /*x*/)
{
    return 1.0;
}

template<>
inline double mypow<1>(const double x)
{
    return x;
}

template<>
inline double mypow<2>(const double x)
{
    return x*x;
}



struct FluidViscosityN2
{
	static double get(double rho, double T)
	{
        const double rho_c = 314;             // [kg/m3]
        const double CVF = 14.058;            // [1e-3 Pa-s]

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

        Omega = std::exp (Omega);

        //eta in [Pa*s]
        const double eta_0 = c1 * std::sqrt(c2 * T) / (sigma * sigma * Omega);

        double sum = loop2_term<2>(rho);
        sum += loop2_term<3>(rho);
        sum += loop2_term<4>(rho);

        //
        const double eta_r = CVF * 1e-6 * (C[0] / (rho - C[1]) + C[0] / C[1] + sum);

        return eta_0 + eta_r;               // [Pa*s]
	}

private:
	template<unsigned i>
	static double loop1_term(double T_star)
	{
		return A[i] * mypow<i>(log(T_star));
	}

	template<unsigned i>
	static double loop2_term(double rho)
	{
		return C[i] * mypow<i-1>(rho);
	}

	static constexpr double A[5] = { 0.46649, -0.57015, 0.19164, -0.03708, 0.00241 };
	static constexpr double C[5] = { -20.09997, 3.4376416, -1.4470051, -0.027766561, -0.21662362 };
};


struct FluidViscosityH2O
{
	static double get(double rho, double T)
	{
		double my,my_0,my_1;
		double H[4];

		T = T / 647.096;
		rho = rho / 322.0;

		H[0] = 1.67752;
		H[1] = 2.20462;
		H[2] = 0.6366564;
		H[3] = -0.241605;

		double h[6][7] = { 0.0 };
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
	static double inner_loop(const double rho, const double T, const double (&h)[6][7])
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

static double fluid_viscosity(const double p, const double T, const double x)
{
	// OGS 5 viscosity model 26

	const double M0 = ProcessLib::TES::M_N2;
	const double M1 = ProcessLib::TES::M_H2O;

	//reactive component
	const double x0 = M0*x/(M0*x + M1*(1.0-x)); //mass in mole fraction
	const double V0 = FluidViscosityH2O::get(M1*p/(GAS_CONST*T), T);
	//inert component
	const double x1 = 1.0 - x0;
	const double V1 = FluidViscosityN2::get(M0*p/(GAS_CONST*T), T);

	const double M0_over_M1 (M1/M0); //reactive over inert
	const double V0_over_V1 (V0/V1);

	const double phi_12 = mypow<2>(1.0 + std::sqrt(V0_over_V1) * std::pow(1.0/M0_over_M1, 0.25))
						  / std::sqrt(8.0*(1.0+M0_over_M1));

	const double phi_21 = phi_12 * M0_over_M1 / V0_over_V1;

	return V0*x0 / (x0 + x1 * phi_12)
			+ V1*x1 / (x1 + x0 * phi_21);
}


struct FluidHeatConductivityN2
{
	static double get(double rho, double T)
	{
		const double X1 = 0.95185202;
		const double X2 = 1.0205422;

		const double rho_c = 314;             // [kg/m3]
		const double M = 28.013;
		const double k = 1.38062e-23;
		const double eps = 138.08483e-23;
		const double N_A = 6.02213E26;
		const double R = 8.31434;
		// const double R = GAS_CONST;
		const double CCF = 4.173;             //mW/m/K

		const double c1 = 0.3125;
		const double c2 = 2.0442e-49;
		const double sigma = 0.36502496e-09;

		rho /= rho_c;

		// dilute heat conductivity
		const double sum1 = loop1_term<0>(T) + loop1_term<1>(T) + loop1_term<2>(T)
							+ loop1_term<3>(T) + loop1_term<4>(T) + loop1_term<5>(T)
							+ loop1_term<6>(T);
		const double temp (std::exp ((f[8] / T)) - 1);
		const double c_v0
				= R * (sum1 + ((f[7] * (f[8] / T) * (f[8] / T) * (std::exp((f[8] / T))))
				                / (temp * temp) - 1));

		double cvint;
		cvint = c_v0 * 1000 / N_A;

		// dilute gas viscosity
		const double log_T_star = std::log(T * k / eps);

		const double Omega
				= std::exp(
					  loop2_term<0>(log_T_star) + loop2_term<1>(log_T_star) + loop2_term<2>(log_T_star) +
					  loop2_term<3>(log_T_star) + loop2_term<4>(log_T_star)
					  );

		//eta in [Pa*s]
		const double eta_0 = 1e6 * (c1 * std::sqrt(c2 * T) / (sigma * sigma * Omega));

		const double F = eta_0 * k * N_A / (M * 1000);

		const double lambda_tr = 2.5 * (1.5 - X1);
		const double lambda_in = X2 * (cvint / k + X1);

		const double lambda_0 = F * (lambda_tr + lambda_in);

		const double sum2 = loop3_term<0>(rho) + loop3_term<1>(rho)
							+ loop3_term<2>(rho) + loop3_term<3>(rho);
		const double lambda_r = sum2 * CCF;

		return (lambda_0 + lambda_r) / 1000;   //lambda in [W/m/K]
	}

private:
	template<int i>
	static double loop1_term(const double T)
	{
		return f[i] * mypow<i-3>(T);
	}

	template<int i>
	static double loop2_term(const double log_T_star)
	{
		return A[i] * mypow<i>(log_T_star);
	}

	template<int i>
	static double loop3_term(const double rho)
	{
		return C[i] * mypow<i+1>(rho);
	}

	constexpr static double A[5] = {
		0.46649, -0.57015, 0.19164, -0.03708, 0.00241
	};

	constexpr static double f[9] = {
		-0.837079888737e3,   0.37914711487e2,  -0.601737844275,
		 0.350418363823e1,  -0.874955653028e-5, 0.148968607239e-7,
		-0.256370354277e-11, 0.100773735767e1,  0.335340610e4
	};

	constexpr static double C[4] = {
		3.3373542, 0.37098251, 0.89913456, 0.16972505
	};
};


struct FluidHeatConductivityH2O
{
	static double get(double rho, double T)
	{
		double S, Q;
		double b[3], B[2], d[4], C[6];

		T   /= 647.096;
		rho /= 317.11;

		b[0] = -0.397070;
		b[1] =  0.400302;
		b[2] =  1.060000;

		B[0] = -0.171587;
		B[1] =  2.392190;

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

		const double sum1 = loop_term<0>(T) + loop_term<1>(T)
							+ loop_term<2>(T) + loop_term<3>(T);

		const double lambda_0 = std::sqrt(T) * sum1;
		const double lambda_1 = b[0] + b[1] * rho + b[2] * std::exp(B[0] * (rho + B[1]) * (rho + B[1]));

		const double dT = fabs(T - 1) + C[3];
		const double dT_pow_3_5 = std::pow(dT, 3./5.);
		Q = 2 + (C[4] / dT_pow_3_5);

		if (T >= 1)
			S = 1 / dT;
		else
			S = C[5] / dT_pow_3_5;

		const double rho_pow_9_5 = std::pow(rho, 9./5.);
		const double rho_pow_Q   = std::pow(rho, Q);
		const double T_pow_3_2   = T * std::sqrt(T);
		const double lambda_2 =
				(d[0] /
				 mypow<10>(T) + d[1]) * rho_pow_9_5 * std::exp(C[0] * (1 - rho * rho_pow_9_5))
				+ d[2]* S * rho_pow_Q * std::exp((Q / (1. + Q)) * (1 - rho * rho_pow_Q))
				+ d[3] * std::exp(C[1] * T_pow_3_2 + C[2] / mypow<5>(rho));

		return lambda_0 + lambda_1 + lambda_2; // lambda in [W/m/K]
	}

private:
	template<unsigned i>
	static double loop_term(const double T)
	{
		return a[i] * mypow<i>(T);
	}

	static constexpr double a[4] = { 0.0102811, 0.0299621, 0.0156146, -0.00422464 };
};


static double fluid_heat_conductivity(const double p, const double T, const double x)
{
	// OGS 5 fluid heat conductivity model 11

	const double M0 = ProcessLib::TES::M_N2;
	const double M1 = ProcessLib::TES::M_H2O;

	// TODO [CL] max() is redundant if the fraction is guaranteed to be between 0 and 1.
	//reactive component
	const double x0 = std::max(M0*x/(M0*x + M1*(1.0-x)), 0.); // convert mass to mole fraction
	const double k0 = FluidHeatConductivityH2O::get(M1*p/(GAS_CONST*T), T);
	//inert component
	const double x1 = 1.0 - x0;
	const double k1 = FluidHeatConductivityN2::get(M0*p/(GAS_CONST*T), T);

	const double M1_over_M2 = M1/M0; //reactive over inert
	const double V1_over_V2 = FluidViscosityH2O::get(M1*p/(GAS_CONST*T), T)
							/ FluidViscosityN2::get(M0*p/(GAS_CONST*T), T);
	const double L1_over_L2 = V1_over_V2 / M1_over_M2;

	const double M12_pow_mquarter = std::pow(M1_over_M2, -0.25);
	const double phi_12 =   (1.0 + std::sqrt(L1_over_L2) * M12_pow_mquarter)
						  * (1.0 + std::sqrt(V1_over_V2) * M12_pow_mquarter)
						  / std::sqrt(8.0 * (1.0 + M1_over_M2));
	const double phi_21 = phi_12 * M1_over_M2 / V1_over_V2;

	return k0*x0 / (x0+x1*phi_12) + k1*x1 / (x1+x0*phi_21);
}


namespace ProcessLib
{

namespace TES
{

template<typename Traits>
Eigen::Matrix3d
LADataNoTpl<Traits>::
getMassCoeffMatrix(const unsigned int_pt)
{
	// TODO: Dalton's law property
	const double dxn_dxm = Ads::Adsorption::d_molar_fraction(
							   _vapour_mass_fraction, _AP->_M_react, _AP->_M_inert);

	const double M_pp = _AP->_poro/_p * _rho_GR;
	const double M_pT = -_AP->_poro/_T *  _rho_GR;
	const double M_px = (_AP->_M_react-_AP->_M_inert) * _p
						/ (GAS_CONST * _T) * dxn_dxm * _AP->_poro;

	const double M_Tp = -_AP->_poro;
	const double M_TT =
			_AP->_poro * _rho_GR * _AP->_cpG // TODO: vapour heat capacity
			+ (1.0-_AP->_poro) * _solid_density[int_pt] * _AP->_cpS; // TODO: adsorbate heat capacity
	const double M_Tx = 0.0;

	const double M_xp = 0.0;
	const double M_xT = 0.0;
	const double M_xx = _AP->_poro * _rho_GR;


	Eigen::Matrix3d M;
	M << M_pp, M_pT, M_px,
		 M_Tp, M_TT, M_Tx,
		 M_xp, M_xT, M_xx;

	return M;
}


template<typename Traits>
typename Traits::LaplaceMatrix
LADataNoTpl<Traits>::
getLaplaceCoeffMatrix(const unsigned /*int_pt*/, const unsigned dim)
{
	const double eta_GR = fluid_viscosity(_p, _T, _vapour_mass_fraction);

	const double lambda_F = fluid_heat_conductivity(_p, _T, _vapour_mass_fraction);
	const double lambda_S = _AP->_solid_heat_cond;

	using Mat = typename Traits::MatrixDimDim;

	typename Traits::LaplaceMatrix L
			= Traits::LaplaceMatrix::Zero(dim*NODAL_DOF, dim*NODAL_DOF);

	// TODO: k_rel
	// L_pp
	Traits::blockDimDim(L,     0,     0, dim, dim)
			= Traits::blockDimDim(_AP->_solid_perm_tensor, 0,0,dim,dim) * _rho_GR / eta_GR;

	// TODO: add zeolite part
	// L_TT
	Traits::blockDimDim(L,   dim,   dim, dim, dim)
			= Mat::Identity(dim, dim)
			  * ( _AP->_poro * lambda_F + (1.0 - _AP->_poro) * lambda_S);

	// L_xx
	Traits::blockDimDim(L, 2*dim, 2*dim, dim, dim)
			= Mat::Identity(dim, dim)
			  * (_AP->_tortuosity * _AP->_poro * _rho_GR
				 * _AP->_diffusion_coefficient_component
				 );

	return L;
}


template<typename Traits>
Eigen::Matrix3d
LADataNoTpl<Traits>::
getAdvectionCoeffMatrix(const unsigned /*int_pt*/)
{
	const double A_pp = 0.0;
	const double A_pT = 0.0;

	const double A_px = 0.0;

	const double A_Tp = 0.0;

	const double A_TT = _rho_GR * _AP->_cpG; // porosity?
	const double A_Tx = 0.0;

	const double A_xp = 0.0;
	const double A_xT = 0.0;
	const double A_xx = _rho_GR; // porosity?


	Eigen::Matrix3d A;
	A << A_pp, A_pT, A_px,
		 A_Tp, A_TT, A_Tx,
		 A_xp, A_xT, A_xx;

	return A;
}


template<typename Traits>
Eigen::Matrix3d
LADataNoTpl<Traits>::
getContentCoeffMatrix(const unsigned /*int_pt*/)
{
	const double C_pp = 0.0;
	const double C_pT = 0.0;

	const double C_px = 0.0;

	const double C_Tp = 0.0;

	const double C_TT = 0.0;
	const double C_Tx = 0.0;

	const double C_xp = 0.0;
	const double C_xT = 0.0;
	const double C_xx = (_AP->_poro - 1.0) * _qR;

	Eigen::Matrix3d C;
	C << C_pp, C_pT, C_px,
		 C_Tp, C_TT, C_Tx,
		 C_xp, C_xT, C_xx;

	return C;
}


template<typename Traits>
Eigen::Vector3d
LADataNoTpl<Traits>::
getRHSCoeffVector(const unsigned int_pt)
{
	const double reaction_enthalpy = _AP->_reaction_system->get_enthalpy(_p_V, _T, _AP->_M_react);

	const double rhs_p = (_AP->_poro - 1.0) * _qR; // TODO [CL] body force term

	const double rhs_T = _rho_GR * _AP->_poro * _AP->_fluid_specific_heat_source
						 + (1.0 - _AP->_poro) * _qR * reaction_enthalpy
						 + _solid_density[int_pt] * (1.0 - _AP->_poro) * _AP->_solid_specific_heat_source;
						 // TODO [CL] momentum production term

	const double rhs_x = (_AP->_poro - 1.0) * _qR; // TODO [CL] what if x < 0.0


	Eigen::Vector3d rhs;
	rhs << rhs_p,
		 rhs_T,
		 rhs_x;

	return rhs;
}


template<typename Traits>
void
LADataNoTpl<Traits>::
initReaction(const unsigned int_pt)
{
    _reaction_adaptor->initReaction(int_pt);
}


template<typename Traits>
void
LADataNoTpl<Traits>::
preEachAssembleIntegrationPoint(
        const unsigned int_pt,
        const std::vector<double> &localX,
        typename Traits::ShapeMatrices::ShapeType const& smN,
        typename Traits::ShapeMatrices::DxShapeType const& /*smDNdx*/,
        typename Traits::ShapeMatrices::JacobianType const& /*smJ*/,
        const double /*smDetJ*/)
{
#ifndef NDEBUG
    // fill local data with garbage to aid in debugging
    _p = _T   = _vapour_mass_fraction
       = _p_V = _rho_GR
       = _qR
       = std::numeric_limits<double>::quiet_NaN();
#endif

    std::array<double*, NODAL_DOF> int_pt_val = { &_p, &_T, &_vapour_mass_fraction };

    NumLib::shapeFunctionInterpolate(localX, smN, int_pt_val);

    // pre-compute certain properties
    _p_V = _p * Ads::Adsorption::get_molar_fraction(_vapour_mass_fraction, _AP->_M_react, _AP->_M_inert);

    initReaction(int_pt);

    assert(_p > 0.0);
    assert(_T > 0.0);
    assert(0.0 <= _vapour_mass_fraction && _vapour_mass_fraction <= 1.0);

    _rho_GR = fluid_density(_p, _T, _vapour_mass_fraction);
}


template<typename Mat>
void
ogs5OutMat(const Mat& mat)
{
    for (unsigned r=0; r<mat.rows(); ++r)
    {
        switch (MATRIX_OUTPUT_FORMAT)
        {
        case MatOutType::OGS5:
            if (r!=0) std::printf("\n");
            std::printf("|");
            break;
        case MatOutType::PYTHON:
            if (r!=0) std::printf(",\n");
            std::printf("[");
            break;
        }

        for (unsigned c=0; c<mat.cols(); ++c)
        {
            switch (MATRIX_OUTPUT_FORMAT)
            {
            case MatOutType::OGS5:
                std::printf(" %.16e", mat(r, c));
                break;
            case MatOutType::PYTHON:
                if (c!=0) std::printf(",");
                std::printf(" %23.16g", mat(r, c));
                break;
            }

        }

        switch (MATRIX_OUTPUT_FORMAT)
        {
        case MatOutType::OGS5:
            std::printf(" | ");
            break;
        case MatOutType::PYTHON:
            std::printf(" ]");
            break;
        }
    }
    std::printf("\n");
}


template<typename Vec>
void
ogs5OutVec(const Vec& vec)
{
    for (unsigned r=0; r<vec.size(); ++r)
    {
        switch (MATRIX_OUTPUT_FORMAT)
        {
        case MatOutType::OGS5:
            if (r!=0) std::printf("\n");
            std::printf("| %.16e | ", vec[r]);
            break;
        case MatOutType::PYTHON:
            if (r!=0) std::printf(",\n");
            std::printf("[ %23.16g ]", vec[r]);
            break;
        }
    }
    std::printf("\n");
}


template<typename Traits>
std::vector<double> const&
LADataNoTpl<Traits>::
getIntegrationPointValues(SecondaryVariables var, std::vector<double>& cache) const
{
    switch (var)
    {
    case SecondaryVariables::REACTION_RATE:
        return _reaction_rate;
    case SecondaryVariables::SOLID_DENSITY:
        return _solid_density;
    case SecondaryVariables::VELOCITY_X:
        return _velocity[0];
    case SecondaryVariables::VELOCITY_Y:
        assert(_velocity.size() >= 2);
        return _velocity[1];
    case SecondaryVariables::VELOCITY_Z:
        assert(_velocity.size() >= 3);
        return _velocity[2];

    case SecondaryVariables::LOADING:
    {
        auto& Cs = cache;
        Cs.clear();
        Cs.reserve(_solid_density.size());

        for (auto rho_SR : _solid_density) {
            Cs.push_back(rho_SR / _AP->_rho_SR_dry - 1.0);
        }

        return Cs;
    }

    case SecondaryVariables::VAPOUR_PARTIAL_PRESSURE:
    case SecondaryVariables::RELATIVE_HUMIDITY:
    case SecondaryVariables::EQUILIBRIUM_LOADING:
    case SecondaryVariables::REACTION_DAMPING_FACTOR:
        // not handled in this method
        break;
    }

    cache.clear();
    return cache;
}


template<typename Traits>
void
LADataNoTpl<Traits>::
assembleIntegrationPoint(unsigned integration_point,
                         typename Traits::LocalMatrix const* /*localA*/,
                         typename Traits::LocalVector const* /*localRhs*/,
                         std::vector<double> const& localX,
                         const typename Traits::ShapeMatrices::ShapeType& smN,
                         const typename Traits::ShapeMatrices::DxShapeType& smDNdx,
                         const typename Traits::ShapeMatrices::JacobianType& smJ,
                         const double smDetJ,
                         const double weight)
{
    preEachAssembleIntegrationPoint(integration_point, localX, smN, smDNdx, smJ, smDetJ);

    auto const N = smDNdx.cols(); // number of integration points
    auto const D = smDNdx.rows(); // global dimension: 1, 2 or 3

    assert(N*NODAL_DOF == _Mas.cols());

    auto const laplaceCoeffMat = getLaplaceCoeffMatrix(integration_point, D);
    assert(laplaceCoeffMat.cols() == D*NODAL_DOF);
    auto const massCoeffMat    = getMassCoeffMatrix(integration_point);
    auto const advCoeffMat     = getAdvectionCoeffMatrix(integration_point);
    auto const contentCoeffMat = getContentCoeffMatrix(integration_point);


    // calculate velocity
    assert((unsigned) smDNdx.rows() == _velocity.size()
           && (unsigned) smDNdx.cols() == _velocity[0].size());

    auto const velocity = (Traits::blockDimDim(laplaceCoeffMat, 0, 0, D, D)
                           * (
                               smDNdx * Eigen::Map<const typename Traits::Vector1Comp>(localX.data(), N) // grad_p
                               / -_rho_GR
                               )
                           ).eval();
    assert(velocity.size() == D);

    for (unsigned d=0; d<D; ++d)
    {
        _velocity[d][integration_point] = velocity[d];
    }

    auto const detJ_w_N = (smDetJ * weight * smN).eval();
    auto const detJ_w_N_NT = (detJ_w_N * smN.transpose()).eval();
    assert(detJ_w_N_NT.rows() == N && detJ_w_N_NT.cols() == N);

    auto const detJ_w_N_vT_dNdx = (detJ_w_N
                                   * velocity.transpose() * smDNdx
                                   ).eval();
    assert(detJ_w_N_vT_dNdx.rows() == N && detJ_w_N_vT_dNdx.cols() == N);

    for (unsigned r=0; r<NODAL_DOF; ++r)
    {
        for (unsigned c=0; c<NODAL_DOF; ++c)
        {
            Traits::blockShpShp(_Lap_Adv_Cnt, N*r, N*c, N, N).noalias() +=
                    smDetJ * weight * smDNdx.transpose()
                    * Traits::blockDimDim(laplaceCoeffMat, D*r, D*c, D, D)
                    * smDNdx                                      // end Laplacian part
                    + detJ_w_N_NT      * contentCoeffMat(r, c)
                    + detJ_w_N_vT_dNdx * advCoeffMat(r, c);
            Traits::blockShpShp(_Mas, N*r, N*c, N, N).noalias() +=
                    detJ_w_N_NT      * massCoeffMat(r, c);
        }
    }

    auto const rhsCoeffVector = getRHSCoeffVector(integration_point);

    for (unsigned r=0; r<NODAL_DOF; ++r)
    {
        Traits::blockShp(_rhs, N*r, N).noalias() +=
                rhsCoeffVector(r) * smN * smDetJ * weight;
    }
}


template<typename Traits>
void
LADataNoTpl<Traits>::init(const unsigned num_int_pts, const unsigned dimension)
{
    _solid_density.resize(num_int_pts, _AP->_initial_solid_density);
    _solid_density_prev_ts.resize(num_int_pts, _AP->_initial_solid_density);

    _reaction_rate.resize(num_int_pts);
    _reaction_rate_prev_ts.resize(num_int_pts);

    _velocity.resize(dimension);
    for (auto& v : _velocity) v.resize(num_int_pts);

    _reaction_adaptor = std::move(TESFEMReactionAdaptor<Traits>::newInstance(*this));

    _Mas         = Traits::LocalMatrix::Zero(num_int_pts*NODAL_DOF, num_int_pts*NODAL_DOF);
    _Lap_Adv_Cnt = Traits::LocalMatrix::Zero(num_int_pts*NODAL_DOF, num_int_pts*NODAL_DOF);
    _rhs         = Traits::LocalVector::Zero(num_int_pts*NODAL_DOF);
}


template<typename Traits>
void
LADataNoTpl<Traits>::preEachAssemble()
{
    if (_AP->_iteration_in_current_timestep == 0)
    {
        if (_AP->_number_of_try_of_iteration == 0)
        {
            _solid_density_prev_ts = _solid_density;
            _reaction_rate_prev_ts = _reaction_rate;

            _reaction_adaptor->preZerothTryAssemble();
        }
        else
        {
            _solid_density = _solid_density_prev_ts;
        }
    }

    _Mas.setZero();
    _Lap_Adv_Cnt.setZero();
    _rhs.setZero();
}


template<typename Traits>
void
LADataNoTpl<Traits>
::postEachAssemble(typename Traits::LocalMatrix& localA,
                   typename Traits::LocalVector& localRhs,
                   typename Traits::LocalVector const& oldX)
{
    localA.noalias() += _Lap_Adv_Cnt + _Mas/_AP->_delta_t;
    localRhs.noalias() += _rhs
                           + _Mas * oldX/_AP->_delta_t;

    if (_AP->_output_element_matrices)
    {
        std::puts("### Element: ?");

        std::puts("---Velocity of water");
        for (auto const& vs : _velocity)
        {
            std::printf("| ");
            for (auto v : vs)
            {
                std::printf("%23.16e ", v);
            }
            std::printf("|\n");
        }

        std::printf("\nStiffness: \n");
        ogs5OutMat(localA);
        std::printf("\n");

        std::printf("\n---Mass matrix: \n");
        ogs5OutMat(_Mas);
        std::printf("\n");

        std::printf("---Laplacian + Advective + Content matrix: \n");
        ogs5OutMat(_Lap_Adv_Cnt);
        std::printf("\n");

        std::printf("---RHS: \n");
        ogs5OutVec(localRhs);
        std::printf("\n");
    }
}

} // namespace TES

} // namespace ProcessLib

#endif // PROCESS_LIB_TESFEM_DATA_IMPL_H_
