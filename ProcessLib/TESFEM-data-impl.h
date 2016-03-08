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

#include "MathLib/Nonlinear/Root1D.h"
#include "NumLib/Function/Interpolation.h"

namespace
{
const double GAS_CONST = 8.3144621;

enum class MatOutType { OGS5, PYTHON };

const MatOutType MATRIX_OUTPUT_FORMAT = MatOutType::PYTHON;

const double EQ_LOADING_FIRST_TS = -1.0;

const double RATE_CONSTANT = 6.e-3;
}


#if 0
static double fluid_specific_isobaric_heat_capacity(
		const double /*p*/, const double /*T*/, const double /*x*/)
{
	// OGS-5 heat cap model 1 value 1000
	// constant cp

	return 1000.0;
}
#endif


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
    const double p = mypow<(i>>1)>(x);
    return (i&1) ? p*p*x : p*p;
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

	/*
	const double phi_12 =   (1.0 + std::sqrt(V0_over_V1) * std::pow(1.0/M0_over_M1, 0.25))
						  * (1.0 + std::sqrt(V0_over_V1) * std::pow(1.0/M0_over_M1, 0.25))
						  / std::pow(8.0*(1.0+M0_over_M1),0.5);
						  */
	const double phi_21 = phi_12 * M0_over_M1 / V0_over_V1;

	return V0*x0 / (x0 + x1 * phi_12)
			+ V1*x1 / (x1 + x0 * phi_21);
}

static double fluid_heat_conductivity_N2(double rho, double T)
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

	double F;
	double A[5],f[9],C[4];
	double sum = 0,eta_0,c_v0,T_star,Omega = 0;
	double lamda_tr,lamda_in,lamda_r,lamda_0,lamda;

	int i;

	T_star = T * k / eps;
	rho = rho / rho_c;

	A[0] = 0.46649;
	A[1] = -0.57015;
	A[2] = 0.19164;
	A[3] = -0.03708;
	A[4] = 0.00241;

	f[0] = -0.837079888737e3;
	f[1] = 0.37914711487e2;
	f[2] = -0.601737844275;
	f[3] = 0.350418363823e1;
	f[4] = -0.874955653028e-5;
	f[5] = 0.148968607239e-7;
	f[6] = -0.256370354277e-11;
	f[7] = 0.100773735767e1;
	f[8] = 0.335340610e4;

	C[0] = 3.3373542;
	C[1] = 0.37098251;
	C[2] = 0.89913456;
	C[3] = 0.16972505;

	// dilute heat conductivity
	for (i = 0; i < 7; i++)
		sum = sum + f[i] * std::pow(T,(i - 3));
	const double temp (std::exp ((f[8] / T)) - 1);
	c_v0 = R * (sum + ((f[7] * (f[8] / T) * (f[8] / T) * (std::exp((f[8] / T)))) / (temp * temp) - 1));
	sum = 0;

	double cvint;
	cvint = c_v0 * 1000 / N_A;

	// dilute gas viscosity
	for (i = 0; i < 5; i++)
		Omega = Omega + A[i] * std::pow(log(T_star),i);
	Omega = std::exp (Omega);

	//eta in [Pa*s]
	eta_0 = 1e6 * (c1 * std::sqrt(c2 * T) / (sigma * sigma * Omega));

	F = eta_0 * k * N_A / (M * 1000);

	lamda_tr = 2.5 * (1.5 - X1);
	lamda_in = X2 * (cvint / k + X1);

	lamda_0 = F * (lamda_tr + lamda_in);
	sum = 0;
	for (i = 0; i < 4; i++)
		sum = sum + C[i] * std::pow(rho,(i + 1));

	lamda_r = sum * CCF;

	lamda = (lamda_0 + lamda_r) / 1000;   //lamda in [W/m/K]

	return lamda;
}

static double fluid_heat_conductivity_H2O(double rho, double T)
{
	double lamda,lamda_0,lamda_1,lamda_2;
	double sum1 = 0;
	double S,Q,dT;
	double a[4],b[3],B[2],d[4],C[6];
	int i;

	T = T / 647.096;
	rho = rho / 317.11;

	a[0] =  0.0102811;
	a[1] =  0.0299621;
	a[2] =  0.0156146;
	a[3] = -0.00422464;

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

	for (i = 0; i < 4; i++)
		sum1 = sum1 + a[i] * std::pow(T,i);

	lamda_0 = std::sqrt(T) * sum1;
	lamda_1 = b[0] + b[1] * rho + b[2] * std::exp(B[0] * (rho + B[1]) * (rho + B[1]));

	dT = fabs(T - 1) + C[3];
	Q = 2 + (C[4] / std::pow(dT,3. / 5.));

	if (T >= 1)
		S = 1 / dT;
	else
		S = C[5] / std::pow(dT,3. / 5.);

	lamda_2 =
	        (d[0] /
	         mypow<10>(T) + d[1]) * std::pow(rho,9. / 5.) * std::exp(C[0] * (1 - std::pow(rho,14. / 5.)))
	        + d[2]* S * std::pow(rho,Q) * std::exp((Q / (1. + Q)) * (1 - std::pow(rho,(1. + Q))))
	        + d[3] * std::exp(C[1] * std::pow(T,3. / 2.) + C[2] / mypow<5>(rho));

	lamda = (lamda_0 + lamda_1 + lamda_2); // lamda in [W/m/K]

	return lamda;
}

static double fluid_heat_conductivity(const double p, const double T, const double x)
{
	// OGS 5 fluid heat conductivity model 11

	const double M0 = ProcessLib::TES::M_N2;
	const double M1 = ProcessLib::TES::M_H2O;

	// TODO [CL] max() is redundant if the fraction is guaranteed to be between 0 and 1.
	//reactive component
	const double x0 = std::max(M0*x/(M0*x + M1*(1.0-x)), 0.); // convert mass to mole fraction
	const double k0 = fluid_heat_conductivity_H2O(M1*p/(GAS_CONST*T), T);
	//inert component
	const double x1 = 1.0 - x0;
	const double k1 = fluid_heat_conductivity_N2(M0*p/(GAS_CONST*T), T);

	const double M1_over_M2 = M1/M0; //reactive over inert
	const double V1_over_V2 = FluidViscosityH2O::get(M1*p/(GAS_CONST*T), T)
							/ FluidViscosityN2::get(M0*p/(GAS_CONST*T), T);
	const double L1_over_L2 = V1_over_V2 / M1_over_M2;

	const double phi_12 =   (1.0 + std::sqrt(L1_over_L2) * std::pow(M1_over_M2, -0.25))
						  * (1.0 + std::sqrt(V1_over_V2) * std::pow(M1_over_M2, -0.25))
						  / std::sqrt(8.0 * (1.0 + M1_over_M2));
	const double phi_21 = phi_12 * M1_over_M2 / V1_over_V2;

	return k0*x0 / (x0+x1*phi_12) + k1*x1 / (x1+x0*phi_21);
}

#if 0
static double solid_specific_isobaric_heat_capacity(const double /*rho_SR*/)
{
	// OGS 5 heat capacity model 1 (constant) value 620.573
	return 620.573;
}
#endif


namespace ProcessLib
{

namespace TES
{

template<typename Traits>
Eigen::Matrix3d
LADataNoTpl<Traits>::
getMassCoeffMatrix(const unsigned int_pt)
{
	const double dxn_dxm = _AP->_adsorption->d_molar_fraction(
							   _vapour_mass_fraction, _AP->_M_react, _AP->_M_inert);

	const double M_pp = _AP->_poro/_p * _rho_GR;
	const double M_pT = -_AP->_poro/_T *  _rho_GR;
	const double M_px = (_AP->_M_react-_AP->_M_inert) * _p
						/ (GAS_CONST * _T) * dxn_dxm * _AP->_poro
						* Trafo::dxdy(_vapour_mass_fraction);

	const double M_Tp = -_AP->_poro;
	const double M_TT = _AP->_poro * _rho_GR * _AP->_cpG // TODO: vapour heat capacity
						+ (1.0-_AP->_poro) * _solid_density[int_pt] * _AP->_cpS; // TODO: adsorbate heat capacity
	const double M_Tx = 0.0;

	const double M_xp = 0.0;
	const double M_xT = 0.0;
	const double M_xx = _AP->_poro * _rho_GR
						* Trafo::dxdy(_vapour_mass_fraction);


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

	// TODO: k_rel
	Mat L_pp = _AP->_solid_perm_tensor.block(0,0,dim,dim) * _rho_GR / eta_GR;

	Mat L_pT = Mat::Zero(dim, dim);
	Mat L_px = Mat::Zero(dim, dim);

	Mat L_Tp = Mat::Zero(dim, dim);

	// TODO: add zeolite part
	Mat L_TT = Mat::Identity(dim, dim)
					  * ( _AP->_poro * lambda_F + (1.0 - _AP->_poro) * lambda_S);

	Mat L_Tx = Mat::Zero(dim, dim);

	Mat L_xp = Mat::Zero(dim, dim);
	Mat L_xT = Mat::Zero(dim, dim);

	Mat L_xx = Mat::Identity(dim, dim)
			   * (_AP->_tortuosity * _AP->_poro * _rho_GR
				  * _AP->_diffusion_coefficient_component
				  * Trafo::dxdy(_vapour_mass_fraction)
				  );

	typename Traits::LaplaceMatrix L(dim*3, dim*3);

	L.block(    0,     0, dim, dim) = L_pp;
	L.block(    0,   dim, dim, dim) = L_pT;
	L.block(    0, 2*dim, dim, dim) = L_px;

	L.block(  dim,     0, dim, dim) = L_Tp;
	L.block(  dim,   dim, dim, dim) = L_TT;
	L.block(  dim, 2*dim, dim, dim) = L_Tx;

	L.block(2*dim,     0, dim, dim) = L_xp;
	L.block(2*dim,   dim, dim, dim) = L_xT;
	L.block(2*dim, 2*dim, dim, dim) = L_xx;

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
	const double A_xx = _rho_GR
						* Trafo::dxdy(_vapour_mass_fraction); // porosity?


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
	const double reaction_enthalpy = _AP->_adsorption->get_enthalpy(_p_V, _T, _AP->_M_react);

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
initReaction(const unsigned int_pt, const std::vector<double>& /*localX*/,
             const typename Traits::ShapeMatrices::DxShapeType& /*smDNdx*/,
             const typename Traits::ShapeMatrices::JacobianType& /*smJ*/,
             const double /*smDetJ*/)
{
    // initReaction_localDiffusionStrategy(int_pt, localX, smDNdx, smJ, smDetJ);
    // initReaction_simpleStrategy(int_pt);
    // initReaction_readjustEquilibriumLoadingStrategy(int_pt);
    initReaction_slowDownUndershootStrategy(int_pt);
}


template<typename Traits>
void LADataNoTpl<Traits>::
initReaction_localDiffusionStrategy(
        const unsigned int_pt,
        const std::vector<double> &localX,
        const typename Traits::ShapeMatrices::DxShapeType& smDNdx,
        const typename Traits::ShapeMatrices::JacobianType& smJ,
        const double smDetJ
        )
{
    if (_AP->_iteration_in_current_timestep == 0)
    {
        bool is_rate_set = false;

        // loading "at the beginning" of this timestep
        auto const loading = Ads::Adsorption::get_loading(_solid_density_prev_ts[int_pt], _AP->_rho_SR_dry);

        auto const react_rate_kinR = _AP->_adsorption->get_reaction_rate(_p_V, _T, _AP->_M_react, loading)
                                  * _AP->_rho_SR_dry;

        // calculate density change for the kinetic reaction
        auto const delta_rhoS_kin = react_rate_kinR * _AP->_delta_t * (1.0 - _AP->_poro);
        auto const delta_rhoV_kin = - delta_rhoS_kin;
        auto const rho_V = _AP->_M_react * _p_V / GAS_CONST / _T * _AP->_poro;

        bool do_equilibrium_reaction = false;
        if (delta_rhoV_kin > rho_V) // there would be more vapour released than there currently is
        {
            do_equilibrium_reaction = true;
        }
        else if (-delta_rhoV_kin > rho_V) // there would be more vapour sucked up than there currently is
        {
            // estimate how much vapour would be added to the system by diffusion

            auto const dim = smDNdx.rows();
            auto const nnodes = smDNdx.cols();

            assert(smJ.rows() == dim && smJ.cols() == dim);
            std::array<typename Traits::template Vec<0>, 3> gradients;
            for (auto& v : gradients) { v.resize(dim); }

            NumLib::shapeFunctionInterpolateGradient(localX, smDNdx, gradients);

            auto const xn = _AP->_adsorption->get_molar_fraction(
                                _vapour_mass_fraction, _AP->_M_react, _AP->_M_inert);

            auto const dxn_dxm = _AP->_adsorption->d_molar_fraction(
                                     _vapour_mass_fraction, _AP->_M_react, _AP->_M_inert);

            // TODO [CL] shouldn't nnodes be compensated for in smJDetJ?
            auto const elem_volume = smDetJ * nnodes;
            auto const elem_linear_extension = std::pow(elem_volume, 1.0/dim);
            // DBUG("elem volume: %g, ext: %g", elem_volume, elem_linear_extension);

            // diffusion current associated with p_V if temperature is assumed to be constant
            typename Traits::template Vec<0> const j_pV
                    = - _AP->_diffusion_coefficient_component *
                      ( _p * dxn_dxm * gradients[2]
                      + xn * gradients[0] );
            auto const j_pV_norm = j_pV.norm();

            auto const delta_pV_diffusion = j_pV_norm / elem_linear_extension * _AP->_delta_t;

            // DBUG("estimated delta_pV_diff: %14.7g, j_pV: %g",
            //      delta_pV_diffusion, j_pV_norm);


            auto const delta_rhoV_diffusion = delta_pV_diffusion * _AP->_M_react / GAS_CONST / _T * _AP->_poro;
            assert (delta_rhoV_diffusion >= 0.0);

            auto const rho_V_missing = - (rho_V + delta_rhoV_kin + delta_rhoV_diffusion);

            if (rho_V_missing > 0.0)
            {
                // there is not enough vapour there, even when taking diffusion into account

                // interpolate between kinetic reaction and equlibrium reaction
                auto const pV0 = _p_V;
                auto const pV = estimateAdsorptionEquilibrium(pV0, loading);

                auto const delta_pV = pV - pV0;

                // set solid density
                auto const delta_rhoS_eq  = -delta_pV * _AP->_M_react / GAS_CONST / _T * _AP->_poro;
                assert(delta_rhoS_eq > 0.0);
                auto theta = 1.0 - rho_V_missing / (delta_rhoS_kin - delta_rhoS_eq);
                theta *= 0.5 * theta;

                assert(0.0 <= theta && theta <= 1.0); // first assertion is guaranteed by above if-clause

                auto const react_rate_eqR = delta_rhoS_eq / (1.0 - _AP->_poro);
                auto const react_rate_R = theta * react_rate_kinR + (1.0 - theta) * react_rate_eqR;

                _reaction_rate[int_pt] = react_rate_R;
                _solid_density[int_pt] = _solid_density_prev_ts[int_pt] + react_rate_R * _AP->_delta_t;

                _reaction_rate_indicator[int_pt] = 100.0 * (1.0 - theta);
                _is_equilibrium_reaction[int_pt] = false;

                is_rate_set = true;

                // assert((1.0-theta)*rho_missing )
            }
        }


        if (do_equilibrium_reaction)
        {
            // in this case more water would be adsorbed than is there, even when considering diffusion
            // try equilibrium reaction

            auto const pV0 = _p_V;
            auto const pV = estimateAdsorptionEquilibrium(pV0, loading);

            auto const delta_pV = pV - pV0;
            _p += delta_pV;
            _p_V = pV;
            // set vapour mass fraction accordingly
            _vapour_mass_fraction = Ads::Adsorption::get_mass_fraction(_p_V/_p, _AP->_M_react, _AP->_M_inert);

            // set solid density
            auto const delta_rhoV = delta_pV * _AP->_M_react / GAS_CONST / _T * _AP->_poro;
            auto const delta_rhoSR = delta_rhoV / (_AP->_poro - 1.0);
            _reaction_rate[int_pt] = delta_rhoSR / _AP->_delta_t;
            _solid_density[int_pt] = _solid_density_prev_ts[int_pt] + delta_rhoSR;

            // DBUG("too much reaction: %14.7g + %14.7g + %14.7g < 0.0; dpV_diff: %14.7g", rho_V, delta_rhoV, delta_rhoV_diffusion, delta_pV_diffusion);

            _reaction_rate_indicator[int_pt] = 100.0;
            _is_equilibrium_reaction[int_pt] = false;

            is_rate_set = true;
        }
        else if (!is_rate_set) // default case, used if reaction rate from the kinetic model should be used unmodified
        {
            _reaction_rate[int_pt] = react_rate_kinR;
            _solid_density[int_pt] = _solid_density_prev_ts[int_pt] + react_rate_kinR * _AP->_delta_t;

            _reaction_rate_indicator[int_pt] = 0.0;
            _is_equilibrium_reaction[int_pt] = false;

            // TODO [CL] maybe also update p_V and p
        }
    } // first iteration

    _qR = _reaction_rate[int_pt];
}


template<typename Traits>
double LADataNoTpl<Traits>::
estimateAdsorptionEquilibrium(const double p_V0, const double C0) const
{
    auto f = [this, p_V0, C0](double pV) -> double
    {
        // pV0 := _p_V
        const double C_eq = _AP->_adsorption->get_equilibrium_loading(pV, _T, _AP->_M_react);
        return (pV - p_V0) * _AP->_M_react / GAS_CONST / _T * _AP->_poro
                + (1.0-_AP->_poro) * (C_eq - C0) * _AP->_rho_SR_dry;
    };

    // range where to search for roots of f
    const double C_eq0 = _AP->_adsorption->get_equilibrium_loading(p_V0, _T, _AP->_M_react);
    const double limit = (C_eq0 > C0)
                         ? 1e-8
                         : Ads::Adsorption::get_equilibrium_vapour_pressure(_T);

    // search for roots
    auto rf = MathLib::Nonlinear::makeRegulaFalsi<MathLib::Nonlinear::Pegasus>(f, p_V0, limit);
    rf.step(3);

    // set vapour pressure
    return rf.get_result();
}


template<typename Traits>
void LADataNoTpl<Traits>::
initReaction_simpleStrategy(const unsigned int_pt)
{
    if (_AP->_iteration_in_current_timestep == 0)
    {
        const double loading = Ads::Adsorption::get_loading(_solid_density_prev_ts[int_pt], _AP->_rho_SR_dry);

        double react_rate_R = _AP->_adsorption->get_reaction_rate(_p_V, _T, _AP->_M_react, loading)
                              * _AP->_rho_SR_dry;

        if (
            (_p_V < 100.0 || _p_V < 0.025 * Ads::Adsorption::get_equilibrium_vapour_pressure(_T))
            && react_rate_R > 0.0)
        {
            react_rate_R = 0.0;
            _reaction_rate_indicator[int_pt] = 100.0;
        }
        else
        {
            _reaction_rate_indicator[int_pt] = 0.0;
        }

        _reaction_rate[int_pt] = react_rate_R;
        _solid_density[int_pt] = _solid_density_prev_ts[int_pt] + react_rate_R * _AP->_delta_t;

        _is_equilibrium_reaction[int_pt] = false;
    }

    _qR = _reaction_rate[int_pt];
}


template<typename Traits>
void LADataNoTpl<Traits>::
initReaction_readjustEquilibriumLoadingStrategy(const unsigned int_pt)
{
    const double C_eq_curr = _AP->_adsorption->get_equilibrium_loading(_p_V, _T, _AP->_M_react);
    _equilibrium_loading[int_pt] = C_eq_curr;

    double C_eq_prev = _equilibrium_loading_prev_ts[int_pt];
    if (_AP->_iteration_in_current_timestep == 0 && C_eq_prev != EQ_LOADING_FIRST_TS) {
        // first iteration in first timestep: init equilibrium loading at previous timestep
        C_eq_prev = C_eq_curr;
        _equilibrium_loading_prev_ts[int_pt] = C_eq_prev;
    }

    const double C_0 = Ads::Adsorption::get_loading(_solid_density_prev_ts[int_pt], _AP->_rho_SR_dry);

    const double dt   = _AP->_delta_t;
    const double c_1  = C_0 - C_eq_prev + (C_eq_curr - C_eq_prev) / RATE_CONSTANT / dt;
    const double C    = c_1 * std::exp(-RATE_CONSTANT*dt) + C_eq_prev
                        + (1.0 - 1.0/RATE_CONSTANT/dt) * (C_eq_curr - C_eq_prev);
    const double Cdot = - RATE_CONSTANT * c_1 * std::exp(-RATE_CONSTANT*dt) + (C_eq_curr - C_eq_prev) / dt;

#if 0
    if (_p_V < 0.025 * Ads::Adsorption::get_equilibrium_vapour_pressure(_T))
    {
        // _solid_density[int_pt] stays constant
        _reaction_rate[int_pt] = 0.0;
    }
    else
#endif
    if ((C_eq_prev < C && C < C_eq_curr)
        || (C_eq_prev > C && C > C_eq_curr))
    {
        _solid_density[int_pt] = (1.0 + C) * _AP->_rho_SR_dry;
        _reaction_rate[int_pt] = 0.5 * (C - C_0) * _AP->_rho_SR_dry / dt
                                 + 0.5 * _reaction_rate[int_pt];
    }
    else
    {
        DBUG("C_eq_prev %14.7g, C_prev %14.7g, C_eq_curr %14.7g, C %14.7g",
             C_eq_prev, _solid_density_prev_ts[int_pt]/_AP->_rho_SR_dry - 1.0,
             C_eq_curr, C);

        _solid_density[int_pt] = (1.0 + C) * _AP->_rho_SR_dry;
        _reaction_rate[int_pt] = Cdot * _AP->_rho_SR_dry;

    }

    _qR = _reaction_rate[int_pt];
}


template<typename Traits>
void
LADataNoTpl<Traits>::
initReaction_slowDownUndershootStrategy(const unsigned int_pt)
{
    assert(_AP->_number_of_try_of_iteration < 20);

    const double loading = Ads::Adsorption::get_loading(_solid_density_prev_ts[int_pt], _AP->_rho_SR_dry);

    double react_rate_R = _AP->_adsorption->get_reaction_rate(_p_V, _T, _AP->_M_react, loading)
                          * _AP->_rho_SR_dry;

    _reaction_rate_indicator[int_pt] = 0.0;

    // set reaction rate based on current damping factor
    react_rate_R = (reaction_damping_factor > 1e-3)
                   ? reaction_damping_factor * react_rate_R
                   : 0.0;

    if (_p_V < 0.01 * Ads::Adsorption::get_equilibrium_vapour_pressure(_T)
        && react_rate_R > 0.0)
    {
        react_rate_R = 0.0;
        _reaction_rate_indicator[int_pt] = -100.0;
    }
    else if (_p_V < 100.0 || _p_V < 0.05 * Ads::Adsorption::get_equilibrium_vapour_pressure(_T))
    {
        // use equilibrium reaction for dry regime

        // in the case of zeroth try in zeroth iteration: _p_V and loading are the values
        // at the end of the previous timestep

        const double pV_eq = estimateAdsorptionEquilibrium(_p_V, loading);
        // TODO [CL]: it would be more correct to subtract pV from the previous timestep here
        const double delta_pV = pV_eq - _p_V;
        const double delta_rhoV = delta_pV * _AP->_M_react / GAS_CONST / _T * _AP->_poro;
        const double delta_rhoSR = delta_rhoV / (_AP->_poro - 1.0);
        double react_rate_R2 = delta_rhoSR/_AP->_delta_t;

        if (bounds_violation[int_pt])
        {
            react_rate_R2 *= 0.5;
        }

        // 0th try: make sure reaction is not slower than allowed by local estimation
        // nth try: make sure reaction is not made faster by local estimation
        if (
            (_AP->_number_of_try_of_iteration == 0
             && std::abs(react_rate_R2) > std::abs(react_rate_R))
            ||
            (_AP->_number_of_try_of_iteration != 0
             && std::abs(react_rate_R2) < std::abs(react_rate_R))
            )
        {
            react_rate_R = react_rate_R2;
            _reaction_rate_indicator[int_pt] = 100.0;
        }
    }

    // smooth out readjustment of reaction rate
    if (_AP->_iteration_in_current_timestep > 3)
    {
        if (_AP->_iteration_in_current_timestep <= 8)
        {
            // update reaction rate for for five iterations
            const auto N = _AP->_iteration_in_current_timestep - 3;

            // take average s.t. does not oscillate so much
            react_rate_R = 1.0 / (1.0+N) * (N*_reaction_rate[int_pt] + 1.0 * react_rate_R);
        }
        else
        {
            // afterwards no update anymore
            react_rate_R = _reaction_rate[int_pt];
        }
    }

    if (_AP->_number_of_try_of_iteration > 0)
    {
        // assert that within tries reaction does not get faster
        // (e.g. due to switch equilibrium reaction <--> kinetic reaction)

        // factor of 0.9*N: in fact, even slow down reaction over tries
        const double r = std::pow(0.9, _AP->_number_of_try_of_iteration)
                         *_reaction_rate[int_pt];
        if (std::abs(react_rate_R) > std::abs(r)) {
            react_rate_R = r;
        }
    }

    _reaction_rate[int_pt] = react_rate_R;
    _solid_density[int_pt] = _solid_density_prev_ts[int_pt] + react_rate_R * _AP->_delta_t;

    _is_equilibrium_reaction[int_pt] = false;

    _qR = _reaction_rate[int_pt];
}


template<typename Traits>
void LADataNoTpl<Traits>::
initReaction_localVapourUptakeStrategy(
        const unsigned int_pt)
{
    if (_AP->_iteration_in_current_timestep == 0)
    {
        // CAUTION: this procedure calculates the reaction rate
        // for the current timestep from the solution of the last timestep.
        // i.e. the local variables in this calculation have to be the solution of the last timestep!!!

        // first get reaction rate from the given kinetics
        // this does not consider that vapour is actually sucked up into the zeolite

        // loading "at the beginning" of this timestep
        const double loading = Ads::Adsorption::get_loading(_solid_density_prev_ts[int_pt], _AP->_rho_SR_dry);

        const double react_rate_R = _AP->_adsorption->get_reaction_rate(_p_V, _T, _AP->_M_react, loading)
                                  * _AP->_rho_SR_dry;

        // calculate density change
        const double delta_rhoS = react_rate_R * _AP->_delta_t * (1.0 - _AP->_poro);
        const double delta_rhoV   = - delta_rhoS;
        const double rho_V = _AP->_M_react * _p_V / GAS_CONST / _T * _AP->_poro;

        if (
            _p_V < 0.05 * Ads::Adsorption::get_equilibrium_vapour_pressure(_T) && (
            -delta_rhoV > rho_V   // there would be more vapour sucked up than there currently is
            || delta_rhoV > rho_V // there would be more vapour released than there currently is
        ))
        {
            // these are the corner cases where the water adsorption/desorption behaviour of the zeolite
            // controls the equilibrium to a great extent

            // in this case with the model only considering adsorption kinetics, the zeolite will adsorb (or release) more water than
            // there currently is ==> limit adsorption, use equilibrium reaction

            // function describing local equilibrium between vapour and zeolite loading
            // temperature is assumed to be constant
            auto f = [this, loading](double pV) -> double
            {
                // pV0 := _p_V
                const double C_eq = _AP->_adsorption->get_equilibrium_loading(pV, _T, _AP->_M_react);
                return (pV - _p_V) * _AP->_M_react / GAS_CONST / _T * _AP->_poro
                        + (1.0-_AP->_poro) * (C_eq - loading) * _AP->_rho_SR_dry;
            };

            // range where to search for roots of f
            const double C_eq0 = _AP->_adsorption->get_equilibrium_loading(_p_V, _T, _AP->_M_react);
            const double limit = (C_eq0 > loading)
                                 ? 1e-8
                                 : Ads::Adsorption::get_equilibrium_vapour_pressure(_T);

            // search for roots
            auto rf = MathLib::Nonlinear::makeRegulaFalsi<MathLib::Nonlinear::Pegasus>(f, _p_V, limit);
            rf.step(3);

            // set vapour pressure
            const double pV = rf.get_result();
            const double delta_pV = pV - _p_V;
            _p += delta_pV;
            _p_V = pV;
            // set vapour mass fraction accordingly
            _vapour_mass_fraction = Ads::Adsorption::get_mass_fraction(_p_V/_p, _AP->_M_react, _AP->_M_inert);

            // set solid density
            const double delta_rhoV = delta_pV * _AP->_M_react / GAS_CONST / _T * _AP->_poro;
            const double delta_rhoSR = delta_rhoV / (_AP->_poro - 1.0);
            _reaction_rate[int_pt] = delta_rhoSR / _AP->_delta_t;
            _solid_density[int_pt] = _solid_density_prev_ts[int_pt] + delta_rhoSR;

            // assert that reaction rate from this model is always slower than from the kinetic model
            // in the corner cases handled here
            assert(react_rate_R <= 0.0 || _reaction_rate[int_pt] < react_rate_R);
            assert(react_rate_R >= 0.0 || _reaction_rate[int_pt] > react_rate_R);

            DBUG("pV: %14.7g, delta pV: %14.7g, rhoSR: %14.7g, delta_rhoSR: %14.7g"
                 ", xm: %14.7g, react_rate_R: %14.7g, react_rate_GG: %14.7g"
                 ", kin rr: %14.7g",
                 _p_V, delta_pV, _solid_density[int_pt], delta_rhoSR,
                 _vapour_mass_fraction,
                 react_rate_R, _reaction_rate[int_pt],
                 _AP->_rho_SR_dry * _AP->_adsorption->get_reaction_rate(
                     _p_V, _T, _AP->_M_react,
                     Ads::Adsorption::get_loading(_solid_density_prev_ts[int_pt], _AP->_rho_SR_dry)
                     )
                 );

            _reaction_rate_indicator[int_pt] = 100.0;
            _is_equilibrium_reaction[int_pt] = true;
            _estimated_vapour_pressure[int_pt] = _p_V;
        }
        else if (_p_V < 50.0 && react_rate_R > 0.0)
        {
            _reaction_rate[int_pt] = 0.0;
            _solid_density[int_pt] = _solid_density_prev_ts[int_pt];

            _reaction_rate_indicator[int_pt] = 0.0;
            _is_equilibrium_reaction[int_pt] = false;
        }
        else
        {
            // in this case considering only the adsorption kintetics is sufficient,
            // and may be also more correct, since it really describes a kinetic.

            // TODO [CL] maybe also alter _p, _p_V, _vapour_mass_fraction

            _reaction_rate[int_pt] = react_rate_R;
            _solid_density[int_pt] = _solid_density_prev_ts[int_pt] + react_rate_R * _AP->_delta_t;

            _reaction_rate_indicator[int_pt] = 0.0;
            _is_equilibrium_reaction[int_pt] = false;
        }
    }
    else if (
             _is_equilibrium_reaction[int_pt]              // only correct the result of the equilibrium reaction
             && _AP->_iteration_in_current_timestep == 1
             )
    {
        // try to correct the reaction rate calculated in the first timestep

        const double damping = 0.5;

        // current loading in this timestep
        const double loading = Ads::Adsorption::get_loading(_solid_density[int_pt], _AP->_rho_SR_dry);

        const double react_rate_R = _AP->_adsorption->get_reaction_rate(_p_V, _T, _AP->_M_react, loading)
                                  * _AP->_rho_SR_dry;

        // calculate density change
        const double delta_rhoS = react_rate_R * _AP->_delta_t * (1.0 - _AP->_poro);
        const double delta_rhoV = - delta_rhoS;
        const double rho_V = _AP->_M_react * _p_V / GAS_CONST / _T * _AP->_poro;

        if (
            _p_V < 0.05 * Ads::Adsorption::get_equilibrium_vapour_pressure(_T) && (
            -delta_rhoV > rho_V   // there would be more vapour sucked up than there currently is
            || delta_rhoV > rho_V // there would be more vapour released than there currently is
        ))
        {
            // try equilibrium reaction again

            // function describing local equilibrium between vapour and zeolite loading
            // temperature is assumed to be constant
            auto f = [this, loading](double pV) -> double
            {
                // pV0 := _p_V
                const double C_eq = _AP->_adsorption->get_equilibrium_loading(pV, _T, _AP->_M_react);
                return (pV - _p_V) * _AP->_M_react / GAS_CONST / _T * _AP->_poro
                        + (1.0-_AP->_poro) * (C_eq - loading) * _AP->_rho_SR_dry;
            };

            // range where to search for roots of f
            const double C_eq0 = _AP->_adsorption->get_equilibrium_loading(_p_V, _T, _AP->_M_react);
            const double limit = (C_eq0 > loading)
                                 ? 1e-8
                                 : Ads::Adsorption::get_equilibrium_vapour_pressure(_T);

            // search for roots
            auto rf = MathLib::Nonlinear::makeRegulaFalsi<MathLib::Nonlinear::Pegasus>(f, _p_V, limit);
            rf.step(3);

            // set vapour pressure
            const double pV = rf.get_result();
            const double delta_pV = damping * (pV - _p_V);
            _p   += delta_pV;
            _p_V += delta_pV;
            // set vapour mass fraction accordingly
            _vapour_mass_fraction = Ads::Adsorption::get_mass_fraction(_p_V/_p, _AP->_M_react, _AP->_M_inert);

            // set solid density
            const double delta_rhoV  = delta_pV * _AP->_M_react / GAS_CONST / _T * _AP->_poro;
            const double delta_rhoSR = delta_rhoV / (_AP->_poro - 1.0);

            // add values to the reaction rate and solid density
            _reaction_rate[int_pt] += delta_rhoSR / _AP->_delta_t;
            _solid_density[int_pt] += delta_rhoSR;

            _reaction_rate_indicator[int_pt] += 50.0;
        }
        else if (_p_V < 50.0 && react_rate_R > 0.0)
        {
            // do not correct in this case
        }
        else
        {
            // in this case considering only the adsorption kintetics is sufficient,
            // and may be also more correct, since it really describes a kinetic.

            // TODO [CL] maybe also alter _p, _p_V, _vapour_mass_fraction

            _reaction_rate[int_pt] += damping * react_rate_R;
            _solid_density[int_pt] += damping * react_rate_R * _AP->_delta_t;

            _reaction_rate_indicator[int_pt] -= 50.0;
        }
    }
    else
    {
        // reaction rate does not change within a timestep
    }

    _qR = _reaction_rate[int_pt];
}


template<typename Traits>
void
LADataNoTpl<Traits>::
preEachAssembleIntegrationPoint(
        const unsigned int_pt,
        const std::vector<double> &localX,
        typename Traits::ShapeMatrices::ShapeType const& smN,
        typename Traits::ShapeMatrices::DxShapeType const& smDNdx,
        typename Traits::ShapeMatrices::JacobianType const& smJ,
        const double smDetJ)
{
#ifndef NDEBUG
    // fill local data with garbage to aid in debugging
    _p = _T = _vapour_mass_fraction = std::numeric_limits<double>::quiet_NaN();
    _p_V = _rho_GR = std::numeric_limits<double>::quiet_NaN();
    _qR = std::numeric_limits<double>::quiet_NaN();
#endif

    std::array<double*, NODAL_DOF> int_pt_val = { &_p, &_T, &_vapour_mass_fraction };

    NumLib::shapeFunctionInterpolate(localX, smN, int_pt_val);

    _vapour_mass_fraction = Trafo::x(_vapour_mass_fraction);

    // pre-compute certain properties
    _p_V = _p * Ads::Adsorption::get_molar_fraction(_vapour_mass_fraction, _AP->_M_react, _AP->_M_inert);

    initReaction(int_pt, localX, smDNdx, smJ, smDetJ);

    /*
    if (_p < 1.0) _p = 1.0;
    if (_T < 274.0) _T = 274.0;
    else if (_T > 600.0) _T = 600.0;
    if (_vapour_mass_fraction < 1e-6) _vapour_mass_fraction = 1e-6;
    else if (_vapour_mass_fraction > 1.0 - 1e-6) _vapour_mass_fraction = 1.0 - 1e-6;
    //*/

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
    case SecondaryVariables::REACTION_KINETIC_INDICATOR:
        return _reaction_rate_indicator;

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

    // assert(N*NODAL_DOF == localA.cols());
    assert(N*NODAL_DOF == _Lap->cols());

    auto const laplaceCoeffMat = getLaplaceCoeffMatrix(integration_point, D);
    assert(laplaceCoeffMat.cols() == D*NODAL_DOF);
    auto const massCoeffMat    = getMassCoeffMatrix(integration_point);
    auto const advCoeffMat     = getAdvectionCoeffMatrix(integration_point);
    auto const contentCoeffMat = getContentCoeffMatrix(integration_point);


    // calculate velocity
    assert((unsigned) smDNdx.rows() == _velocity.size() && (unsigned) smDNdx.cols() == _velocity[0].size());

    // using auto for the type went terribly wrong!
    // calculating grad_p not separately also went wrong!
    auto const grad_p = (smDNdx * Eigen::Map<const typename Traits::Vector1Comp>(localX.data(), N)).eval();
    assert(grad_p.size() == D);
    auto const velocity = (laplaceCoeffMat.block(0, 0, D, D) * grad_p
                                     / (-_rho_GR)).eval();
    assert(velocity.size() == D);

    for (unsigned d=0; d<D; ++d)
    {
        _velocity[d][integration_point] = velocity[d];
    }

    auto const detJ_w_N = (smDetJ * weight * smN).eval();
    auto const detJ_w_N_NT = (detJ_w_N * smN.transpose()).eval();
    assert(detJ_w_N_NT.rows() == N && detJ_w_N_NT.cols() == N);

    auto const vT_dNdx = (velocity.transpose() * smDNdx).eval();
    assert(vT_dNdx.cols() == N && vT_dNdx.rows() == 1);
    auto const detJ_w_N_vT_dNdx = (detJ_w_N * vT_dNdx).eval();
    assert(detJ_w_N_vT_dNdx.rows() == N && detJ_w_N_vT_dNdx.cols() == N);

    for (unsigned r=0; r<NODAL_DOF; ++r)
    {
        for (unsigned c=0; c<NODAL_DOF; ++c)
        {
            auto tmp = (smDetJ * weight * smDNdx.transpose()).eval();
            assert(tmp.cols() == D && tmp.rows() == N);
            tmp *= laplaceCoeffMat.block(D*r, D*c, D, D);
            assert(tmp.cols() == D && tmp.rows() == N);
            auto const tmp2 = (tmp * smDNdx).eval();
            assert(tmp2.cols() == N && tmp2.rows() == N);

            _Lap->block(N*r, N*c, N, N).noalias() += tmp2;
            _Mas->block(N*r, N*c, N, N).noalias() += detJ_w_N_NT      * massCoeffMat(r, c);
            _Adv->block(N*r, N*c, N, N).noalias() += detJ_w_N_vT_dNdx * advCoeffMat(r, c);
            _Cnt->block(N*r, N*c, N, N).noalias() += detJ_w_N_NT      * contentCoeffMat(r, c);
        }
    }

    auto const rhsCoeffVector = getRHSCoeffVector(integration_point);

    for (unsigned r=0; r<NODAL_DOF; ++r)
    {
        _rhs->block(N*r, 0, N, 1).noalias() +=
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

    // _velocity.resize(num_int_pts, dimension);
    _velocity.resize(dimension);
    for (auto& v : _velocity) v.resize(num_int_pts);

    _reaction_rate_indicator.resize(num_int_pts);

    _is_equilibrium_reaction.resize(num_int_pts);
    _estimated_vapour_pressure.resize(num_int_pts);

    _equilibrium_loading.resize(num_int_pts);
    _equilibrium_loading_prev_ts.resize(
                num_int_pts, EQ_LOADING_FIRST_TS); // TODO [CL] provide separate "first assembly" method

    bounds_violation.resize(num_int_pts, false);

    _Lap.reset(new typename Traits::LocalMatrix(num_int_pts*NODAL_DOF, num_int_pts*NODAL_DOF));
    _Mas.reset(new typename Traits::LocalMatrix(num_int_pts*NODAL_DOF, num_int_pts*NODAL_DOF));
    _Adv.reset(new typename Traits::LocalMatrix(num_int_pts*NODAL_DOF, num_int_pts*NODAL_DOF));
    _Cnt.reset(new typename Traits::LocalMatrix(num_int_pts*NODAL_DOF, num_int_pts*NODAL_DOF));
    _rhs.reset(new typename Traits::LocalVector(num_int_pts*NODAL_DOF));

    _Lap->setZero();
    _Mas->setZero();
    _Adv->setZero();
    _Cnt->setZero();
    _rhs->setZero();
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
            _equilibrium_loading_prev_ts = _equilibrium_loading;

            reaction_damping_factor = std::min(
                std::sqrt(reaction_damping_factor),
                10.0*reaction_damping_factor);
        }
        else
        {
            _solid_density = _solid_density_prev_ts;
        }
    }

    _Lap->setZero();
    _Mas->setZero();
    _Adv->setZero();
    _Cnt->setZero();
    _rhs->setZero();
}


template<typename Traits>
void
LADataNoTpl<Traits>
::postEachAssemble(typename Traits::LocalMatrix& localA,
                   typename Traits::LocalVector& localRhs,
                   typename Traits::LocalVector const& oldX)
{
    localA.noalias() += *_Lap + *_Mas/_AP->_delta_t + *_Adv + *_Cnt;
    localRhs.noalias() += *_rhs
                           + *_Mas * oldX/_AP->_delta_t;

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
        ogs5OutMat(*_Mas);
        std::printf("\n");

        std::printf("---Laplacian matrix: \n");
        ogs5OutMat(*_Lap);
        std::printf("\n");

        std::printf("---Advective matrix: \n");
        ogs5OutMat(*_Adv);
        std::printf("\n");

        std::printf("---Content: \n");
        ogs5OutMat(*_Cnt);
        std::printf("\n");

        std::printf("---RHS: \n");
        ogs5OutVec(localRhs);
        std::printf("\n");
    }
}

} // namespace TES

} // namespace ProcessLib

#endif // PROCESS_LIB_TESFEM_DATA_IMPL_H_
