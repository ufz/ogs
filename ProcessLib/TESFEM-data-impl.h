/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
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

#include <logog/include/logog.hpp>

#include "NumLib/Function/Interpolation.h"

#include "TESFEM-data-fwd.h"
#include "TESFEMReactionAdaptor.h"
#include "TESOGS5MaterialModels.h"

namespace
{
enum class MatOutType { OGS5, PYTHON };

const MatOutType MATRIX_OUTPUT_FORMAT = MatOutType::PYTHON;
}

namespace ProcessLib
{

namespace TES
{

template<typename Traits>
Eigen::Matrix3d
TESLocalAssemblerInner<Traits>::
getMassCoeffMatrix(const unsigned int_pt)
{
	// TODO: Dalton's law property
	const double dxn_dxm = Ads::Adsorption::d_molar_fraction(
							   _d.vapour_mass_fraction, _d.ap->M_react, _d.ap->M_inert);

	const double M_pp = _d.ap->poro/_d.p * _d.rho_GR;
	const double M_pT = -_d.ap->poro/_d.T *  _d.rho_GR;
	const double M_px = (_d.ap->M_react-_d.ap->M_inert) * _d.p
						/ (GAS_CONST * _d.T) * dxn_dxm * _d.ap->poro;

	const double M_Tp = -_d.ap->poro;
	const double M_TT =
			_d.ap->poro * _d.rho_GR * _d.ap->cpG // TODO: vapour heat capacity
			+ (1.0-_d.ap->poro) * _d.solid_density[int_pt] * _d.ap->cpS; // TODO: adsorbate heat capacity
	const double M_Tx = 0.0;

	const double M_xp = 0.0;
	const double M_xT = 0.0;
	const double M_xx = _d.ap->poro * _d.rho_GR;


	Eigen::Matrix3d M;
	M << M_pp, M_pT, M_px,
		 M_Tp, M_TT, M_Tx,
		 M_xp, M_xT, M_xx;

	return M;
}


template<typename Traits>
typename Traits::LaplaceMatrix
TESLocalAssemblerInner<Traits>::
getLaplaceCoeffMatrix(const unsigned /*int_pt*/, const unsigned dim)
{
	const double eta_GR = fluid_viscosity(_d.p, _d.T, _d.vapour_mass_fraction);

	const double lambda_F = fluid_heat_conductivity(_d.p, _d.T, _d.vapour_mass_fraction);
	const double lambda_S = _d.ap->solid_heat_cond;

	using Mat = typename Traits::MatrixDimDim;

	typename Traits::LaplaceMatrix L
			= Traits::LaplaceMatrix::Zero(dim*NODAL_DOF, dim*NODAL_DOF);

	// TODO: k_rel
	// L_pp
	Traits::blockDimDim(L,     0,     0, dim, dim)
			= Traits::blockDimDim(_d.ap->solid_perm_tensor, 0,0,dim,dim) * _d.rho_GR / eta_GR;

	// TODO: add zeolite part
	// L_TT
	Traits::blockDimDim(L,   dim,   dim, dim, dim)
			= Mat::Identity(dim, dim)
			  * ( _d.ap->poro * lambda_F + (1.0 - _d.ap->poro) * lambda_S);

	// L_xx
	Traits::blockDimDim(L, 2*dim, 2*dim, dim, dim)
			= Mat::Identity(dim, dim)
			  * (_d.ap->tortuosity * _d.ap->poro * _d.rho_GR
				 * _d.ap->diffusion_coefficient_component
				 );

	return L;
}


template<typename Traits>
Eigen::Matrix3d
TESLocalAssemblerInner<Traits>::
getAdvectionCoeffMatrix(const unsigned /*int_pt*/)
{
	const double A_pp = 0.0;
	const double A_pT = 0.0;

	const double A_px = 0.0;

	const double A_Tp = 0.0;

	const double A_TT = _d.rho_GR * _d.ap->cpG; // porosity?
	const double A_Tx = 0.0;

	const double A_xp = 0.0;
	const double A_xT = 0.0;
	const double A_xx = _d.rho_GR; // porosity?


	Eigen::Matrix3d A;
	A << A_pp, A_pT, A_px,
		 A_Tp, A_TT, A_Tx,
		 A_xp, A_xT, A_xx;

	return A;
}


template<typename Traits>
Eigen::Matrix3d
TESLocalAssemblerInner<Traits>::
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
	const double C_xx = (_d.ap->poro - 1.0) * _d.qR;

	Eigen::Matrix3d C;
	C << C_pp, C_pT, C_px,
		 C_Tp, C_TT, C_Tx,
		 C_xp, C_xT, C_xx;

	return C;
}


template<typename Traits>
Eigen::Vector3d
TESLocalAssemblerInner<Traits>::
getRHSCoeffVector(const unsigned int_pt)
{
	const double reaction_enthalpy = _d.ap->react_sys->get_enthalpy(_d.p_V, _d.T, _d.ap->M_react);

	const double rhs_p = (_d.ap->poro - 1.0) * _d.qR; // TODO [CL] body force term

	const double rhs_T = _d.rho_GR * _d.ap->poro * _d.ap->fluid_specific_heat_source
						 + (1.0 - _d.ap->poro) * _d.qR * reaction_enthalpy
						 + _d.solid_density[int_pt] * (1.0 - _d.ap->poro) * _d.ap->solid_specific_heat_source;
						 // TODO [CL] momentum production term

	const double rhs_x = (_d.ap->poro - 1.0) * _d.qR; // TODO [CL] what if x < 0.0


	Eigen::Vector3d rhs;
	rhs << rhs_p,
		 rhs_T,
		 rhs_x;

	return rhs;
}


template<typename Traits>
void
TESLocalAssemblerInner<Traits>::
initReaction(const unsigned int_pt)
{
    auto const& rate = _d.reaction_adaptor->initReaction(int_pt);
    _d.qR = rate.reaction_rate;
    _d.reaction_rate[int_pt] = rate.reaction_rate;
    _d.solid_density[int_pt] = rate.solid_density;
}


template<typename Traits>
void
TESLocalAssemblerInner<Traits>::
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
    _d.p = _d.T   = _d.vapour_mass_fraction
       = _d.p_V = _d.rho_GR
       = _d.qR
       = std::numeric_limits<double>::quiet_NaN();
#endif

    std::array<double*, NODAL_DOF> int_pt_val = { &_d.p, &_d.T, &_d.vapour_mass_fraction };

    NumLib::shapeFunctionInterpolate(localX, smN, int_pt_val);

    // pre-compute certain properties
    _d.p_V = _d.p * Ads::Adsorption::get_molar_fraction(_d.vapour_mass_fraction, _d.ap->M_react, _d.ap->M_inert);

    initReaction(int_pt);

    assert(_d.p > 0.0);
    assert(_d.T > 0.0);
    assert(0.0 <= _d.vapour_mass_fraction && _d.vapour_mass_fraction <= 1.0);

    _d.rho_GR = fluid_density(_d.p, _d.T, _d.vapour_mass_fraction);
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
TESLocalAssemblerInner<Traits>::
getIntegrationPointValues(SecondaryVariables var, std::vector<double>& cache) const
{
    switch (var)
    {
    case SecondaryVariables::REACTION_RATE:
        return _d.reaction_rate;
    case SecondaryVariables::SOLID_DENSITY:
        return _d.solid_density;
    case SecondaryVariables::VELOCITY_X:
        return _d.velocity[0];
    case SecondaryVariables::VELOCITY_Y:
        assert(_d.velocity.size() >= 2);
        return _d.velocity[1];
    case SecondaryVariables::VELOCITY_Z:
        assert(_d.velocity.size() >= 3);
        return _d.velocity[2];

    case SecondaryVariables::LOADING:
    {
        auto& Cs = cache;
        Cs.clear();
        Cs.reserve(_d.solid_density.size());

        for (auto rho_SR : _d.solid_density) {
            Cs.push_back(rho_SR / _d.ap->rho_SR_dry - 1.0);
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
TESLocalAssemblerInner<Traits>::
assembleIntegrationPoint(unsigned integration_point,
                         std::vector<double> const& localX,
                         const typename Traits::ShapeMatrices::ShapeType& smN,
                         const typename Traits::ShapeMatrices::DxShapeType& smDNdx,
                         const typename Traits::ShapeMatrices::JacobianType& smJ,
                         const double smDetJ,
                         const double weight,
                         typename Traits::LocalMatrix& local_M,
                         typename Traits::LocalMatrix& local_K,
                         typename Traits::LocalVector& local_b)
{
    preEachAssembleIntegrationPoint(integration_point, localX, smN, smDNdx, smJ, smDetJ);

    auto const N = smDNdx.cols(); // number of integration points
    auto const D = smDNdx.rows(); // global dimension: 1, 2 or 3

    assert(N*NODAL_DOF == local_M.cols());

    auto const laplaceCoeffMat = getLaplaceCoeffMatrix(integration_point, D);
    assert(laplaceCoeffMat.cols() == D*NODAL_DOF);
    auto const massCoeffMat    = getMassCoeffMatrix(integration_point);
    auto const advCoeffMat     = getAdvectionCoeffMatrix(integration_point);
    auto const contentCoeffMat = getContentCoeffMatrix(integration_point);


    // calculate velocity
    assert((unsigned) smDNdx.rows() == _d.velocity.size()
           && (unsigned) smDNdx.cols() == _d.velocity[0].size());

    auto const velocity = (Traits::blockDimDim(laplaceCoeffMat, 0, 0, D, D)
                           * (
                               smDNdx * Eigen::Map<const typename Traits::Vector1Comp>(localX.data(), N) // grad_p
                               / -_d.rho_GR
                               )
                           ).eval();
    assert(velocity.size() == D);

    for (unsigned d=0; d<D; ++d)
    {
        _d.velocity[d][integration_point] = velocity[d];
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
            Traits::blockShpShp(local_K, N*r, N*c, N, N).noalias() +=
                    smDetJ * weight * smDNdx.transpose()
                    * Traits::blockDimDim(laplaceCoeffMat, D*r, D*c, D, D)
                    * smDNdx                                      // end Laplacian part
                    + detJ_w_N_NT      * contentCoeffMat(r, c)
                    + detJ_w_N_vT_dNdx * advCoeffMat(r, c);
            Traits::blockShpShp(local_M, N*r, N*c, N, N).noalias() +=
                    detJ_w_N_NT      * massCoeffMat(r, c);
        }
    }

    auto const rhsCoeffVector = getRHSCoeffVector(integration_point);

    for (unsigned r=0; r<NODAL_DOF; ++r)
    {
        Traits::blockShp(local_b, N*r, N).noalias() +=
                rhsCoeffVector(r) * smN * smDetJ * weight;
    }
}


template<typename Traits>
void
TESLocalAssemblerInner<Traits>::
init(const unsigned num_int_pts, const unsigned dimension)
{
    _d.solid_density.resize(num_int_pts, _d.ap->initial_solid_density);
    _d.solid_density_prev_ts.resize(num_int_pts, _d.ap->initial_solid_density);

    _d.reaction_rate.resize(num_int_pts);
    _d.reaction_rate_prev_ts.resize(num_int_pts);

    _d.velocity.resize(dimension);
    for (auto& v : _d.velocity) v.resize(num_int_pts);

    _d.reaction_adaptor = std::move(TESFEMReactionAdaptor::newInstance(_d));
}


template<typename Traits>
void
TESLocalAssemblerInner<Traits>::preEachAssemble()
{
    if (_d.ap->iteration_in_current_timestep == 0)
    {
        if (_d.ap->number_of_try_of_iteration == 0)
        {
            _d.solid_density_prev_ts = _d.solid_density;
            _d.reaction_rate_prev_ts = _d.reaction_rate;

            _d.reaction_adaptor->preZerothTryAssemble();
        }
        else
        {
            _d.solid_density = _d.solid_density_prev_ts;
        }
    }
}


#if 0
template<typename Traits>
void
TESLocalAssemblerInner<Traits>
::postEachAssemble()
{
    if (_d.ap->_output_element_matrices)
    {
        std::puts("### Element: ?");

        std::puts("---Velocity of water");
        for (auto const& vs : _d.velocity)
        {
            std::printf("| ");
            for (auto v : vs)
            {
                std::printf("%23.16e ", v);
            }
            std::printf("|\n");
        }

        // TODO meaning has changed. Not the same as in OGS5 anymore!
        std::printf("\nStiffness: \n");
        ogs5OutMat(local_K);
        std::printf("\n");

        std::printf("\n---Mass matrix: \n");
        ogs5OutMat(_d.Mas);
        std::printf("\n");

        std::printf("---Laplacian + Advective + Content matrix: \n");
        ogs5OutMat(_d.Lap_Adv_Cnt);
        std::printf("\n");

        std::printf("---RHS: \n");
        ogs5OutVec(local_b);
        std::printf("\n");
    }
}
#endif

} // namespace TES

} // namespace ProcessLib

#endif // PROCESS_LIB_TESFEM_DATA_IMPL_H_
