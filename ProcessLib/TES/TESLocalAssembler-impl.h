/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef PROCESS_LIB_TES_FEM_IMPL_H_
#define PROCESS_LIB_TES_FEM_IMPL_H_


#include "MaterialsLib/adsorption/adsorption.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/ProcessUtil.h"

#include "TESLocalAssembler.h"
#include "TESReactionAdaptor.h"

namespace ProcessLib
{

namespace TES
{

template <typename ShapeFunction_,
          typename IntegrationMethod_,
          typename GlobalMatrix,
          typename GlobalVector,
          unsigned GlobalDim>
void
TESLocalAssembler<ShapeFunction_,
    IntegrationMethod_,
    GlobalMatrix,
    GlobalVector,
    GlobalDim>::
init(MeshLib::Element const& e,
     std::size_t const /*local_matrix_size*/,
     unsigned const integration_order, TESProcessInterface* process)
{
    _integration_order = integration_order;

    _shape_matrices =
        initShapeMatrices<ShapeFunction, ShapeMatricesType, IntegrationMethod_, GlobalDim>(
            e, integration_order);

    constexpr unsigned MAT_SIZE = ShapeFunction::NPOINTS * NODAL_DOF;

    // Resize will only do something for dynamically allocated matrices.
    _local_M.resize(MAT_SIZE, MAT_SIZE);
    _local_K.resize(MAT_SIZE, MAT_SIZE);
    _local_b.resize(MAT_SIZE);

    _data.setAssemblyParameters(process->getAssemblyParams());

    auto const n_integration_points = _shape_matrices.front().dNdx.rows();
    _data.init(n_integration_points, GlobalDim);

    _integration_point_values_cache.reset(new std::vector<double>);
}


template <typename ShapeFunction_,
          typename IntegrationMethod_,
          typename GlobalMatrix,
          typename GlobalVector,
          unsigned GlobalDim>
void
TESLocalAssembler<ShapeFunction_,
    IntegrationMethod_,
    GlobalMatrix,
    GlobalVector,
    GlobalDim>::
assemble(const double /*t*/, std::vector<double> const& local_x)
{
    _local_M.setZero();
    _local_K.setZero();
    _local_b.setZero();

    IntegrationMethod_ integration_method(_integration_order);
    unsigned const n_integration_points = integration_method.getNPoints();

    _data.preEachAssemble();

    for (std::size_t ip(0); ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];
        auto const& wp = integration_method.getWeightedPoint(ip);
        auto const weight = wp.getWeight();

        _data.assembleIntegrationPoint(ip, local_x,
                                       sm.N, sm.dNdx, sm.J, sm.detJ, weight,
                                       _local_M, _local_K, _local_b);
    }
}


template <typename ShapeFunction_,
          typename IntegrationMethod_,
          typename GlobalMatrix,
          typename GlobalVector,
          unsigned GlobalDim>
void
TESLocalAssembler<ShapeFunction_,
    IntegrationMethod_,
    GlobalMatrix,
    GlobalVector,
    GlobalDim>::
addToGlobal(
        AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const& indices,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) const
{
    M.add(indices, _local_M);
    K.add(indices, _local_K);
    b.add(indices.rows, _local_b);
}


template <typename ShapeFunction_,
          typename IntegrationMethod_,
          typename GlobalMatrix,
          typename GlobalVector,
          unsigned GlobalDim>
std::vector<double> const&
TESLocalAssembler<ShapeFunction_,
    IntegrationMethod_,
    GlobalMatrix,
    GlobalVector,
    GlobalDim>::
getIntegrationPointValues(SecondaryVariables var, NumLib::LocalNodalDOF& nodal_dof) const
{

    switch (var)
    {
    case SecondaryVariables::REACTION_RATE:
    case SecondaryVariables::SOLID_DENSITY:
    case SecondaryVariables::VELOCITY_X:
    case SecondaryVariables::VELOCITY_Y:
    case SecondaryVariables::VELOCITY_Z:
    case SecondaryVariables::LOADING:
        // These cases do not need access to nodal values
        // Thus, they can be handled inside _data
        return _data.getIntegrationPointValues(var, *_integration_point_values_cache);

    // TODO [CL] the following cases could be better provided directly using nodal values without extrapolation
    case SecondaryVariables::VAPOUR_PARTIAL_PRESSURE:
    {
        IntegrationMethod_ integration_method(_integration_order);
        auto const n_integration_points = integration_method.getNPoints();

        auto& pVs = *_integration_point_values_cache;
        pVs.clear();
        pVs.reserve(n_integration_points);

        auto const& ps = nodal_dof.getElementNodalValues(0); // TODO [CL] use constants for DOF indices
        auto const& xs = nodal_dof.getElementNodalValues(2);

        auto const& AP = _data.getAssemblyParameters();

        for (auto const& sm : _shape_matrices)
        {
            double p, xm;

            using Array = std::array<double*, 1>;
            NumLib::shapeFunctionInterpolate(ps, sm.N, Array{ &p  });
            NumLib::shapeFunctionInterpolate(xs, sm.N, Array{ &xm });

            // TODO: Dalton's law method
            auto const xn = Ads::Adsorption::get_molar_fraction(xm, AP.M_react, AP.M_inert);
            pVs.push_back(p * xn);
        }

        return pVs;
    }
    case SecondaryVariables::RELATIVE_HUMIDITY:
    {
        IntegrationMethod_ integration_method(_integration_order);
        auto const n_integration_points = integration_method.getNPoints();

        auto& rhs = *_integration_point_values_cache;
        rhs.clear();
        rhs.reserve(n_integration_points);

        auto const& nodal_vals = nodal_dof.getElementNodalValues();

        auto const& AP = _data.getAssemblyParameters();

        for (auto const& sm : _shape_matrices)
        {
            double p, T, xm;

            using Array = std::array<double*, 3>;
            NumLib::shapeFunctionInterpolate(nodal_vals, sm.N, Array{ &p, &T, &xm });

            // TODO: Dalton's law method
            auto const xn = Ads::Adsorption::get_molar_fraction(xm, AP.M_react, AP.M_inert);
            auto const pS = Ads::Adsorption::get_equilibrium_vapour_pressure(T);
            rhs.push_back(p * xn / pS);
        }

        return rhs;
    }
    case SecondaryVariables::EQUILIBRIUM_LOADING:
    {
        IntegrationMethod_ integration_method(_integration_order);
        auto const n_integration_points = integration_method.getNPoints();

        auto& Cs = *_integration_point_values_cache;
        Cs.clear();
        Cs.reserve(n_integration_points);

        auto const nodal_vals = nodal_dof.getElementNodalValues();

        auto const& AP = _data.getAssemblyParameters();

        for (auto const& sm : _shape_matrices)
        {
            double p, T, xm;

            using Array = std::array<double*, 3>;
            NumLib::shapeFunctionInterpolate(nodal_vals, sm.N, Array{ &p, &T, &xm });

            // TODO: Dalton's law method
            auto const xn = Ads::Adsorption::get_molar_fraction(xm, AP.M_react, AP.M_inert);
            auto const pV = p * xn;
            if (pV < 0.0) {
                Cs.push_back(0.0);
            } else {
                Cs.push_back(AP.react_sys->get_equilibrium_loading(pV, T, AP.M_react));
            }
        }

        return Cs;
    }
    case SecondaryVariables::REACTION_DAMPING_FACTOR:
    {
        auto& alphas = *_integration_point_values_cache;
        alphas.clear();
        alphas.resize(_shape_matrices.size(),
                      _data.getReactionAdaptor().getReactionDampingFactor());

        return alphas;
    }
    }

    _integration_point_values_cache->clear();
    return *_integration_point_values_cache;
}


template <typename ShapeFunction_,
          typename IntegrationMethod_,
          typename GlobalMatrix,
          typename GlobalVector,
          unsigned GlobalDim>
bool
TESLocalAssembler<ShapeFunction_,
    IntegrationMethod_,
    GlobalMatrix,
    GlobalVector,
    GlobalDim>::
checkBounds(std::vector<double> const& localX,
            const std::vector<double>& localX_pts)
{
    return _data.getReactionAdaptor().checkBounds(localX, localX_pts);
}

}   // namespace TES
}   // namespace ProcessLib

#endif // PROCESS_LIB_TES_FEM_IMPL_H_
