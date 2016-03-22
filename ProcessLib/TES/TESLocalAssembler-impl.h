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
TESLocalAssembler<ShapeFunction_,
    IntegrationMethod_,
    GlobalMatrix,
    GlobalVector,
    GlobalDim>::
TESLocalAssembler(MeshLib::Element const& e,
                  std::size_t const /*local_matrix_size*/,
                  unsigned const integration_order,
                  AssemblyParams const& asm_params)
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

    _d.setAssemblyParameters(asm_params);

    auto const n_integration_points = _shape_matrices.front().N.rows();
    _d.init(n_integration_points, GlobalDim);
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

    _d.preEachAssemble();

    for (std::size_t ip(0); ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];
        auto const& wp = integration_method.getWeightedPoint(ip);
        auto const weight = wp.getWeight();

        _d.assembleIntegrationPoint(ip, local_x,
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
getIntegrationPointValues(TESIntPtVariables const var,
                          const std::size_t /*element_id*/,
                          const GlobalVector& /*x*/,
                          std::vector<double>& cache) const
{
    return _d.getIntegrationPointValues(var, cache);
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
checkBounds(std::vector<double> const& local_x)
{
    return _d.getReactionAdaptor().checkBounds(local_x);
}

}   // namespace TES
}   // namespace ProcessLib

#endif // PROCESS_LIB_TES_FEM_IMPL_H_
