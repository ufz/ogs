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
    : _shape_matrices{
          initShapeMatrices<ShapeFunction, ShapeMatricesType, IntegrationMethod_, GlobalDim>(
              e, integration_order)}
    , _d{asm_params,
         // TODO narrowing conversion
         static_cast<const unsigned>(_shape_matrices.front().N.rows()) /* number of integration points */,
         GlobalDim}
    , _local_M{ShapeFunction::NPOINTS * NODAL_DOF, ShapeFunction::NPOINTS * NODAL_DOF}
    , _local_K{ShapeFunction::NPOINTS * NODAL_DOF, ShapeFunction::NPOINTS * NODAL_DOF}
    , _local_b{ShapeFunction::NPOINTS * NODAL_DOF}
    , _integration_order{integration_order}
{
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
