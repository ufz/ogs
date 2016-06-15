/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef PROCESS_LIB_TES_FEM_IMPL_H_
#define PROCESS_LIB_TES_FEM_IMPL_H_

#include "MaterialsLib/Adsorption/Adsorption.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "TESLocalAssembler.h"
#include "TESReactionAdaptor.h"

namespace
{
enum class MatOutType
{
    OGS5,
    PYTHON
};

const MatOutType MATRIX_OUTPUT_FORMAT = MatOutType::PYTHON;

// TODO move to some location in the OGS core.
template <typename Mat>
void ogs5OutMat(const Mat& mat)
{
    for (unsigned r = 0; r < mat.rows(); ++r)
    {
        switch (MATRIX_OUTPUT_FORMAT)
        {
            case MatOutType::OGS5:
                if (r != 0)
                    std::printf("\n");
                std::printf("|");
                break;
            case MatOutType::PYTHON:
                if (r != 0)
                    std::printf(",\n");
                std::printf("[");
                break;
        }

        for (unsigned c = 0; c < mat.cols(); ++c)
        {
            switch (MATRIX_OUTPUT_FORMAT)
            {
                case MatOutType::OGS5:
                    std::printf(" %.16e", mat(r, c));
                    break;
                case MatOutType::PYTHON:
                    if (c != 0)
                        std::printf(",");
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

template <typename Vec>
void ogs5OutVec(const Vec& vec)
{
    for (unsigned r = 0; r < vec.size(); ++r)
    {
        switch (MATRIX_OUTPUT_FORMAT)
        {
            case MatOutType::OGS5:
                if (r != 0)
                    std::printf("\n");
                std::printf("| %.16e | ", vec[r]);
                break;
            case MatOutType::PYTHON:
                if (r != 0)
                    std::printf(",\n");
                std::printf("[ %23.16g ]", vec[r]);
                break;
        }
    }
    std::printf("\n");
}

}  // anonymous namespace

namespace ProcessLib
{
namespace TES
{
template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
TESLocalAssembler<
    ShapeFunction_, IntegrationMethod_,
    GlobalDim>::TESLocalAssembler(MeshLib::Element const& e,
                                  std::size_t const /*local_matrix_size*/,
                                  unsigned const integration_order,
                                  AssemblyParams const& asm_params)
    : _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                        IntegrationMethod_, GlobalDim>(
          e, integration_order)),
      _d(asm_params,
         // TODO narrowing conversion
         static_cast<const unsigned>(
             _shape_matrices.front()
                 .N.rows()) /* number of integration points */,
         GlobalDim),
      _local_M(ShapeFunction::NPOINTS * NODAL_DOF,
               ShapeFunction::NPOINTS * NODAL_DOF),
      _local_K(ShapeFunction::NPOINTS * NODAL_DOF,
               ShapeFunction::NPOINTS * NODAL_DOF),
      _local_b(ShapeFunction::NPOINTS * NODAL_DOF),
      _integration_order(integration_order)
{
}

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
void TESLocalAssembler<ShapeFunction_, IntegrationMethod_,
                       GlobalDim>::assemble(const double /*t*/,
                                            std::vector<double> const& local_x)
{
    _local_M.setZero();
    _local_K.setZero();
    _local_b.setZero();

    IntegrationMethod_ integration_method(_integration_order);
    unsigned const n_integration_points = integration_method.getNumberOfPoints();

    _d.preEachAssemble();

    for (std::size_t ip(0); ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];
        auto const& wp = integration_method.getWeightedPoint(ip);
        auto const weight = wp.getWeight();

        _d.assembleIntegrationPoint(ip, local_x, sm.N, sm.dNdx, sm.J, sm.detJ,
                                    weight, _local_M, _local_K, _local_b);
    }

    if (_d.getAssemblyParameters().output_element_matrices)
    {
        std::puts("### Element: ?");

        std::puts("---Velocity of water");
        for (auto const& vs : _d.getData().velocity)
        {
            std::printf("| ");
            for (auto v : vs)
            {
                std::printf("%23.16e ", v);
            }
            std::printf("|\n");
        }

        std::printf("\n---Mass matrix: \n");
        ogs5OutMat(_local_M);
        std::printf("\n");

        std::printf("---Laplacian + Advective + Content matrix: \n");
        ogs5OutMat(_local_K);
        std::printf("\n");

        std::printf("---RHS: \n");
        ogs5OutVec(_local_b);
        std::printf("\n");
    }
}

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
void TESLocalAssembler<ShapeFunction_, IntegrationMethod_, GlobalDim>::
    addToGlobal(
        NumLib::LocalToGlobalIndexMap::RowColumnIndices const& indices,
        GlobalMatrix& M, GlobalMatrix& K, GlobalVector& b) const
{
    M.add(indices, _local_M);
    K.add(indices, _local_K);
    b.add(indices.rows, _local_b);
}

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
std::vector<double> const& TESLocalAssembler<
    ShapeFunction_, IntegrationMethod_,
    GlobalDim>::getIntegrationPointValues(TESIntPtVariables const var,
                                          std::vector<double>& cache) const
{
    return _d.getIntegrationPointValues(var, cache);
}

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
bool TESLocalAssembler<
    ShapeFunction_, IntegrationMethod_,
    GlobalDim>::checkBounds(std::vector<double> const& local_x,
                            std::vector<double> const& local_x_prev_ts)
{
    return _d.getReactionAdaptor().checkBounds(local_x, local_x_prev_ts);
}

}  // namespace TES
}  // namespace ProcessLib

#endif  // PROCESS_LIB_TES_FEM_IMPL_H_
