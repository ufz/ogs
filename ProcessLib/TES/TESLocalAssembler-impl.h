/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "MaterialLib/Adsorption/Adsorption.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
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
TESLocalAssembler<ShapeFunction_, IntegrationMethod_, GlobalDim>::
    TESLocalAssembler(MeshLib::Element const& e,
                      std::size_t const /*local_matrix_size*/,
                      bool is_axially_symmetric,
                      unsigned const integration_order,
                      AssemblyParams const& asm_params)
    : _integration_method(integration_order),
      _shape_matrices(initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                        IntegrationMethod_, GlobalDim>(
          e, is_axially_symmetric, _integration_method)),
      _d(asm_params, _integration_method.getNumberOfPoints(), GlobalDim)
{
}

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
void TESLocalAssembler<ShapeFunction_, IntegrationMethod_, GlobalDim>::assemble(
    double const /*t*/, std::vector<double> const& local_x,
    std::vector<double>& local_M_data, std::vector<double>& local_K_data,
    std::vector<double>& local_b_data)
{
    auto const local_matrix_size = local_x.size();
    // This assertion is valid only if all nodal d.o.f. use the same shape matrices.
    assert(local_matrix_size == ShapeFunction::NPOINTS * NODAL_DOF);

    auto local_M = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<NodalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
    auto local_b = MathLib::createZeroedVector<NodalVectorType>(local_b_data,
                                                            local_matrix_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    _d.preEachAssemble();

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = _shape_matrices[ip];
        auto const& wp = _integration_method.getWeightedPoint(ip);
        auto const weight = wp.getWeight();

        _d.assembleIntegrationPoint(ip, local_x, sm, weight, local_M, local_K,
                                    local_b);
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
        ogs5OutMat(local_M);
        std::printf("\n");

        std::printf("---Laplacian + Advective + Content matrix: \n");
        ogs5OutMat(local_K);
        std::printf("\n");

        std::printf("---RHS: \n");
        ogs5OutVec(local_b);
        std::printf("\n");
    }
}

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
std::vector<double> const& TESLocalAssembler<
    ShapeFunction_, IntegrationMethod_,
    GlobalDim>::getIntPtSolidDensity(std::vector<double>& /*cache*/) const
{
    return _d.getData().solid_density;
}

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
std::vector<double> const& TESLocalAssembler<
    ShapeFunction_, IntegrationMethod_,
    GlobalDim>::getIntPtLoading(std::vector<double>& cache) const
{
    auto const rho_SR = _d.getData().solid_density;
    auto const rho_SR_dry = _d.getAssemblyParameters().rho_SR_dry;

    cache.clear();
    cache.reserve(rho_SR.size());

    for (auto const rho : rho_SR) {
        cache.push_back(rho/rho_SR_dry - 1.0);
    }

    return cache;
}

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
std::vector<double> const&
TESLocalAssembler<ShapeFunction_, IntegrationMethod_, GlobalDim>::
    getIntPtReactionDampingFactor(std::vector<double>& cache) const
{
    auto const fac = _d.getData().reaction_adaptor->getReactionDampingFactor();
    auto const num_integration_points = _d.getData().solid_density.size();

    cache.clear();
    cache.resize(num_integration_points, fac);

    return cache;
}

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
std::vector<double> const& TESLocalAssembler<
    ShapeFunction_, IntegrationMethod_,
    GlobalDim>::getIntPtReactionRate(std::vector<double>& /*cache*/) const
{
    return _d.getData().reaction_rate;
}

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
std::vector<double> const& TESLocalAssembler<
    ShapeFunction_, IntegrationMethod_,
    GlobalDim>::getIntPtDarcyVelocityX(std::vector<double>& /*cache*/) const
{
    return _d.getData().velocity[0];
}

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
std::vector<double> const& TESLocalAssembler<
    ShapeFunction_, IntegrationMethod_,
    GlobalDim>::getIntPtDarcyVelocityY(std::vector<double>& /*cache*/) const
{
    assert(_d.getData().velocity.size() > 1);
    return _d.getData().velocity[1];
}

template <typename ShapeFunction_, typename IntegrationMethod_,
          unsigned GlobalDim>
std::vector<double> const& TESLocalAssembler<
    ShapeFunction_, IntegrationMethod_,
    GlobalDim>::getIntPtDarcyVelocityZ(std::vector<double>& /*cache*/) const
{
    assert(_d.getData().velocity.size() > 2);
    return _d.getData().velocity[2];
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
