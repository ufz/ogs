/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_NEUMANN_BC_ASSEMBLER_H_
#define PROCESS_LIB_NEUMANN_BC_ASSEMBLER_H_

#include <memory>
#include <vector>

#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "LocalAssemblerUtil.h"

namespace ProcessLib
{

template <typename GlobalMatrix, typename GlobalVector>
class LocalNeumannBcAsmDataInterface
{
public:
    virtual ~LocalNeumannBcAsmDataInterface() = default;

    virtual void assemble(const double t, std::vector<double>& local_b_data) = 0;
};

template <typename ShapeFunction_,
         typename IntegrationMethod_,
         typename GlobalMatrix,
         typename GlobalVector,
         unsigned GlobalDim>
class LocalNeumannBcAsmData : public LocalNeumannBcAsmDataInterface<GlobalMatrix, GlobalVector>
{
public:
    using ShapeFunction = ShapeFunction_;
    using NodalMatrixType = typename ShapeMatrixPolicyType<ShapeFunction,GlobalDim>::NodalMatrixType;
    using NodalVectorType = typename ShapeMatrixPolicyType<ShapeFunction,GlobalDim>::NodalVectorType;

    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction,GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    /// The neumann_bc_value factor is directly integrated into the local
    /// element matrix.
    LocalNeumannBcAsmData(
            MeshLib::Element const& e,
            std::size_t const local_matrix_size,
            unsigned const integration_order,
            std::function<double (MeshLib::Element const&)> const& value_lookup)
        : _local_matrix_size(local_matrix_size)
        , _integration_order(integration_order)
    {
        using FemType = NumLib::TemplateIsoparametric<
            ShapeFunction, ShapeMatricesType>;

        FemType fe(*static_cast<const typename ShapeFunction::MeshElement*>(&e));

        IntegrationMethod_ integration_method(_integration_order);
        std::size_t const n_integration_points = integration_method.getNPoints();

        _shape_matrices.reserve(n_integration_points);
        for (std::size_t ip(0); ip < n_integration_points; ip++) {
            _shape_matrices.emplace_back(ShapeFunction::DIM, GlobalDim,
                                         ShapeFunction::NPOINTS);
            fe.computeShapeFunctions(
                    integration_method.getWeightedPoint(ip).getCoords(),
                    _shape_matrices[ip]);
        }

        _neumann_bc_value = value_lookup(e);
    }

    void assemble(const double t, std::vector<double>& local_b_data) override
    {
        (void) t; // TODO time-dependent Neumann BCs

        auto local_b = setupLocalVector(local_b_data, _local_matrix_size);

        IntegrationMethod_ integration_method(_integration_order);
        std::size_t const n_integration_points = integration_method.getNPoints();

        for (std::size_t ip(0); ip < n_integration_points; ip++) {
            auto const& sm = _shape_matrices[ip];
            auto const& wp = integration_method.getWeightedPoint(ip);
            local_b.noalias() += sm.N * _neumann_bc_value
                        * sm.detJ * wp.getWeight();
        }
    }

private:
    std::vector<ShapeMatrices> _shape_matrices;
    double _neumann_bc_value;

    std::size_t const _local_matrix_size;
    unsigned const _integration_order = 2;
};

}   // namespace ProcessLib

#endif  // PROCESS_LIB_NEUMANN_BC_ASSEMBLER_H_
