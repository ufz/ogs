/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_GROUNDWATERFLOW_FEM_H_
#define PROCESS_LIB_GROUNDWATERFLOW_FEM_H_

#include <memory>
#include <vector>

#include "Parameter.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"

#include "GroundwaterFlowAssemblyParameters.h"
#include "GroundwaterFlowProcess.h"

namespace ProcessLib
{

namespace GroundwaterFlow
{

template <typename GlobalMatrix, typename GlobalVector>
class LocalAssemblerDataInterface
{
public:
    virtual ~LocalAssemblerDataInterface() = default;

    virtual void init(MeshLib::Element const& e,
            std::size_t const local_matrix_size,
            Parameter<double, MeshLib::Element const&> const& hydraulic_conductivity,
            AssemblyParameters const& assembly_params,
            unsigned const integration_order) = 0;

    virtual void assemble(std::vector<double> const& local_x,
                          std::vector<double> const& local_x_prev_ts) = 0;

    virtual void addToGlobal(GlobalMatrix& A, GlobalVector& rhs,
            AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const&) const = 0;
};

template <typename ShapeFunction_,
         typename IntegrationMethod_,
         typename GlobalMatrix,
         typename GlobalVector,
         unsigned GlobalDim>
class LocalAssemblerData : public LocalAssemblerDataInterface<GlobalMatrix, GlobalVector>
{
public:
    using ShapeFunction = ShapeFunction_;
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;
    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    /// The hydraulic_conductivity factor is directly integrated into the local
    /// element matrix.
    void init(MeshLib::Element const& e,
              std::size_t const local_matrix_size,
              Parameter<double, MeshLib::Element const&> const&
                  hydraulic_conductivity,
              AssemblyParameters const& assembly_params,
              unsigned const integration_order) override
    {
        using FemType = NumLib::TemplateIsoparametric<
            ShapeFunction, ShapeMatricesType>;

        FemType fe(*static_cast<const typename ShapeFunction::MeshElement*>(&e));


        _integration_order = integration_order;
        IntegrationMethod_ integration_method(_integration_order);
        std::size_t const n_integration_points = integration_method.getNPoints();

        _shape_matrices.resize(n_integration_points);
        for (std::size_t ip(0); ip < n_integration_points; ip++) {
            _shape_matrices[ip].resize(ShapeFunction::DIM, ShapeFunction::NPOINTS);
            fe.computeShapeFunctions(
                    integration_method.getWeightedPoint(ip).getCoords(),
                    _shape_matrices[ip]);
        }

        _hydraulic_conductivity = [&hydraulic_conductivity, &e]()
        {
            return hydraulic_conductivity(e);
        };

        _assembly_params = &assembly_params;

        _localA.reset(new NodalMatrixType(local_matrix_size, local_matrix_size));
        _localRhs.reset(new NodalVectorType(local_matrix_size));
    }

    void assemble(std::vector<double> const& /*local_x*/,
                  std::vector<double> const& /*local_x_prev_ts*/) override
    {
        _localA->setZero();
        _localRhs->setZero();

        IntegrationMethod_ integration_method(_integration_order);
        unsigned const n_integration_points = integration_method.getNPoints();

        const double delta_t = _assembly_params->delta_t;

        for (std::size_t ip(0); ip < n_integration_points; ip++) {
            auto const& sm = _shape_matrices[ip];
            auto const& wp = integration_method.getWeightedPoint(ip);
            _localA->noalias() += sm.dNdx.transpose() *
                                  _hydraulic_conductivity() * sm.dNdx *
                                  sm.detJ * wp.getWeight();
        }
    }

    void addToGlobal(
        GlobalMatrix& A, GlobalVector& rhs,
        AssemblerLib::LocalToGlobalIndexMap::RowColumnIndices const& indices)
        const override
    {
        A.add(indices, *_localA);
        rhs.add(indices.rows, *_localRhs);
    }

private:
    std::vector<ShapeMatrices> _shape_matrices;
    std::function<double(void)> _hydraulic_conductivity;

    std::unique_ptr<NodalMatrixType> _localA;
    std::unique_ptr<NodalVectorType> _localRhs;

    unsigned _integration_order = 2;

    AssemblyParameters const* _assembly_params;
};


}   // namespace GroundwaterFlow
}   // namespace ProcessLib

#endif  // PROCESS_LIB_GROUNDWATERFLOW_FEM_H_
