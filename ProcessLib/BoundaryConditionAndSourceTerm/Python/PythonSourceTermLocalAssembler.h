/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "ProcessLib/LocalAssemblerTraits.h"
#include "PythonSourceTerm.h"
#include "PythonSourceTermLocalAssemblerInterface.h"
#include "PythonSourceTermPythonSideInterface.h"

namespace ProcessLib
{
namespace SourceTerms
{
namespace Python
{
//! Can be thrown to indicate that a member function is not overridden in a
//! derived class (in particular, if a Python class inherits from a C++ class).
struct MethodNotOverriddenInDerivedClassException
{
};

template <typename NodalRowVectorType>
struct IntegrationPointData final
{
    IntegrationPointData(NodalRowVectorType const& N_,
                         double const& integration_weight_)
        : N(N_), integration_weight(integration_weight_)
    {
    }

    NodalRowVectorType const N;
    double const integration_weight;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
class PythonSourceTermLocalAssembler final
    : public PythonSourceTermLocalAssemblerInterface
{
public:
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;
    using NodalMatrixType = typename ShapeMatricesType::NodalMatrixType;

    explicit PythonSourceTermLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool is_axially_symmetric,
        unsigned const integration_order,
        PythonSourceTermData const& data)
        : _data(data),
          _element(e),
          _is_axially_symmetric(is_axially_symmetric),
          _integration_method(integration_order)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        auto const shape_matrices =
            NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType,
                                      GlobalDim>(
                _element, _is_axially_symmetric, _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(
                shape_matrices[ip].N,
                _integration_method.getWeightedPoint(ip).getWeight() *
                    shape_matrices[ip].integralMeasure *
                    shape_matrices[ip].detJ);
        }
    }

    void assemble(std::size_t const source_term_element_id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_source_term,
                  double const t, const GlobalVector& x, GlobalVector& b,
                  GlobalMatrix* Jac) override
    {
        using ShapeMatricesType =
            ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
        auto const fe = NumLib::createIsoparametricFiniteElement<
            ShapeFunction, ShapeMatricesType>(_element);

        unsigned const num_integration_points =
            _integration_method.getNumberOfPoints();
        auto const num_var = dof_table_source_term.getNumberOfVariables();
        auto const num_nodes = ShapeFunction::NPOINTS;
        auto const num_comp_total =
            dof_table_source_term.getNumberOfGlobalComponents();

        // gather primary variables
        typename ShapeMatricesType::template MatrixType<ShapeFunction::NPOINTS,
                                                        Eigen::Dynamic>
            primary_variables_mat(num_nodes, num_comp_total);
        for (int var = 0; var < num_var; ++var)
        {
            auto const num_comp =
                dof_table_source_term.getNumberOfVariableComponents(var);
            for (int comp = 0; comp < num_comp; ++comp)
            {
                auto const global_component =
                    dof_table_source_term.getGlobalComponent(var, comp);

                for (unsigned element_node_id = 0; element_node_id < num_nodes;
                     ++element_node_id)
                {
                    auto const boundary_node_id =
                        _element.getNode(element_node_id)->getID();
                    MeshLib::Location loc{_data.source_term_mesh_id,
                                          MeshLib::MeshItemType::Node,
                                          boundary_node_id};
                    auto const dof_idx =
                        dof_table_source_term.getGlobalIndex(loc, var, comp);
                    if (dof_idx == NumLib::MeshComponentMap::nop)
                    {
                        // TODO extend Python BC to mixed FEM ansatz functions
                        OGS_FATAL(
                            "No d.o.f. found for (node={:d}, var={:d}, "
                            "comp={:d}).  "
                            "That might be due to the use of mixed FEM ansatz "
                            "functions, which is currently not supported by "
                            "the implementation of Python BCs. That excludes, "
                            "e.g., the HM process.",
                            boundary_node_id, var, comp);
                    }
                    primary_variables_mat(element_node_id, global_component) =
                        x[dof_idx];
                }
            }
        }

        NodalRowVectorType local_rhs = Eigen::VectorXd::Zero(num_nodes);
        NodalMatrixType local_Jac =
            Eigen::MatrixXd::Zero(num_nodes, num_nodes * num_comp_total);
        std::vector<double> prim_vars_data(num_comp_total);
        auto prim_vars = MathLib::toVector(prim_vars_data);

        for (unsigned ip = 0; ip < num_integration_points; ip++)
        {
            auto const& ip_data = _ip_data[ip];
            auto const& N = ip_data.N;
            auto const& coords = fe.interpolateCoordinates(N);
            // Assumption: all primary variables have same shape functions.
            prim_vars = N * primary_variables_mat;

            auto const& w = ip_data.integration_weight;
            auto const flux_dflux =
                _data.source_term_object->getFlux(t, coords, prim_vars_data);
            auto const& flux = flux_dflux.first;
            auto const& dflux = flux_dflux.second;
            local_rhs.noalias() += N * (flux * w);

            if (static_cast<int>(dflux.size()) != num_comp_total)
            {
                // This strict check is technically mandatory only if a
                // Jacobian is assembled. However, it is done as a
                // consistency check also for cases without Jacobian
                // assembly.
                OGS_FATAL(
                    "The Python source term must return the derivative of "
                    "the flux w.r.t. each primary variable. {:d} components "
                    "expected. {:d} components returned from Python.",
                    num_comp_total, dflux.size());
            }

            if (Jac)
            {
                for (int comp = 0; comp < num_comp_total; ++comp)
                {
                    auto const top = 0;
                    auto const left = comp * num_nodes;
                    auto const width = num_nodes;
                    auto const height = num_nodes;
                    // The assignment -= takes into account the sign convention
                    // of 1st-order in time ODE systems in OpenGeoSys.
                    local_Jac.block(top, left, width, height).noalias() -=
                        ip_data.N.transpose() * (dflux[comp] * w) * N;
                }
            }
        }

        auto const& indices_specific_component =
            dof_table_source_term(source_term_element_id,
                                  _data.global_component_id)
                .rows;
        b.add(indices_specific_component, local_rhs);

        if (Jac)
        {
            // only assemble a block of the Jacobian, not the whole local matrix
            auto const indices_all_components = NumLib::getIndices(
                source_term_element_id, dof_table_source_term);
            MathLib::RowColumnIndices<GlobalIndexType> rci{
                indices_specific_component, indices_all_components};
            Jac->add(rci, local_Jac);
        }
    }

private:
    PythonSourceTermData const& _data;
    MeshLib::Element const& _element;
    bool const _is_axially_symmetric;

    IntegrationMethod const _integration_method;
    std::vector<
        IntegrationPointData<NodalRowVectorType>,
        Eigen::aligned_allocator<IntegrationPointData<NodalRowVectorType>>>
        _ip_data;
};

}  // namespace Python
}  // namespace SourceTerms
}  // namespace ProcessLib
