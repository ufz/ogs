/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "PythonBoundaryCondition.h"

#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/BoundaryCondition/GenericNaturalBoundaryConditionLocalAssembler.h"

#include "PythonBoundaryConditionPythonSideInterface.h"

namespace ProcessLib
{
//! Can be thrown to indicate that a member function is not overridden in a
//! derived class (in particular, if a Python class inherits from a C++ class).
struct MethodNotOverriddenInDerivedClassException
{
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class PythonBoundaryConditionLocalAssembler final
    : public GenericNaturalBoundaryConditionLocalAssembler<
          ShapeFunction, IntegrationMethod, GlobalDim>
{
    using Base = GenericNaturalBoundaryConditionLocalAssembler<
        ShapeFunction, IntegrationMethod, GlobalDim>;

public:
    PythonBoundaryConditionLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        bool is_axially_symmetric,
        unsigned const integration_order,
        PythonBoundaryConditionData const& data)
        : Base(e, is_axially_symmetric, integration_order),
          _data(data),
          _element(e)
    {
    }

    void assemble(std::size_t const boundary_element_id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, const GlobalVector& x, GlobalMatrix& /*K*/,
                  GlobalVector& b, GlobalMatrix* Jac) override
    {
        using ShapeMatricesType =
            ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
        using FemType =
            NumLib::TemplateIsoparametric<ShapeFunction, ShapeMatricesType>;
        FemType fe(*static_cast<const typename ShapeFunction::MeshElement*>(
            &_element));

        unsigned const num_integration_points =
            Base::_integration_method.getNumberOfPoints();
        auto const num_var = _data.dof_table_bulk.getNumberOfVariables();
        auto const num_nodes = _element.getNumberOfNodes();
        auto const num_comp_total =
            _data.dof_table_bulk.getNumberOfComponents();

        auto const& bulk_node_ids_map =
            *_data.boundary_mesh.getProperties()
                 .template getPropertyVector<std::size_t>("bulk_node_ids");

        // gather primary variables
        Eigen::MatrixXd primary_variables_mat(num_nodes, num_comp_total);
        for (int var = 0; var < num_var; ++var)
        {
            auto const num_comp =
                _data.dof_table_bulk.getNumberOfVariableComponents(var);
            for (int comp = 0; comp < num_comp; ++comp)
            {
                auto const global_component =
                    dof_table_boundary.getGlobalComponent(var, comp);

                for (unsigned element_node_id = 0; element_node_id < num_nodes;
                     ++element_node_id)
                {
                    auto const* const node = _element.getNode(element_node_id);
                    auto const boundary_node_id = node->getID();
                    auto const bulk_node_id =
                        bulk_node_ids_map[boundary_node_id];
                    MeshLib::Location loc{_data.bulk_mesh_id,
                                          MeshLib::MeshItemType::Node,
                                          bulk_node_id};
                    auto const dof_idx =
                        _data.dof_table_bulk.getGlobalIndex(loc, var, comp);
                    if (dof_idx == NumLib::MeshComponentMap::nop)
                    {
                        // TODO extend Python BC to mixed FEM ansatz functions
                        OGS_FATAL(
                            "No d.o.f. found for (node=%d, var=%d, comp=%d).  "
                            "That might be due to the use of mixed FEM ansatz "
                            "functions, which is currently not supported by "
                            "the implementation of Python BCs. That excludes, "
                            "e.g., the HM process.",
                            bulk_node_id, var, comp);
                    }
                    primary_variables_mat(element_node_id, global_component) =
                        x[dof_idx];
                }
            }
        }

        Eigen::VectorXd local_rhs = Eigen::VectorXd::Zero(num_nodes);
        Eigen::MatrixXd local_Jac =
            Eigen::MatrixXd::Zero(num_nodes, num_nodes * num_comp_total);
        std::vector<double> prim_vars_data(num_comp_total);
        auto prim_vars = MathLib::toVector(prim_vars_data);

        for (unsigned ip = 0; ip < num_integration_points; ip++)
        {
            auto const& sm = Base::_shape_matrices[ip];
            auto const coords = fe.interpolateCoordinates(sm.N);
            prim_vars =
                sm.N *
                primary_variables_mat;  // Assumption: all primary variables
                                        // have same shape functions.
            auto const flag_flux_dFlux =
                _data.bc_object->getFlux(t, coords, prim_vars_data);
            if (!_data.bc_object->isOverriddenNatural())
            {
                // getFlux() is not overridden in Python, so we can skip the
                // whole BC assembly (i.e., for all boundary elements).
                throw MethodNotOverriddenInDerivedClassException{};
            }

            if (!std::get<0>(flag_flux_dFlux))
            {
                // No flux value for this integration point. Skip assembly of
                // the entire element.
                return;
            }
            auto const flux = std::get<1>(flag_flux_dFlux);
            auto const& dFlux = std::get<2>(flag_flux_dFlux);

            auto const& wp = Base::_integration_method.getWeightedPoint(ip);
            auto const w = sm.detJ * wp.getWeight() * sm.integralMeasure;
            local_rhs.noalias() += sm.N * (flux * w);

            if (static_cast<int>(dFlux.size()) != num_comp_total)
            {
                // This strict check is technically mandatory only if a Jacobian
                // is assembled. However, it is done as a consistency check also
                // for cases without Jacobian assembly.
                OGS_FATAL(
                    "The Python BC must return the derivative of the flux "
                    "w.r.t. each primary variable. %d components expected. %d "
                    "components returned from Python.",
                    num_comp_total, dFlux.size());
            }

            if (Jac)
            {
                for (int comp = 0; comp < num_comp_total; ++comp)
                {
                    auto const top = 0;
                    auto const left = comp * num_nodes;
                    auto const width = num_nodes;
                    auto const height = num_nodes;
                    // The assignement -= takes into account the sign convention
                    // of 1st-order in time ODE systems in OpenGeoSys.
                    local_Jac.block(top, left, width, height).noalias() -=
                        sm.N.transpose() * (dFlux[comp] * w) * sm.N;
                }
            }
        }

        auto const& indices_specific_component =
            dof_table_boundary(boundary_element_id, _data.global_component_id)
                .rows;
        b.add(indices_specific_component, local_rhs);

        if (Jac)
        {
            // only assemble a block of the Jacobian, not the whole local matrix
            auto const indices_all_components =
                NumLib::getIndices(boundary_element_id, dof_table_boundary);
            MathLib::RowColumnIndices<GlobalIndexType> rci{
                indices_specific_component, indices_all_components};
            Jac->add(rci, local_Jac);
        }
    }

private:
    PythonBoundaryConditionData const& _data;
    MeshLib::Element const& _element;
};

}  // namespace ProcessLib
