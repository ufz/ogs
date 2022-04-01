/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "CollectAndInterpolateNodalDof.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NsAndWeight.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/Python/PythonBoundaryConditionPythonSideInterface.h"
#include "ProcessLib/ProcessVariable.h"

namespace ProcessLib::BoundaryConditionAndSourceTerm::Python
{

/*!
 * Common parts of the implementation of Python boundary conditions and source
 * terms.
 *
 * The implementation of Python boundary conditions and source terms are
 * almost the same, however, their interfaces differ slightly. This struct
 * contains all common parts of the local assemblers.
 */
template <typename BcOrStData, typename ShapeFunction,
          typename LowerOrderShapeFunction, typename IntegrationMethod,
          int GlobalDim>
struct BcAndStLocalAssemblerImpl
{
public:
    using Traits =
        NsAndWeightsTraits<ShapeFunction, LowerOrderShapeFunction, GlobalDim>;

private:
    using ShapeMatrixPolicy = typename Traits::ShapeMatrixPolicy;

public:
    BcAndStLocalAssemblerImpl(MeshLib::Element const& e,
                              bool is_axially_symmetric,
                              unsigned const integration_order,
                              BcOrStData const& data)
        : bc_or_st_data(data),
          element(e),
          integration_method(integration_order),
          nss_and_weights(
              computeNsAndWeights<ShapeFunction, LowerOrderShapeFunction,
                                  GlobalDim>(e, is_axially_symmetric,
                                             integration_method))
    {
    }

    //! Assembles the local rhs and (possibly) Jacobian for this BC or ST and
    //! adds them to the global rhs vector and Jacobian matrix, respectively.
    void assemble(std::size_t const boundary_element_id,
                  NumLib::LocalToGlobalIndexMap const& dof_table_boundary,
                  double const t, GlobalVector const& x, GlobalVector& b,
                  GlobalMatrix* const Jac) const
    {
        // Variables for input data
        Eigen::MatrixXd const primary_variables_mat = ProcessLib::
            BoundaryConditionAndSourceTerm::Python::collectDofsToMatrix(
                element, bc_or_st_data.bc_or_st_mesh.getID(),
                dof_table_boundary, x);

        // Variables for output data
        auto const& indices_this_component =
            dof_table_boundary(boundary_element_id,
                               bc_or_st_data.global_component_id)
                .rows;
        std::vector<GlobalIndexType> indices_all_components_for_Jac;

        auto const num_dof_this_component = indices_this_component.size();

        Eigen::VectorXd local_rhs =
            Eigen::VectorXd::Zero(num_dof_this_component);
        Eigen::MatrixXd local_Jac;

        if (Jac)
        {
            indices_all_components_for_Jac =
                NumLib::getIndices(boundary_element_id, dof_table_boundary);

            auto const num_dof_all_components =
                indices_all_components_for_Jac.size();

            local_Jac = Eigen::MatrixXd::Zero(num_dof_this_component,
                                              num_dof_all_components);
        }

        // Variables for temporary data
        auto const num_comp_total =
            dof_table_boundary.getNumberOfGlobalComponents();
        std::vector<double> prim_vars_data(num_comp_total);
        auto prim_vars = MathLib::toVector(prim_vars_data);

        // Variables for the integration loop
        unsigned const num_integration_points =
            integration_method.getNumberOfPoints();
        auto const fe = NumLib::createIsoparametricFiniteElement<
            ShapeFunction, ShapeMatrixPolicy>(element);

        for (unsigned ip = 0; ip < num_integration_points; ip++)
        {
            auto const& ns_and_weight = nss_and_weights[ip];
            auto const coords = interpolateCoords(ns_and_weight, fe);

            ProcessLib::BoundaryConditionAndSourceTerm::Python::interpolate(
                primary_variables_mat,
                bc_or_st_data.all_process_variables_for_this_process,
                nss_and_weights[ip], prim_vars);

            auto const [flag, flux, dFlux] =
                bc_or_st_data.getFlagAndFluxAndDFlux(t, coords, prim_vars_data);

            if (!flag)
            {
                // No flux value for this integration point. Skip assembly of
                // the entire element.
                return;
            }

            assembleLocalRhs(flux, ns_and_weight, local_rhs);

            if (static_cast<int>(dFlux.size()) != num_comp_total)
            {
                // This strict check is technically mandatory only if a Jacobian
                // is assembled. However, it is done as a consistency check also
                // for cases without Jacobian assembly.
                OGS_FATAL(
                    "The Python BC or ST must return the derivative of the "
                    "flux w.r.t. each component of each primary variable. "
                    "{:d} components expected. {:d} components returned from "
                    "Python."
                    "The expected number of components depends on multiple "
                    "factors. In the case of vectorial primary variables (such "
                    "as displacement), it will depend on the space dimension. "
                    "In the case of staggered coupling, the number of "
                    "components can be different for each staggered process "
                    "and different from monolithic coupling.",
                    num_comp_total, dFlux.size());
            }

            if (Jac)
            {
                assembleLocalJacobian(dFlux, ns_and_weight, local_Jac);
            }
        }

        b.add(indices_this_component, local_rhs);

        if (Jac)
        {
            MathLib::RowColumnIndices<GlobalIndexType> rci{
                indices_this_component, indices_all_components_for_Jac};
            Jac->add(rci, local_Jac);
        }
    }

private:
    //! Determines the coordinates that the point associated with the passed
    //! shape matrix by interpolating the element's node coordinates
    //!
    //! \note We use the shape function of higher order to interpolate
    //! coordinates.
    std::array<double, 3> interpolateCoords(
        typename Traits::NsAndWeight const& ns_and_weight,
        NumLib::TemplateIsoparametric<ShapeFunction, ShapeMatrixPolicy> const&
            fe) const
    {
        auto const& N_higher = ns_and_weight.NHigherOrder();
        return fe.interpolateCoordinates(N_higher);
    }

    //! Assembles the term \f$N^T \cdot \mathrm{flux} \cdot \mathrm{weight} \f$
    //! and adds it to the passed local rhs vector.
    void assembleLocalRhs(double const flux,
                          typename Traits::NsAndWeight const& ns_and_weight,
                          Eigen::VectorXd& local_rhs) const
    {
        auto const w = ns_and_weight.weight();
        auto const N = ns_and_weight.N(bc_or_st_data.shape_function_order);

        local_rhs.noalias() += N.transpose() * (flux * w);
    }

    //! Assembles the term \f$N_\mathrm{this}^T \cdot \mathrm{dFlux}_c \cdot
    //! \mathrm{weight} \cdot N_c \f$ for each component \f$c\f$ of each primary
    //! variable and adds it to the passed local Jacobian.
    //!
    //! The shape functions \f$N_\alpha\f$ of different primary variables
    //! \f$\alpha\f$ might differ, e.g., for Taylor-Hood elements.
    //!
    //! The index \em this refers to the primary variable to which this boundary
    //! condition or source term belongs. The index \f$c\f$ (called \em other in
    //! the source code in the method body) runs over all components of all
    //! primary variables of the ProcessLib::Process to which this boundary
    //! condition or source term belongs.
    void assembleLocalJacobian(
        std::vector<double> const& dFlux,
        typename Traits::NsAndWeight const& ns_and_weight,
        Eigen::MatrixXd& local_Jac) const
    {
        auto const& pv_refs =
            bc_or_st_data.all_process_variables_for_this_process;
        auto const weight = ns_and_weight.weight();

        auto const this_shp_fct_order = bc_or_st_data.shape_function_order;
        auto const N_this = ns_and_weight.N(this_shp_fct_order);

        Eigen::Index left = 0;

        for (auto pv_ref : pv_refs)
        {
            auto const& pv = pv_ref.get();
            auto const other_num_comp = pv.getNumberOfGlobalComponents();
            auto const other_shp_fct_order = pv.getShapeFunctionOrder();
            auto const N_other = ns_and_weight.N(other_shp_fct_order);
            auto const width = N_other.size();

            for (decltype(+other_num_comp) other_comp = 0;
                 other_comp < other_num_comp;
                 ++other_comp)
            {
                // The assignment -= takes into account the sign convention
                // of 1st-order in time in ODE systems in OpenGeoSys.
                local_Jac.middleCols(left, width).noalias() -=
                    N_this.transpose() * (dFlux[other_comp] * weight) * N_other;

                left += width;
            }
        }

        assert(left == local_Jac.cols());
    }

public:
    BcOrStData const& bc_or_st_data;
    MeshLib::Element const& element;

    IntegrationMethod const integration_method;
    std::vector<typename Traits::NsAndWeight> const nss_and_weights;
};

}  // namespace ProcessLib::BoundaryConditionAndSourceTerm::Python
