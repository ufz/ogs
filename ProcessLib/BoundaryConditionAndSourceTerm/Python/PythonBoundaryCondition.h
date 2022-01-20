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

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "ProcessLib/BoundaryConditionAndSourceTerm/BoundaryCondition.h"
#include "PythonBoundaryConditionLocalAssemblerInterface.h"
#include "PythonBoundaryConditionPythonSideInterface.h"
#include "Utils/BcOrStData.h"

namespace ProcessLib
{
class ProcessVariable;

//! Can be thrown to indicate that a member function is not overridden in a
//! derived class (in particular, if a Python class inherits from a C++ class).
struct MethodNotOverriddenInDerivedClassException
{
};

struct PythonBcData final
    : ProcessLib::BoundaryConditionAndSourceTerm::Python::BcOrStData<
          PythonBoundaryConditionPythonSideInterface>
{
    //! The interfaces of Python BCs and STs differ slightly.
    //!
    //! This method provides an interface for accessing Python BC and
    //! ST objects in a uniform way.
    ProcessLib::BoundaryConditionAndSourceTerm::Python::FlagAndFluxAndDFlux
    getFlagAndFluxAndDFlux(double const t, std::array<double, 3> const coords,
                           std::vector<double> const& prim_vars_data) const
    {
        auto [flag, flux, dFlux] =
            bc_or_st_object->getFlux(t, coords, prim_vars_data);

        if (!bc_or_st_object->isOverriddenNatural())
        {
            // getFlux() is not overridden in Python, so we can skip the
            // whole BC assembly (i.e., for all boundary elements).
            throw MethodNotOverriddenInDerivedClassException{};
        }

        return {flag, flux, std::move(dFlux)};
    }
};

//! A boundary condition whose values are computed by a Python script.
class PythonBoundaryCondition final : public BoundaryCondition
{
public:
    PythonBoundaryCondition(
        PythonBcData&& bc_data,
        unsigned const integration_order,
        bool const flush_stdout,
        unsigned const bulk_mesh_dimension,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk);

    void getEssentialBCValues(
        const double t, const GlobalVector& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

    void applyNaturalBC(const double t, std::vector<GlobalVector*> const& x,
                        int const process_id, GlobalMatrix& K, GlobalVector& b,
                        GlobalMatrix* Jac) override;

private:
    //! Collects primary variables at the passed node from the passed
    //! GlobalVector to \c primary_variables.
    //!
    //! Primary variables at higher order nodes are interpolated from base nodes
    //! if necessary, e.g., for Taylor-Hood elements.
    //!
    //! \post \c primary_variables will contain the value of each component of
    //! each primary variable. Their order is determined by the d.o.f. table.
    void collectPrimaryVariables(std::vector<double>& primary_variables,
                                 MeshLib::Node const& boundary_node,
                                 GlobalVector const& x) const;

    //! Get the d.o.f. index at the given \c boundary_node_id for this BC's
    //! variable and component.
    GlobalIndexType getDofIdx(std::size_t const boundary_node_id) const;

    //! Get the d.o.f. index at the given \c boundary_node_id for the given
    //! variable and component.
    GlobalIndexType getDofIdx(std::size_t const boundary_node_id, int const var,
                              int const comp) const;

    //! Interpolates the given component of the given variable to the given \c
    //! boundary_node.
    double interpolateToHigherOrderNode(
        GlobalVector const& x, int const var, int const comp,
        MeshLib::Node const& boundary_node) const;

    //! Auxiliary data used by the local assemblers.
    PythonBcData _bc_data;

    //! Local dof table for the boundary mesh.
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> _dof_table_boundary;

    //! Local assemblers for all elements of the boundary mesh.
    std::vector<std::unique_ptr<PythonBoundaryConditionLocalAssemblerInterface>>
        _local_assemblers;

    //! Whether or not to flush standard output before and after each call to
    //! Python code. Ensures right order of output messages and therefore
    //! simplifies debugging.
    bool const _flush_stdout;
};

//! Creates a new PythonBoundaryCondition object.
std::unique_ptr<PythonBoundaryCondition> createPythonBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& boundary_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
    MeshLib::Mesh const& bulk_mesh, int const variable_id,
    int const component_id, unsigned const integration_order,
    unsigned const shapefunction_order,
    std::vector<std::reference_wrapper<ProcessVariable>> const&
        all_process_variables_for_this_process);

}  // namespace ProcessLib
