/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#ifdef OGS_USE_PYTHON

#include <pybind11/pybind11.h>
#include "BoundaryCondition.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"

#include "GenericNaturalBoundaryConditionLocalAssembler.h"

namespace ProcessLib
{
class PythonBoundaryConditionPythonSideInterface;

//! Groups data used by essential and natural BCs, in particular by the
//! local assemblers of the latter.
struct PythonBoundaryConditionData
{
    //! Python object computing BC values.
    PythonBoundaryConditionPythonSideInterface* bc_object;

    //! DOF table of the entire domain.
    NumLib::LocalToGlobalIndexMap const& dof_table_bulk;

    //! Mesh ID of the entire domain.
    std::size_t const bulk_mesh_id;

    //! Global component ID of the (variable, component) to which this BC is
    //! applied.
    int const global_component_id;

    //! The boundary mesh, i.e., the domain of this BC.
    const MeshLib::Mesh& boundary_mesh;
};

//! A boundary condition whose values are computed by a Python script.
class PythonBoundaryCondition final : public BoundaryCondition
{
public:
    PythonBoundaryCondition(PythonBoundaryConditionData&& bc_data,
                            unsigned const integration_order,
                            unsigned const shapefunction_order,
                            unsigned const global_dim,
                            bool const flush_stdout);

    void getEssentialBCValues(
        const double t, const GlobalVector& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

    void applyNaturalBC(const double t, const GlobalVector& x, GlobalMatrix& K,
                        GlobalVector& b, GlobalMatrix* Jac) override;

private:
    //! Auxiliary data.
    PythonBoundaryConditionData _bc_data;

    //! Local dof table for the boundary mesh.
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> _dof_table_boundary;

    //! Local assemblers for all elements of the boundary mesh.
    std::vector<
        std::unique_ptr<GenericNaturalBoundaryConditionLocalAssemblerInterface>>
        _local_assemblers;

    //! Whether or not to flush standard output before and after each call to
    //! Python code. Ensures right order of output messages and therefore
    //! simplifies debugging.
    bool const _flush_stdout;
};

//! Creates a new PythonBoundaryCondition object.
std::unique_ptr<PythonBoundaryCondition> createPythonBoundaryCondition(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& boundary_mesh,
    NumLib::LocalToGlobalIndexMap const& dof_table, std::size_t bulk_mesh_id,
    int const variable_id, int const component_id,
    unsigned const integration_order, unsigned const shapefunction_order,
    unsigned const global_dim);

}  // namespace ProcessLib

#endif
