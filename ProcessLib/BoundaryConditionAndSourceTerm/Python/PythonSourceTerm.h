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

#include "ProcessLib/BoundaryConditionAndSourceTerm/SourceTerm.h"
#include "PythonSourceTermLocalAssemblerInterface.h"
#include "PythonSourceTermPythonSideInterface.h"

namespace ProcessLib
{
class LocalToGlobalIndexMap;
}

namespace ProcessLib
{
namespace SourceTerms
{
namespace Python
{
//! Groups data used by source terms, in particular by the local assemblers.
struct PythonSourceTermData final
{
    //! Python object computing source term values.
    PythonSourceTermPythonSideInterface* source_term_object;

    //! Global component ID of the (variable, component) to which this source
    //! term is applied.
    int const global_component_id;

    //! The source term mesh, i.e., the (sub-) domain of this source term.
    const MeshLib::Mesh& source_term_mesh;

    //! Mesh ID of the entire domain.
    std::size_t const source_term_mesh_id;
};

//! A source term whose values are computed by a Python script.
class PythonSourceTerm final : public ProcessLib::SourceTerm
{
public:
    explicit PythonSourceTerm(
        std::unique_ptr<NumLib::LocalToGlobalIndexMap> source_term_dof_table,
        PythonSourceTermData&& source_term_data,
        unsigned const integration_order, unsigned const shapefunction_order,
        unsigned const global_dim, bool const flush_stdout);

    void integrate(const double t, GlobalVector const& x, GlobalVector& b,
                   GlobalMatrix* jac) const override;

private:
    //! Auxiliary data.
    PythonSourceTermData _source_term_data;

    //! Local assemblers for all elements of the source term mesh.
    std::vector<std::unique_ptr<PythonSourceTermLocalAssemblerInterface>>
        _local_assemblers;

    //! Whether or not to flush standard output before and after each call to
    //! Python code. Ensures right order of output messages and therefore
    //! simplifies debugging.
    bool const _flush_stdout;
};

}  // namespace Python
}  // namespace SourceTerms
}  // namespace ProcessLib
