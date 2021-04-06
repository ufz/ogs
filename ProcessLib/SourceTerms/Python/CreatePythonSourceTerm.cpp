/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreatePythonSourceTerm.h"

#include <pybind11/pybind11.h>

#include "BaseLib/ConfigTree.h"
#include "MeshLib/Mesh.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/SourceTerms/SourceTerm.h"
#include "PythonSourceTerm.h"

namespace ProcessLib
{
std::unique_ptr<SourceTerm> createPythonSourceTerm(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& source_term_mesh,
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> dof_table,
    int const variable_id, int const component_id,
    unsigned const integration_order, unsigned const shapefunction_order,
    unsigned const global_dim)
{
    DBUG("Constructing PythonSourceTerm from config.");
    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__type}
    config.checkConfigParameter("type", "Python");

    auto const source_term_object =
        //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__Python__source_term_object}
        config.getConfigParameter<std::string>("source_term_object");
    //! \ogs_file_param{prj__process_variables__process_variable__source_terms__source_term__Python__flush_stdout}
    auto const flush_stdout = config.getConfigParameter("flush_stdout", false);

    // Evaluate Python code in scope of main module
    pybind11::object scope =
        pybind11::module::import("__main__").attr("__dict__");

    if (!scope.contains(source_term_object))
    {
        OGS_FATAL(
            "Function `{:s}' is not defined in the python script file, or "
            "there was no python script file specified.",
            source_term_object);
    }

    auto* source_term = scope[source_term_object.c_str()]
                            .cast<ProcessLib::SourceTerms::Python::
                                      PythonSourceTermPythonSideInterface*>();

    // In case of partitioned mesh the source_term could be empty, i.e. there is
    // no source_term condition.
#ifdef USE_PETSC
    // This can be extracted to createSourceTerm() but then the config
    // parameters are not read and will cause an error.
    // TODO (naumov): Add a function to ConfigTree for skipping the tags of the
    // subtree and move the code up in createSourceTerm().
    if (source_term_mesh.getDimension() == 0 &&
        source_term_mesh.getNumberOfNodes() == 0 &&
        source_term_mesh.getNumberOfElements() == 0)
    {
        return nullptr;
    }
#endif  // USE_PETSC

    auto const global_component_id =
        dof_table->getGlobalComponent(variable_id, component_id);
    return std::make_unique<ProcessLib::SourceTerms::Python::PythonSourceTerm>(
        std::move(dof_table),
        ProcessLib::SourceTerms::Python::PythonSourceTermData{
            source_term, global_component_id, source_term_mesh,
            source_term_mesh.getID()},
        integration_order, shapefunction_order, global_dim, flush_stdout);
}

}  // namespace ProcessLib
