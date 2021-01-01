/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "PythonSourceTerm.h"

#include <pybind11/pybind11.h>
#include <iostream>

#include "MeshLib/MeshSearch/NodeSearch.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "ProcessLib/Utils/CreateLocalAssemblers.h"
#include "PythonSourceTermLocalAssembler.h"

namespace
{
//! Optionally flushes the standard output upon creation and destruction.
//! Can be used to improve the debug output readability when printing debug
//! messages both from OGS and from Python.
class FlushStdoutGuard final
{
public:
    //! Optionally flushes C++ stdout before running Python code.
    explicit FlushStdoutGuard(bool const flush) : _flush(flush)
    {
        if (!flush)
        {
            return;
        }

        std::cout << std::flush;
    }

    //! Optionally flushes Python's stdout after running Python code.
    ~FlushStdoutGuard()
    {
        if (!_flush)
        {
            return;
        }

        using namespace pybind11::literals;
        pybind11::print("end"_a = "", "flush"_a = true);
    }

private:
    //! To flush or not to flush.
    const bool _flush;
};
}  // anonymous namespace

namespace ProcessLib
{
namespace SourceTerms
{
namespace Python
{
PythonSourceTerm::PythonSourceTerm(
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> source_term_dof_table,
    PythonSourceTermData&& source_term_data, unsigned const integration_order,
    unsigned const shapefunction_order, unsigned const global_dim,
    bool const flush_stdout)
    : SourceTerm(std::move(source_term_dof_table)),
      _source_term_data(std::move(source_term_data)),
      _flush_stdout(flush_stdout)
{
    createLocalAssemblers<PythonSourceTermLocalAssembler>(
        global_dim, _source_term_data.source_term_mesh.getElements(),
        *_source_term_dof_table, shapefunction_order, _local_assemblers,
        _source_term_data.source_term_mesh.isAxiallySymmetric(),
        integration_order, _source_term_data);
}

void PythonSourceTerm::integrate(const double t, const GlobalVector& x,
                                 GlobalVector& b, GlobalMatrix* Jac) const
{
    FlushStdoutGuard guard(_flush_stdout);

    GlobalExecutor::executeMemberOnDereferenced(
        &PythonSourceTermLocalAssemblerInterface::assemble, _local_assemblers,
        *_source_term_dof_table, t, x, b, Jac);
}

}  // namespace Python
}  // namespace SourceTerms
}  // namespace ProcessLib
