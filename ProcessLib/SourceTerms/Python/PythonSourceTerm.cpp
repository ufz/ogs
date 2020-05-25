/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
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
    explicit FlushStdoutGuard(bool const flush) : flush_(flush)
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
        if (!flush_)
        {
            return;
        }

        using namespace pybind11::literals;
        pybind11::print("end"_a = "", "flush"_a = true);
    }

private:
    //! To flush or not to flush.
    const bool flush_;
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
      source_term_data_(std::move(source_term_data)),
      flush_stdout_(flush_stdout)
{
    createLocalAssemblers<PythonSourceTermLocalAssembler>(
        global_dim, source_term_data_.source_term_mesh.getElements(),
        *source_term_dof_table_, shapefunction_order, local_assemblers_,
        source_term_data_.source_term_mesh.isAxiallySymmetric(),
        integration_order, source_term_data_);
}

void PythonSourceTerm::integrate(const double t, const GlobalVector& x,
                                 GlobalVector& b, GlobalMatrix* Jac) const
{
    FlushStdoutGuard guard(flush_stdout_);

    GlobalExecutor::executeMemberOnDereferenced(
        &PythonSourceTermLocalAssemblerInterface::assemble, local_assemblers_,
        *source_term_dof_table_, t, x, b, Jac);
}

}  // namespace Python
}  // namespace SourceTerms
}  // namespace ProcessLib
