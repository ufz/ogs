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
#include "Utils/BcOrStData.h"

namespace ProcessLib
{
class ProcessVariable;

namespace SourceTerms::Python
{
struct PythonStData final
    : ProcessLib::BoundaryConditionAndSourceTerm::Python::BcOrStData<
          PythonSourceTermPythonSideInterface>
{
    //! The interfaces of Python BCs and STs differ slightly.
    //!
    //! This method provides an interface for accessing Python BC and
    //! ST objects in a uniform way.
    ProcessLib::BoundaryConditionAndSourceTerm::Python::FlagAndFluxAndDFlux
    getFlagAndFluxAndDFlux(double const t, std::array<double, 3> const coords,
                           std::vector<double> const& prim_vars_data) const
    {
        auto [flux, dFlux] =
            bc_or_st_object->getFlux(t, coords, prim_vars_data);

        return {true, flux, std::move(dFlux)};
    }
};

//! A source term whose values are computed by a Python script.
class PythonSourceTerm final : public ProcessLib::SourceTerm
{
public:
    explicit PythonSourceTerm(
        std::unique_ptr<NumLib::LocalToGlobalIndexMap> source_term_dof_table,
        PythonStData&& source_term_data, unsigned const integration_order,
        unsigned const global_dim, bool const flush_stdout);

    void integrate(const double t, GlobalVector const& x, GlobalVector& b,
                   GlobalMatrix* jac) const override;

private:
    //! Auxiliary data used by the local assemblers.
    PythonStData _source_term_data;

    //! Local assemblers for all elements of the source term mesh.
    std::vector<std::unique_ptr<PythonSourceTermLocalAssemblerInterface>>
        _local_assemblers;

    //! Whether or not to flush standard output before and after each call to
    //! Python code. Ensures right order of output messages and therefore
    //! simplifies debugging.
    bool const _flush_stdout;
};

}  // namespace SourceTerms::Python
}  // namespace ProcessLib
