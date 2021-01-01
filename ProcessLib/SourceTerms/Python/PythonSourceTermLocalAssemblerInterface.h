
/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

namespace ProcessLib
{
namespace SourceTerms
{
namespace Python
{

class PythonSourceTermLocalAssemblerInterface
{
public:
    virtual void assemble(
        std::size_t const source_term_element_id,
        NumLib::LocalToGlobalIndexMap const& source_term_dof_table,
        double const t, const GlobalVector& x, GlobalVector& b,
        GlobalMatrix* Jac) = 0;

    virtual ~PythonSourceTermLocalAssemblerInterface() = default;
};

}  // namespace Python
}  // namespace SourceTerms
}  // namespace ProcessLib
