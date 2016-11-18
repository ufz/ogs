/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "MeshLib/Mesh.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/NumericsConfig.h"
#include "NumLib/NamedFunctionCaller.h"
#include "SecondaryVariableContext.h"

namespace ProcessLib
{
//! Computes a global vector from a NamedFunction (which can only compute double
//! values).
class GlobalVectorFromNamedFunction final
{
public:
    /*! Constructs a new instance.
     *
     * \param function_caller will provide the individual entries of the
     * GlobalVector to be computed.
     * \param mesh to which the \c dof_table_single is associated
     * \param dof_table_single used for constructing the GlobalVector
     * \param context used by the \c function_caller to access "auxiliary
     * unbound arguments".
     */
    GlobalVectorFromNamedFunction(
        NumLib::SpecificFunctionCaller&& function_caller,
        MeshLib::Mesh const& mesh,
        NumLib::LocalToGlobalIndexMap const& dof_table_single,
        SecondaryVariableContext& context);

    //! Computes the GlobalVector.
    //!
    //! The signature of this method matches
    //! SecondaryVariableFunctions::Function, i.e., this method can be used to
    //! compute a secondary variable.
    GlobalVector const& call(const double t, GlobalVector const& x,
                             NumLib::LocalToGlobalIndexMap const& dof_table,
                             std::unique_ptr<GlobalVector>& result);

private:
    NumLib::SpecificFunctionCaller _function_caller;
    MeshLib::Mesh const& _mesh;
    NumLib::LocalToGlobalIndexMap const& _dof_table_single;
    SecondaryVariableContext& _context;
};
}  // namespace ProcessLib
