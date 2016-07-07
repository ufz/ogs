/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_GLOBALVECTORFROMNAMEDFUNCTION_H
#define PROCESSLIB_GLOBALVECTORFROMNAMEDFUNCTION_H

#include "MeshLib/Mesh.h"
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/NumericsConfig.h"
#include "NumLib/NamedFunctionCaller.h"
#include "SecondaryVariableContext.h"

namespace ProcessLib
{
class GlobalVectorFromNamedFunction
{
public:
    GlobalVectorFromNamedFunction(
        NumLib::SpecialFunctionCaller&& function_caller,
        MeshLib::Mesh const& mesh,
        NumLib::LocalToGlobalIndexMap const& dof_table_single,
        SecondaryVariableContext& context);

    GlobalVector const& call(GlobalVector const& x,
                             NumLib::LocalToGlobalIndexMap const& dof_table,
                             std::unique_ptr<GlobalVector>& result_cache);

private:
    NumLib::SpecialFunctionCaller _function_caller;
    MeshLib::Mesh const& _mesh;
    NumLib::LocalToGlobalIndexMap const& _dof_table_single;
    SecondaryVariableContext& _context;
};
}  // namespace ProcessLib

#endif  // PROCESSLIB_GLOBALVECTORFROMNAMEDFUNCTION_H
