/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/IndexValueVector.h"
#include "ProcessLib/BoundaryCondition/BoundaryCondition.h"
#include "ProcessLib/BoundaryCondition/GenericNaturalBoundaryConditionLocalAssembler.h"

#include "BHEInflowPythonBoundaryConditionPythonSideInterface.h"

#include "ProcessLib/HeatTransportBHE/BHE/BHEAbstract.h"
#include "ProcessLib/Parameter/Parameter.h"

namespace ProcessLib
{
//! Groups data used by essential and natural BCs, in particular by the
//! local assemblers of the latter.

//! A boundary condition whose values are computed by a Python script.
class BHEInflowPythonBoundaryCondition final : public BoundaryCondition
{
public:
    BHEInflowPythonBoundaryCondition(std::pair<GlobalIndexType, GlobalIndexType>&& in_out_global_indices,
                            MeshLib::Mesh const& bc_mesh,
                            std::vector<MeshLib::Node*> const& vec_inflow_bc_nodes,
                            int const variable_id,
                            int const component_id,
                            std::unique_ptr<ProcessLib::HeatTransportBHE::BHE::BHEAbstract> const&
                            pt_bhe,
                            BHEInflowPythonBoundaryConditionPythonSideInterface* py_bc_object);

    void getEssentialBCValues(
        const double t, const GlobalVector& x,
        NumLib::IndexValueVector<GlobalIndexType>& bc_values) const override;

    void applyNaturalBC(const double t, const GlobalVector& x, GlobalMatrix& K,
                        GlobalVector& b, GlobalMatrix* Jac) override;

    void preTimestep(const double t, const GlobalVector& x) override;

private:

    // TODO (haibing) re-organize as the bottom BC data structure
    MeshLib::Mesh const& _bc_mesh;

    /// Stores the results of the outflow temperatures per boundary node.
    std::vector<double> _T_out_values;
    std::vector<GlobalIndexType> _T_out_indices;

    NumLib::IndexValueVector<GlobalIndexType> _bc_values;

    std::unique_ptr<ProcessLib::HeatTransportBHE::BHE::BHEAbstract> const&
        _pt_bhe;

    BHEInflowPythonBoundaryConditionPythonSideInterface* _py_bc_object;
};

//! Creates a new PythonBoundaryCondition object.
std::unique_ptr<BHEInflowPythonBoundaryCondition>createBHEInflowPythonBoundaryCondition(
    std::pair<GlobalIndexType, GlobalIndexType>&& in_out_global_indices,
    MeshLib::Mesh const& bc_mesh,
    std::vector<MeshLib::Node*> const& vec_outflow_bc_nodes,
    int const variable_id, int const component_id,
    std::unique_ptr<ProcessLib::HeatTransportBHE::BHE::BHEAbstract> const&
        pt_bhe,
    BHEInflowPythonBoundaryConditionPythonSideInterface* py_bc_object);

}  // namespace ProcessLib
