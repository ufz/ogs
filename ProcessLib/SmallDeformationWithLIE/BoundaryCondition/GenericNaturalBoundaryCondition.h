/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_SMALLDEFORMATION_WITH_LIE_GENERICNATURALBOUNDARYCONDITION_H_
#define PROCESSLIB_SMALLDEFORMATION_WITH_LIE_GENERICNATURALBOUNDARYCONDITION_H_

#include "MeshLib/MeshSubset.h"
#include "ProcessLib/BoundaryCondition/BoundaryCondition.h"
#include "ProcessLib/SmallDeformationWithLIE/Common/FractureProperty.h"

namespace ProcessLib
{
class GenericNaturalBoundaryConditionLocalAssemblerInterface;

namespace SmallDeformationWithLIE
{

template <typename BoundaryConditionData,
          template <typename, typename, unsigned>
          class LocalAssemblerImplementation>
class GenericNaturalBoundaryCondition final : public BoundaryCondition
{
public:
    /// Create a boundary condition process from given config,
    /// DOF-table, and a mesh subset for a given variable and its component.
    /// A local DOF-table, a subset of the given one, is constructed.
    template <typename Data>
    GenericNaturalBoundaryCondition(
        typename std::enable_if<
            std::is_same<typename std::decay<BoundaryConditionData>::type,
                         typename std::decay<Data>::type>::value,
            bool>::type is_axially_symmetric,
        unsigned const integration_order,
        NumLib::LocalToGlobalIndexMap const& dof_table_bulk,
        int const variable_id, int const component_id,
        unsigned const global_dim,
        std::vector<MeshLib::Element*>&& elements, Data&& data,
        FractureProperty const& fracture_prop);

    ~GenericNaturalBoundaryCondition() override;

    /// Calls local assemblers which calculate their contributions to the global
    /// matrix and the right-hand-side.
    void applyNaturalBC(const double t,
                        GlobalVector const& x,
                        GlobalMatrix& K,
                        GlobalVector& b) override;

private:
    /// Data used in the assembly of the specific boundary condition.
    BoundaryConditionData _data;

    /// Vector of lower-dimensional elements on which the boundary condition is
    /// defined.
    std::vector<MeshLib::Element*> _elements;

    std::unique_ptr<MeshLib::MeshSubset const> _mesh_subset_all_nodes;

    /// Local dof table, a subset of the global one restricted to the
    /// participating #_elements of the boundary condition.
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> _dof_table_boundary;

    /// Integration order for integration over the lower-dimensional elements
    unsigned const _integration_order;

    /// Local assemblers for each element of #_elements.
    std::vector<
        std::unique_ptr<GenericNaturalBoundaryConditionLocalAssemblerInterface>>
        _local_assemblers;
};

}  // SmallDeformationWithLIE
}  // ProcessLib

#include "GenericNaturalBoundaryCondition-impl.h"

#endif // PROCESSLIB_SMALLDEFORMATION_WITH_LIE_GENERICNATURALBOUNDARYCONDITION_H_
