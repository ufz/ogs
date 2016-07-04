/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_NEUMANN_BC_H_
#define PROCESS_LIB_NEUMANN_BC_H_

#include <memory>
#include <vector>

#include "MeshLib/MeshSubset.h"
#include "NumLib/NumericsConfig.h"
#include "NeumannBcConfig.h"
#include "NeumannBcAssembler.h"

namespace ProcessLib
{

/// The NeumannBc class is a process which is integrating a single, Neumann type
/// boundary condition in to the global matrix and the right-hand-side.
///
/// The process operates on a set of elements and a subset of the DOF-table (the
/// local to global index map). For each element a local assembler is created
/// (NeumannBcAssembler).
///
/// In the construction phase the NeumannBcConfig together with global DOF-table
/// and mesh subset are used.
/// The creation of the local assemblers and binding to the global matrix and
/// right-hand-sides happen in the initialize() function.
/// The integration() function provides calls then the actual integration of the
/// Neumann boundary condition.
class NeumannBc
{
public:
    /// Create a Neumann boundary condition process from given config,
    /// DOF-table, and a mesh subset for a given variable and its component.
    /// A local DOF-table, a subset of the given one, is constructed.
    NeumannBc(
        NeumannBcConfig const& bc,
        unsigned const integration_order,
        NumLib::LocalToGlobalIndexMap const& local_to_global_index_map,
        int const variable_id,
        int const component_id);

    ~NeumannBc();

    /// Calls local assemblers which calculate their contributions to the global
    /// matrix and the right-hand-side.
    void integrate(const double t, GlobalVector& b);

    void initialize(unsigned global_dim);

private:
    /// The right-hand-side function of the Neumann boundary condition given as
    /// \f$ \alpha(x) \, \partial u(x) / \partial n = \text{_function}(x)\f$.
    MathLib::ConstantFunction<double> const _function;

    /// Vector of lower-dimensional elements on which the boundary condition is
    /// defined.
    std::vector<MeshLib::Element*> _elements;

    MeshLib::MeshSubset const* _mesh_subset_all_nodes = nullptr;

    /// Local dof table, a subset of the global one restricted to the
    /// participating #_elements of the boundary condition.
    std::unique_ptr<NumLib::LocalToGlobalIndexMap> _local_to_global_index_map;

    /// Integration order for integration over the lower-dimensional elements of
    /// the #_function.
    unsigned const _integration_order;

    /// Local assemblers for each element of #_elements.
    std::vector<std::unique_ptr<LocalNeumannBcAsmDataInterface>> _local_assemblers;

};

}   // namespace ProcessLib

#endif  // PROCESS_LIB_NEUMANN_BC_H_
