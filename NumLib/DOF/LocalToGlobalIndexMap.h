/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#ifndef NDEBUG
#include <iosfwd>
#endif  // NDEBUG

#include <vector>

#include <Eigen/Dense>

#include "MathLib/LinAlg/RowColumnIndices.h"

#include "MeshComponentMap.h"

namespace MeshLib
{
class Node;
class Element;
}  // namespace MeshLib

namespace NumLib
{

/// Row and column indices in global linear algebra objects for each mesh item.
///
/// The row and column indices in global matrix and rhs vector for example,
/// required for addition of local contributions from every mesh item (like node
/// or cell) to global objects.
///
/// The number of rows should be equal to the number of mesh items and the
/// number of columns should be equal to the number of the components on that
/// mesh item.
class LocalToGlobalIndexMap final
{
public:
    using RowColumnIndices = MathLib::RowColumnIndices<GlobalIndexType>;
    using LineIndex = RowColumnIndices::LineIndex;

public:
    /// Creates a MeshComponentMap internally and stores the global indices for
    /// each mesh element of the given mesh_subsets.
    ///
    /// \attention This constructor assumes the number of the given mesh subsets
    /// is equal to the number of variables, i.e. every variable has a single component.
    LocalToGlobalIndexMap(std::vector<MeshLib::MeshSubset>&& mesh_subsets,
                          NumLib::ComponentOrder const order);

    /// Creates a MeshComponentMap internally and stores the global indices for
    /// each mesh element of the given mesh_subsets.
    ///
    /// \param mesh_subsets  a vector of components
    /// \param vec_var_n_components  a vector of the number of variable components.
    /// The size of the vector should be equal to the number of variables. Sum of the entries
    /// should be equal to the size of the mesh_subsets.
    /// \param order  type of ordering values in a vector
    LocalToGlobalIndexMap(std::vector<MeshLib::MeshSubset>&& mesh_subsets,
                          std::vector<int> const& vec_var_n_components,
                          NumLib::ComponentOrder const order);

    /// Creates a MeshComponentMap internally and stores the global indices for
    /// the given mesh elements
    ///
    /// \param mesh_subsets  a vector of components
    /// \param vec_var_n_components  a vector of the number of variable components.
    /// The size of the vector should be equal to the number of variables. Sum of the entries
    /// should be equal to the size of the mesh_subsets.
    /// \param vec_var_elements  a vector of active mesh elements for each variable.
    /// \param order  type of ordering values in a vector
    LocalToGlobalIndexMap(
        std::vector<MeshLib::MeshSubset>&& mesh_subsets,
        std::vector<int> const& vec_var_n_components,
        std::vector<std::vector<MeshLib::Element*> const*> const&
            vec_var_elements,
        NumLib::ComponentOrder const order);

    /// Derive a LocalToGlobalIndexMap constrained to the mesh subset and mesh
    /// subset's elements. A new mesh component map will be constructed using
    /// the passed mesh_subset for the given variable and component ids.
    LocalToGlobalIndexMap* deriveBoundaryConstrainedMap(
        int const variable_id,
        std::vector<int> const& component_ids,
        MeshLib::MeshSubset&& mesh_subset) const;

    /// Returns total number of degrees of freedom including those located in
    /// the ghost nodes.
    std::size_t dofSizeWithGhosts() const;

    /// Returns total number of local degrees of freedom of the present rank,
    /// which does not count the unknowns associated with ghost nodes (for DDC
    /// with node-wise mesh partitioning).
    std::size_t dofSizeWithoutGhosts() const;

    std::size_t size() const;

    int getNumberOfVariables() const;

    int getNumberOfVariableComponents(int variable_id) const;

    int getNumberOfComponents() const;

    RowColumnIndices operator()(std::size_t const mesh_item_id,
                                const int component_id) const;

    std::size_t getNumberOfElementDOF(std::size_t const mesh_item_id) const;

    std::size_t getNumberOfElementComponents(std::size_t const mesh_item_id) const;

    std::vector<int> getElementVariableIDs(
        std::size_t const mesh_item_id) const;

    GlobalIndexType getGlobalIndex(MeshLib::Location const& l,
                                   int const variable_id,
                                   int const component_id) const;

    GlobalIndexType getGlobalIndex(MeshLib::Location const& l,
                                   int const global_component_id) const;

    /// Forwards the respective method from MeshComponentMap.
    std::vector<GlobalIndexType> getGlobalIndices(
        const MeshLib::Location& l) const;

    /// Get ghost indices, forwarded from MeshComponentMap.
    std::vector<GlobalIndexType> const& getGhostIndices() const;

    /// Computes the index in a local (for DDC) vector for a given location and
    /// component; forwarded from MeshComponentMap.
    GlobalIndexType getLocalIndex(MeshLib::Location const& l,
                                  std::size_t const comp_id,
                                  std::size_t const range_begin,
                                  std::size_t const range_end) const;

    MeshLib::MeshSubset const& getMeshSubset(int const variable_id,
                                             int const component_id) const;

    MeshLib::MeshSubset const& getMeshSubset(
        int const global_component_id) const;

private:
    /// Private constructor used by internally created local-to-global index
    /// maps. The mesh_component_map is passed as argument instead of being
    /// created by the constructor.
    /// \attention The passed mesh_component_map is in undefined state after
    /// this construtor.
    explicit LocalToGlobalIndexMap(
        std::vector<MeshLib::MeshSubset>&& mesh_subsets,
        std::vector<int> const& global_component_ids,
        std::vector<MeshLib::Element*> const& elements,
        NumLib::MeshComponentMap&& mesh_component_map);

    template <typename ElementIterator>
    void findGlobalIndices(ElementIterator first, ElementIterator last,
                           std::vector<MeshLib::Node*> const& nodes,
                           std::size_t const mesh_id, const int component_id,
                           const int comp_id_write);

    template <typename ElementIterator>
    void findGlobalIndicesWithElementID(
        ElementIterator first, ElementIterator last,
        std::vector<MeshLib::Node*> const& nodes, std::size_t const mesh_id,
        const int component_id, const int comp_id_write);

    /// The global component id for the specific variable (like velocity) and a
    /// component (like x, or y, or z).
    int getGlobalComponent(int const variable_id, int const component_id) const;

private:
    /// A vector of mesh subsets for each process variables' components.
    std::vector<MeshLib::MeshSubset> _mesh_subsets;
    NumLib::MeshComponentMap _mesh_component_map;

    using Table = Eigen::Matrix<LineIndex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    /// Table contains for each element (first index) and each component (second index)
    /// a vector (\c LineIndex) of indices in the global stiffness matrix or vector
    Table _rows;

    /// \see _rows
    Table const& _columns = _rows;

    std::vector<int> const _variable_component_offsets;
#ifndef NDEBUG
    /// Prints first rows of the table, every line, and the mesh component map.
    friend std::ostream& operator<<(std::ostream& os, LocalToGlobalIndexMap const& map);
#endif  // NDEBUG

};

}   // namespace NumLib
