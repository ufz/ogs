/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_LOCALTOGLOBALINDEXMAP_H_
#define NUMLIB_LOCALTOGLOBALINDEXMAP_H_

#ifndef NDEBUG
#include <iosfwd>
#endif  // NDEBUG

#include <vector>

#include <Eigen/Dense>

#include "MathLib/LinAlg/RowColumnIndices.h"
#include "MeshLib/MeshSubsets.h"

#include "MeshComponentMap.h"

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
    typedef MathLib::RowColumnIndices<GlobalIndexType> RowColumnIndices;
    typedef RowColumnIndices::LineIndex LineIndex;

public:
    /// Creates a MeshComponentMap internally and stores the global indices for
    /// each mesh element of the given mesh_subsets.
    explicit LocalToGlobalIndexMap(
        std::vector<std::unique_ptr<MeshLib::MeshSubsets>>&& mesh_subsets,
        NumLib::ComponentOrder const order);

    /// Derive a LocalToGlobalIndexMap constrained to a set of mesh subsets and
    /// elements. A new mesh component map will be constructed using the passed
    /// mesh_subsets for the given variable and component ids.
    ///
    /// \note The elements are not necessarily those used in the mesh_subsets.
    LocalToGlobalIndexMap* deriveBoundaryConstrainedMap(
        int const variable_id,
        int const component_id,
        std::unique_ptr<MeshLib::MeshSubsets>&& mesh_subsets,
        std::vector<MeshLib::Element*> const& elements) const;

    /// Returns total number of degrees of freedom including those located in
    /// the ghost nodes.
    std::size_t dofSizeWithGhosts() const;

    /// Returns total number of local degrees of freedom of the present rank,
    /// which does not count the unknowns associated with ghost nodes (for DDC
    /// with node-wise mesh partitioning).
    std::size_t dofSizeWithoutGhosts() const
    {
        return _mesh_component_map.dofSizeWithoutGhosts();
    }

    std::size_t size() const;

    std::size_t getNumberOfComponents() const { return _mesh_subsets.size(); }

    RowColumnIndices operator()(std::size_t const mesh_item_id, const unsigned component_id) const;

    std::size_t getNumberOfElementDOF(std::size_t const mesh_item_id) const;

    GlobalIndexType getGlobalIndex(MeshLib::Location const& l,
                                   int const variable_id,
                                   int const component_id) const
    {
        auto const c = getGlobalComponent(variable_id, component_id);
        return _mesh_component_map.getGlobalIndex(l, c);
    }

    /// Forwards the respective method from MeshComponentMap.
    std::vector<GlobalIndexType> getGlobalIndices(const MeshLib::Location &l) const
    {
        return _mesh_component_map.getGlobalIndices(l);
    }

    /// Get ghost indices, forwarded from MeshComponentMap.
    std::vector<GlobalIndexType> const& getGhostIndices() const
    {
        return _mesh_component_map.getGhostIndices();
    }

    /// Computes the index in a local (for DDC) vector for a given location and
    /// component; forwarded from MeshComponentMap.
    GlobalIndexType getLocalIndex(MeshLib::Location const& l, std::size_t const comp_id,
                                  std::size_t const range_begin,
                                  std::size_t const range_end) const
    {
        return _mesh_component_map.getLocalIndex(l, comp_id, range_begin,
                                                 range_end);
    }

    MeshLib::MeshSubsets const& getMeshSubsets(int const variable_id,
                                               int const component_id) const
    {
        return *_mesh_subsets[getGlobalComponent(variable_id, component_id)];
    }

private:
    /// Private constructor used by internally created local-to-global index
    /// maps. The mesh_component_map is passed as argument instead of being
    /// created by the constructor.
    /// \attention The passed mesh_component_map is in undefined state after
    /// this construtor.
    explicit LocalToGlobalIndexMap(
        std::vector<std::unique_ptr<MeshLib::MeshSubsets>>&& mesh_subsets,
        int const component_id,
        std::vector<MeshLib::Element*> const& elements,
        NumLib::MeshComponentMap&& mesh_component_map);

    template <typename ElementIterator>
    void
    findGlobalIndices(ElementIterator first, ElementIterator last,
        std::vector<MeshLib::Node*> const& nodes,
        std::size_t const mesh_id,
        const unsigned component_id, const unsigned comp_id_write);

    /// The global component id for the specific variable (like velocity) and a
    /// component (like x, or y, or z).
    std::size_t getGlobalComponent(int const variable_id,
                                   int const component_id) const
    {
        return _variable_component_offsets[variable_id] + component_id;
    }

private:
    /// A vector of mesh subsets for each process variables' components.
    std::vector<std::unique_ptr<MeshLib::MeshSubsets>> const _mesh_subsets;
    NumLib::MeshComponentMap _mesh_component_map;

    using Table = Eigen::Matrix<LineIndex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

    /// Table contains for each element (first index) and each component (second index)
    /// a vector (\c LineIndex) of indices in the global stiffness matrix or vector
    Table _rows;

    /// \see _rows
    Table const& _columns = _rows;

    std::vector<int> _variable_component_offsets;
#ifndef NDEBUG
    /// Prints first rows of the table, every line, and the mesh component map.
    friend std::ostream& operator<<(std::ostream& os, LocalToGlobalIndexMap const& map);
#endif  // NDEBUG

};

}   // namespace NumLib

#endif  // NUMLIB_LOCALTOGLOBALINDEXMAP_H_
