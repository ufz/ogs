/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ASSEMBLERLIB_LOCALASSEMBLERBUILDER_H_
#define ASSEMBLERLIB_LOCALASSEMBLERBUILDER_H_

#include "LocalToGlobalIndexMap.h"

namespace AssemblerLib
{

/// The class calls a builder to initialize local assemblers for each mesh
/// element item MeshItem passing the latter as an argument to the builder
/// object.
/// The indices in global objects are provided by the LocalToGlobalIndexMap in
/// the construction.
///
/// \tparam MeshItem is usually a mesh element, but could be a mesh's edge or
/// node depending on discretization method.
/// \tparam Builder is the local assembler builder type, used for initialisation
/// of the local assembler data depending on the passed mesh item.
template<
    typename MeshItem,
    typename Builder>
class LocalAssemblerBuilder
{
public:
    /// \param builder is the local assembler builder for initialization of the
    /// local assembler data.
    /// \param local_to_global_index_map is used to determine the size of local
    /// matrix/vector for a mesh item.
    LocalAssemblerBuilder(Builder &builder,
        LocalToGlobalIndexMap const& local_to_global_index_map)
    : _builder(builder), _local_to_global_index_map(local_to_global_index_map)
    {}

    ~LocalAssemblerBuilder() {}

    /// Executes given builder for the given mesh item and a data item.
    /// The item_data is initialized depending on the mesh item and the passed
    /// arguments.
    ///
    /// The positions in the global matrix/vector are taken from
    /// the LocalToGlobalIndexMap provided in the constructor at index \c id.
    /// \attention The index \c id is not necesserily the mesh item's id.
    template <typename ItemData, typename ...Args>
    void operator()(
            std::size_t const id,
            MeshItem const* item,
            ItemData& item_data,
            Args&&... args) const
    {
        assert(_local_to_global_index_map.size() > id);

        LocalToGlobalIndexMap::RowColumnIndices const& indices =
            _local_to_global_index_map[id];

        assert(indices.rows.size() >= indices.columns.size());
        _builder(*item, item_data, indices.rows.size(),
            std::forward<Args>(args)...);
    }

protected:
    Builder &_builder;
    LocalToGlobalIndexMap const& _local_to_global_index_map;
};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_LOCALASSEMBLERBUILDER_H_
