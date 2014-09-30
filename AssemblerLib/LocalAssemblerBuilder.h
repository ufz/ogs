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

namespace AssemblerLib
{

/// The class calls a builder to initialize local assemblers for each mesh
/// element item MeshItem passing the latter as an argument to the builder
/// object.
///
/// \tparam MeshItem is usually a mesh element, but could be a mesh's edge or
/// node depending on discretization method.
/// \tparam Builder is the local assembler builder, initializing the data
/// depending on the passed mesh item.
template<
    typename MeshItem,
    typename Builder>
class LocalAssemblerBuilder
{
public:
    LocalAssemblerBuilder(Builder &builder)
    : _builder(builder) {}

    ~LocalAssemblerBuilder() {}

    /// Executes given builder for the given mesh item and a data item.
    /// The item_data is initialized depending on the mesh item and the passed
    /// arguments.
    ///
    /// The index id is not used here and is only required for the
    /// SerialExecutor::execute() interface.
    template <typename ItemData, typename ...Args>
    void operator()(
            std::size_t const /*id*/,
            MeshItem const* item,
            ItemData& item_data,
            Args&&... args) const
    {
        _builder(*item, item_data, std::forward<Args>(args)...);
    }

protected:
    Builder &_builder;
};

}   // namespace AssemblerLib

#endif  // ASSEMBLERLIB_LOCALASSEMBLERBUILDER_H_
