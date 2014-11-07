/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef ASSEMBLER_LIB_LOCALDATAINITIALIZER_H_
#define ASSEMBLER_LIB_LOCALDATAINITIALIZER_H_

#include <functional>
#include <typeindex>
#include <typeinfo>
#include <unordered_map>

#include "MeshLib/Elements/Line.h"
#include "MeshLib/Elements/Quad.h"
#include "MeshLib/Elements/Hex.h"

#include "NumLib/Fem/ShapeFunction/ShapeLine2.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri3.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex8.h"

namespace AssemblerLib
{

/// The LocalDataInitializer is a functor creating a local assembler data with
/// corresponding to the mesh element type shape functions and calling
/// initialization of the new local assembler data.
/// For example for MeshLib::Line a local assembler data with template argument
/// NumLib::ShapeLine2 is created.
template <
    template <typename, typename> class LocalAssemblerDataInterface_,
    template <typename, typename, typename, typename> class LocalAssemblerData_,
    template <typename> class IntegrationPolicy_,
    typename GlobalMatrix_,
    typename GlobalVector_>
class LocalDataInitializer
{
    template <typename ShapeFunction_>
        using LAData = LocalAssemblerData_<
                ShapeFunction_, IntegrationPolicy_<ShapeFunction_>,
                GlobalMatrix_, GlobalVector_>;


public:
    LocalDataInitializer()
    {
        _builder[std::type_index(typeid(MeshLib::Line))] =
                [](){ return new LAData<NumLib::ShapeLine2>; };
        _builder[std::type_index(typeid(MeshLib::Tri))] =
                [](){ return new LAData<NumLib::ShapeTri3>; };
        _builder[std::type_index(typeid(MeshLib::Quad))] =
                [](){ return new LAData<NumLib::ShapeQuad4>; };
        _builder[std::type_index(typeid(MeshLib::Hex))] =
                [](){ return new LAData<NumLib::ShapeHex8>; };
    }

    /// Sets the provided data_ptr to the newly created local assembler data and
    /// calls init() forwarding all remaining arguments.
    template <typename ...Args_>
    void operator()(const MeshLib::Element& e,
        LocalAssemblerDataInterface_<GlobalMatrix_, GlobalVector_>*& data_ptr, Args_&&... args)
    {
        data_ptr = _builder[std::type_index(typeid(e))]();
        data_ptr->init(e, std::forward<Args_>(args)...);
    }

private:
    /// Mapping of element types to local assembler constructors.
    std::unordered_map<
        std::type_index,
        std::function<LocalAssemblerDataInterface_<GlobalMatrix_, GlobalVector_>*()>
            > _builder;
};

}   // namespace AssemblerLib

#endif  // ASSEMBLER_LIB_LOCALDATAINITIALIZER_H_
