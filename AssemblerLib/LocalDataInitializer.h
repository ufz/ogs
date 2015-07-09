/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
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

#include "MeshLib/Elements/Elements.h"

#include "NumLib/Fem/Integration/GaussIntegrationPolicy.h"

#include "NumLib/Fem/ShapeFunction/ShapeHex20.h"
#include "NumLib/Fem/ShapeFunction/ShapeHex8.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine2.h"
#include "NumLib/Fem/ShapeFunction/ShapeLine3.h"
#include "NumLib/Fem/ShapeFunction/ShapePrism15.h"
#include "NumLib/Fem/ShapeFunction/ShapePrism6.h"
#include "NumLib/Fem/ShapeFunction/ShapePyra13.h"
#include "NumLib/Fem/ShapeFunction/ShapePyra5.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad4.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad8.h"
#include "NumLib/Fem/ShapeFunction/ShapeQuad9.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet10.h"
#include "NumLib/Fem/ShapeFunction/ShapeTet4.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri3.h"
#include "NumLib/Fem/ShapeFunction/ShapeTri6.h"

namespace AssemblerLib
{

/// The LocalDataInitializer is a functor creating a local assembler data with
/// corresponding to the mesh element type shape functions and calling
/// initialization of the new local assembler data.
/// For example for MeshLib::Line a local assembler data with template argument
/// NumLib::ShapeLine2 is created.
template <
    template <typename, typename> class LocalAssemblerDataInterface_,
    template <typename, typename, typename, typename, unsigned> class LocalAssemblerData_,
    typename GlobalMatrix_,
    typename GlobalVector_,
    unsigned GlobalDim>
class LocalDataInitializer
{
    template <typename ShapeFunction_>
    using IntegrationMethod = typename NumLib::GaussIntegrationPolicy<
                typename ShapeFunction_::MeshElement>::IntegrationMethod;

    template <typename ShapeFunction_>
        using LAData = LocalAssemblerData_<
                ShapeFunction_,
                IntegrationMethod<ShapeFunction_>,
                GlobalMatrix_, GlobalVector_, GlobalDim>;


public:
    LocalDataInitializer()
    {
        _builder[std::type_index(typeid(MeshLib::Hex20))] =
            [](){ return new LAData<NumLib::ShapeHex20>; };
        _builder[std::type_index(typeid(MeshLib::Hex))] =
            [](){ return new LAData<NumLib::ShapeHex8>; };
        _builder[std::type_index(typeid(MeshLib::Line))] =
            [](){ return new LAData<NumLib::ShapeLine2>; };
        _builder[std::type_index(typeid(MeshLib::Line3))] =
            [](){ return new LAData<NumLib::ShapeLine3>; };
        _builder[std::type_index(typeid(MeshLib::Prism15))] =
            [](){ return new LAData<NumLib::ShapePrism15>; };
        _builder[std::type_index(typeid(MeshLib::Prism))] =
            [](){ return new LAData<NumLib::ShapePrism6>; };
        _builder[std::type_index(typeid(MeshLib::Pyramid13))] =
            [](){ return new LAData<NumLib::ShapePyra13>; };
        _builder[std::type_index(typeid(MeshLib::Pyramid))] =
            [](){ return new LAData<NumLib::ShapePyra5>; };
        _builder[std::type_index(typeid(MeshLib::Quad))] =
            [](){ return new LAData<NumLib::ShapeQuad4>; };
        _builder[std::type_index(typeid(MeshLib::Quad8))] =
            [](){ return new LAData<NumLib::ShapeQuad8>; };
        _builder[std::type_index(typeid(MeshLib::Quad9))] =
            [](){ return new LAData<NumLib::ShapeQuad9>; };
        _builder[std::type_index(typeid(MeshLib::Tet10))] =
            [](){ return new LAData<NumLib::ShapeTet10>; };
        _builder[std::type_index(typeid(MeshLib::Tet))] =
            [](){ return new LAData<NumLib::ShapeTet4>; };
        _builder[std::type_index(typeid(MeshLib::Tri))] =
            [](){ return new LAData<NumLib::ShapeTri3>; };
        _builder[std::type_index(typeid(MeshLib::Tri6))] =
            [](){ return new LAData<NumLib::ShapeTri6>; };
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
