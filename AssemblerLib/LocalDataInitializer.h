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

template <
    template <typename, typename> class LocalAssemblerDataBase_,
    template <typename, typename, unsigned, typename, typename> class LocalAssemblerData_,
    template <typename> class IntegrationPolicy_,
    unsigned IntegrationOrder_,
    typename GlobalMatrix,
    typename GlobalVector>
class LocalDataInitializer
{
    template <typename ShapeFunction_>
        using LAData = LocalAssemblerData_<
                ShapeFunction_, IntegrationPolicy_<ShapeFunction_>, IntegrationOrder_,
                GlobalMatrix, GlobalVector>;


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

    template <typename ...Args_>
    void operator()(const MeshLib::Element& e,
        LocalAssemblerDataBase_<GlobalMatrix, GlobalVector>*& data_ptr, Args_&&... args)
    {
        data_ptr = _builder[std::type_index(typeid(e))]();
        data_ptr->init(e, std::forward<Args_>(args)...);
    }

private:
    /// Mapping of element types to local assembler constructors.
    std::unordered_map<
        std::type_index,
        std::function<LocalAssemblerDataBase_<GlobalMatrix, GlobalVector>*()>
            > _builder;
};

}   // namespace AssemblerLib

#endif  // ASSEMBLER_LIB_LOCALDATAINITIALIZER_H_
