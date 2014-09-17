/**
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_PROCESS_H_
#define PROCESS_LIB_PROCESS_H_

#include "MeshLib/Mesh.h"

namespace ProcessLib
{

class Process
{
public:
    Process(MeshLib::Mesh const& mesh, std::size_t const numberOfNodeComponents)
        : _mesh(mesh), _numberOfNodeComponents(numberOfNodeComponents)
    { }

    virtual ~Process() = default;

    virtual void initialize() = 0;

protected:
    MeshLib::Mesh const& _mesh;
    std::size_t const _numberOfNodeComponents;
};

}   // namespace ProcessLib

//
// Include all known processes here.
//
#include "GroundwaterFlow.h"

#endif  // PROCESS_LIB_PROCESS_H_
