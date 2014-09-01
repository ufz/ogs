/**
 * \copyright
 * Copyright (c) 2014, OpenGeoSys Community (http://www.opengeosys.org)
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
    Process(MeshLib::Mesh const& mesh)
        : _mesh(mesh)
    { }

    virtual ~Process() = default;

private:
    MeshLib::Mesh const& _mesh;
};

}   // namespace ProcessLib

//
// Include all known processes here.
//
#include "GroundwaterFlow.h"

#endif  // PROCESS_LIB_PROCESS_H_
