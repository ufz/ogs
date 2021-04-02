/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "GenerateFaceSetMeshes.h"

#include "MeshLib/Elements/Quad.h"
#include "MeshLib/IO/writeMeshToFile.h"

namespace FileIO
{
namespace Gocad
{
void generateFaceSets(GocadSGridReader const& reader, std::string const& path)
{
    for (std::size_t l(0); l < 128; l++)
    {
        std::unique_ptr<MeshLib::Mesh> face_set(reader.getFaceSetMesh(l));

        if (!face_set)
        {
            continue;
        }
        INFO("Face set mesh created. #nodes: {:d}, #elements: {:d}",
             face_set->getNumberOfNodes(),
             face_set->getNumberOfElements());

        std::string const mesh_out_fname(path + face_set->getName() + ".vtu");
        INFO("Writing face set mesh to '{:s}'.", mesh_out_fname);
        MeshLib::IO::writeMeshToFile(*face_set, mesh_out_fname);
    }
}

}  //  namespace Gocad
}  //  namespace FileIO
