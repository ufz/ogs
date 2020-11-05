/**
 * \file
 * \author Tobias Meisel
 * \date   2020-10-06
 * \brief  Writes MeshLib:Mesh to a Xdmf formated file
 *
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <filesystem>

namespace MeshLib
{
class Mesh;

namespace IO
{
/// Writes mesh to XDMF file.
/// \return True on success, false on error
/// \param mesh           Mesh holds all data to be written.
/// \param file_name      File name.
bool writeXdmf3(MeshLib::Mesh const& mesh, std::filesystem::path const& file_path);

}  // end namespace IO
}  // end namespace MeshLib