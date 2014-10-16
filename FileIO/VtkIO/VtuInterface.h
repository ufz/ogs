/**
 * \file
 * \author Lars Bilke
 * \date   2014-09-25
 * \brief  Implementation of the VtuInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2014, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef VTUINTERFACE_H_
#define VTUINTERFACE_H_

#include <string>

#include "Writer.h"

namespace MeshLib {
	class Mesh;
}

namespace FileIO
{

/**
 * \brief Reads and writes VtkXMLUnstructuredGrid-files (vtu) to and from OGS data structures.
 */
class VtuInterface : public Writer
{
public:
	VtuInterface();
	~VtuInterface();

	/// Read an unstructured grid from a VTU file
	static MeshLib::Mesh* readVTUFile(const std::string &file_name);

	/// Decide if the mesh data should be written compressed (default is false).
	void setCompressData(bool flag=true) { _use_compressor = flag; }

	/// Sets the mesh to write.
	void setMesh(const MeshLib::Mesh* mesh);

private:
	bool write();

	MeshLib::Mesh* _mesh;
	bool _use_compressor;
};

}

#endif /* VTUINTERFACE_H_ */
