/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file VTKInterface.h
 *
 *  Created on 2012-08-30 by Karsten Rink
 */

#ifndef VTKINTERFACE_H_
#define VTKINTERFACE_H_

#include "Writer.h"

#include <string>
#include <vector>

namespace MeshLib {
	class Mesh;
	class Node;
	class Element;
}

namespace FileIO
{

/**
 * \brief Reads and writes VTK-files to and from OGS data structures.
 */
class VTKInterface : public Writer
{
public:
	VTKInterface();
	~VTKInterface();

	static MeshLib::Mesh* readVTUFile(const std::string &file_name);

protected:
	int write(std::ostream& stream);

private:
	static MeshLib::Element* readElement(std::stringstream &iss, const std::vector<MeshLib::Node*> &nodes, unsigned material, unsigned type);
};
}

#endif /* VTKINTERFACE_H_ */
