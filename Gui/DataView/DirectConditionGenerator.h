/**
 * Copyright (c) 2012, OpenGeoSys Community (http://www.opengeosys.net)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.net/LICENSE.txt
 *
 * \file DirectConditionGenerator.h
 *
 * Created on 2012-01-04 by Karsten Rink
 */

#ifndef DIRECTCONDITIONGENERATOR_H
#define DIRECTCONDITIONGENERATOR_H

#include <vector>

namespace MeshLib {
	class Mesh;
}

class DirectConditionGenerator
{
public:
	DirectConditionGenerator() {};
	~DirectConditionGenerator() {};

	const std::vector< std::pair<std::size_t,double> >& directToSurfaceNodes(const MeshLib::Mesh &mesh, const std::string &filename);

	const std::vector< std::pair<std::size_t,double> >& directWithSurfaceIntegration(MeshLib::Mesh &mesh, const std::string &filename, double scaling);

	int writeToFile(const std::string &name) const;

private:
	std::vector< std::pair<std::size_t,double> > _direct_values;

};

#endif // DIRECTCONDITIONGENERATOR_H

