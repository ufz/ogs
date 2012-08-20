/**
 * \file DirectConditionGenerator.h
 * 2012/01/04 KR Initial implementation
 *
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

	const std::vector< std::pair<size_t,double> >& directToSurfaceNodes(const MeshLib::Mesh &mesh, const std::string &filename);

	//TODO6
	//const std::vector< std::pair<size_t,double> >& directWithSurfaceIntegration(MeshLib::Mesh &mesh, const std::string &filename, double scaling);

	int writeToFile(const std::string &name) const;

private:
	std::vector< std::pair<size_t,double> > _direct_values;

};

#endif // DIRECTCONDITIONGENERATOR_H

