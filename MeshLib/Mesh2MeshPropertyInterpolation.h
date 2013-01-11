/**
 * \file
 * \author Thomas Fischer
 * \date   Oct 12, 2012
 * \brief  Implementation of the Mesh2MeshPropertyInterpolation class.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef MESH2MESHPROPERTYINTERPOLATION_H_
#define MESH2MESHPROPERTYINTERPOLATION_H_

namespace MeshLib {

class Mesh;
/**
 * Class Mesh2MeshPropertyInterpolation transfers properties of
 * mesh elements of a (source) mesh to mesh elements of another
 * (destination) mesh deploying weighted interpolation. The two
 * meshes must have the same dimension.
 */
class Mesh2MeshPropertyInterpolation {
public:
	/**
	 * Constructor taking the source or input mesh and properties.
	 * @param source_mesh the mesh the given property information is
	 * assigned to.
	 * @param source_properties a vector containing property information
	 * assigned to the mesh. The number of entries in the property vector
	 * must be at least the number of different properties stored in mesh.
	 * For instance if mesh has \f$n\f$ (pairwise) different property
	 * indices the vector of properties must have \f$\ge n\f$ entries.
	 */
	Mesh2MeshPropertyInterpolation(Mesh const*const source_mesh, std::vector<double> const*const source_properties);
	virtual ~Mesh2MeshPropertyInterpolation();

	/**
	 * Calculates entries for the property vector and sets appropriate indices in the
	 * mesh elements.
	 * @param mesh the mesh the property information will be calculated and set via
	 * weighted interpolation
	 * @param properties at input a vector of length equal to the number of elements,
	 * 	at output interpolated property values
	 * @return true if the operation was successful, false on error
	 */
	bool setPropertiesForMesh(Mesh *mesh, std::vector<double>& properties) const;

private:
	/**
	 *
	 * @param dest_mesh
	 * @param dest_properties
	 */
	void interpolatePropertiesForMesh(Mesh *dest_mesh, std::vector<double>& dest_properties) const;
	/**
	 * Method interpolates the element wise given properties to the nodes of the element
	 * @param interpolated_node_properties the vector must have the same number of entries as
	 * the source mesh has number of nodes, the content of the particular entries will be overwritten
	 */
	void interpolateElementPropertiesToNodeProperties(std::vector<double> &interpolated_node_properties) const;

	Mesh const*const _src_mesh;
	std::vector<double> const*const _src_properties;
};

} // end namespace MeshLib

#endif /* MESH2MESHPROPERTYINTERPOLATION_H_ */
