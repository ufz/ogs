/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#ifndef PROCESSLIB_SMALLDEFORMATION_WITH_LIE_COMMON_LEVELSETFUNCTION_H_
#define PROCESSLIB_SMALLDEFORMATION_WITH_LIE_COMMON_LEVELSETFUNCTION_H_

namespace ProcessLib
{
namespace SmallDeformationWithLIE
{

struct FractureProperty;

/// calculate the level set function
/// \f$ \psi(\mathbf{x}) = H(|\mathbf{x} - \mathbf{x_d}| \mathrm{sign}[\mathbf{n_d} \cdot (\mathbf{x}-\mathbf{x_d}]) \f$
/// where \f$H(u)\f$ is the Heaviside step function, \f$\mathbf{x_d}\f$ is a point on the fracture plane,
/// and \f$\mathbf{n_d}\f$ is the normal vector of a fracture plane
double calculateLevelSetFunction(
        FractureProperty const& fracture_property, double const* x);

} // SmallDeformationWithLIE
} // ProcessLib

#endif // PROCESSLIB_SMALLDEFORMATION_WITH_LIE_COMMON_LEVELSETFUNCTION_H_
