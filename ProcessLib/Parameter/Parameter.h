/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_PARAMETER_H_
#define PROCESS_LIB_PARAMETER_H_

#include <memory>
#include <vector>
#include "SpatialPosition.h"

namespace BaseLib
{
class ConfigTree;
}  // BaseLib

namespace MeshLib
{
class Mesh;
}  // MeshLib

namespace ProcessLib
{
/// Base class for all parameters, not an interface class. This avoids using of
/// void* when storing parameters and convenient destruction.
/// Its property name helps addressing the right parameter.
struct ParameterBase
{
    virtual ~ParameterBase() = default;

    std::string name;
};

/*! A Parameter is a function \f$ (t, x) \mapsto f(t, x) \in T^n \f$.
 *
 * Where \f$ t \f$ is the time and \f$ x \f$ is the SpatialPosition.
 * \f$ n \f$ is the number of components of \f$f\f$'s results, i.e., 1 for a
 * scalar parameter and >1 for a vectorial or tensorial parameter.
 */
template <typename T>
struct Parameter : public ParameterBase
{
    virtual ~Parameter() = default;

    //! Returns the number of components this Parameter has at every position and
    //! point in time.
    virtual unsigned getNumberOfComponents() const = 0;

    //! Returns the parameter value at the given time and position.
    virtual std::vector<T> const& getTuple(
        double const t, SpatialPosition const& pos) const = 0;
};

//! Constructs a new ParameterBase from the given configuration.
//!
//! The \c meshes vector is used to set up parameters from mesh input data.
std::unique_ptr<ParameterBase> createParameter(
    BaseLib::ConfigTree const& config,
    const std::vector<MeshLib::Mesh*>& meshes);

}  // namespace ProcessLib

#endif  // PROCESS_LIB_PARAMETER_H_
