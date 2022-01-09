/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Dense>
#include <map>
#include <memory>
#include <optional>
#include <utility>
#include <vector>

#include "BaseLib/Error.h"
#include "CoordinateSystem.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Node.h"
#include "SpatialPosition.h"

namespace BaseLib
{
class ConfigTree;
}  // namespace BaseLib

namespace MathLib
{
class PiecewiseLinearInterpolation;
}  // namespace MathLib

namespace MeshLib
{
class Mesh;
}  // namespace MeshLib

namespace ParameterLib
{
/// Base class for all parameters, not an interface class. This avoids using of
/// void* when storing parameters and convenient destruction.
/// Its property name helps addressing the right parameter.
struct ParameterBase
{
    explicit ParameterBase(std::string name_,
                           MeshLib::Mesh const* mesh = nullptr)
        : name(std::move(name_)), _mesh(mesh)
    {
    }

    virtual ~ParameterBase() = default;

    virtual bool isTimeDependent() const = 0;

    void setCoordinateSystem(CoordinateSystem const& coordinate_system)
    {
        _coordinate_system = coordinate_system;
    }

    /// Parameters might depend on each other; this method allows to set up the
    /// dependencies between parameters after they have been constructed.
    virtual void initialize(
        std::vector<std::unique_ptr<ParameterBase>> const& /*parameters*/)
    {
    }

    MeshLib::Mesh const* mesh() const { return _mesh; }

    std::string const name;

protected:
    std::vector<double> rotateWithCoordinateSystem(
        std::vector<double> const& values, SpatialPosition const& pos) const
    {
        assert(!!_coordinate_system);  // It is checked before calling this
                                       // function.

        // Don't rotate isotropic/scalar values.
        if (values.size() == 1)
        {
            return values;
        }
        if (values.size() == 2)
        {
            auto const result =
                _coordinate_system->rotateDiagonalTensor<2>(values, pos);
            return {result(0, 0), result(0, 1), result(1, 0), result(1, 1)};
        }
        if (values.size() == 3)
        {
            auto const result =
                _coordinate_system->rotateDiagonalTensor<3>(values, pos);
            return {
                result(0, 0), result(0, 1), result(0, 2),
                result(1, 0), result(1, 1), result(1, 2),
                result(2, 0), result(2, 1), result(2, 2),
            };
        }
        if (values.size() == 4)
        {
            auto const result =
                _coordinate_system->rotateTensor<2>(values, pos);
            return {result(0, 0), result(0, 1), result(1, 0), result(1, 1)};
        }
        if (values.size() == 9)
        {
            auto const result =
                _coordinate_system->rotateTensor<3>(values, pos);
            return {
                result(0, 0), result(0, 1), result(0, 2),
                result(1, 0), result(1, 1), result(1, 2),
                result(2, 0), result(2, 1), result(2, 2),
            };
        }
        OGS_FATAL(
            "Coordinate transformation for a {:d}-component parameter is not "
            "implemented.",
            values.size());
    }

protected:
    std::optional<CoordinateSystem> _coordinate_system;

    /// A mesh on which the parameter is defined. Some parameters might be
    /// mesh-independent.
    MeshLib::Mesh const* _mesh;
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
    using ParameterBase::ParameterBase;

    ~Parameter() override = default;

    //! Returns the number of components this Parameter has at every position
    //! and point in time.
    virtual int getNumberOfGlobalComponents() const = 0;

    //! Returns the parameter value at the given time and position.
    virtual std::vector<T> operator()(double const t,
                                      SpatialPosition const& pos) const = 0;

    //! Returns a matrix of values for all nodes of the given element.
    //
    // The matrix is of the shape NxC, where N is the number of nodes and C is
    // the number of components, such that subsequent multiplication with shape
    // functions matrix from left (which is a row vector) results in a row
    // vector of length C.
    //
    // The default implementation covers all cases, but the derived classes may
    // provide faster implementations.
    virtual Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>
    getNodalValuesOnElement(MeshLib::Element const& element,
                            double const t) const
    {
        auto const n_nodes = static_cast<int>(element.getNumberOfNodes());
        auto const n_components = getNumberOfGlobalComponents();
        Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> result(n_nodes,
                                                                n_components);

        // Column vector of values, copied for each node.
        SpatialPosition x_position;
        auto const nodes = element.getNodes();
        for (int i = 0; i < n_nodes; ++i)
        {
            x_position.setAll(
                nodes[i]->getID(), element.getID(), std::nullopt, *nodes[i]);
            auto const& values = this->operator()(t, x_position);
            auto const row_values =
                Eigen::Map<Eigen::Matrix<T, Eigen::Dynamic, 1> const>(
                    values.data(), values.size());
            result.row(i) = row_values;
        }

        return result;
    }
};

//! Constructs a new ParameterBase from the given configuration.
//!
//! The \c meshes vector is used to set up parameters from mesh input data.
std::unique_ptr<ParameterBase> createParameter(
    BaseLib::ConfigTree const& config,
    const std::vector<std::unique_ptr<MeshLib::Mesh>>& meshes,
    std::map<std::string,
             std::unique_ptr<MathLib::PiecewiseLinearInterpolation>> const&
        curves);

//! Checks whether the parameter can be used on the given mesh. The parameter's
//! domain of definition can be arbitrary (like for constant parameters), or the
//! parameter is defined on a mesh. In the latter case that mesh must be equal
//! to the given mesh.
//! \returns nothing if the parameter can be used on the given mesh, or an error
//! string otherwise.
std::optional<std::string> isDefinedOnSameMesh(ParameterBase const& parameter,
                                               MeshLib::Mesh const& mesh);

}  // namespace ParameterLib
