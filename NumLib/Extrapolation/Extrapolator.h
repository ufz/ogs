/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Core>
#include <vector>

#include "ExtrapolatableElementCollection.h"
#include "NumLib/NumericsConfig.h"

namespace NumLib
{
class LocalToGlobalIndexMap;

//! Interface for classes that extrapolate integration point values to nodal
//! values.
class Extrapolator
{
public:
    //! Extrapolates the given \c property from the given local assemblers.
    virtual void extrapolate(
        const unsigned num_components,
        ExtrapolatableElementCollection const& extrapolatables,
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table) = 0;

    /*! Computes residuals from the extrapolation of the given \c property.
     *
     * The residuals are computed as element values.
     *
     * \pre extrapolate() must have been called before with the same arguments.
     */
    virtual void calculateResiduals(
        const unsigned num_components,
        ExtrapolatableElementCollection const& extrapolatables,
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table) = 0;

    //! Returns the extrapolated nodal values.
    //! \todo Maybe write directly to a MeshProperty.
    virtual GlobalVector const& getNodalValues() const = 0;

    //! Returns the extrapolation residuals.
    //! \todo Maybe write directly to a MeshProperty.
    virtual GlobalVector const& getElementResiduals() const = 0;

    virtual ~Extrapolator() = default;
};

}  // namespace NumLib
