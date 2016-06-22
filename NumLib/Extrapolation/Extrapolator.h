/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_EXTRAPOLATOR_H
#define NUMLIB_EXTRAPOLATOR_H

#include <vector>

#include <Eigen/Eigen>
#include "NumLib/NumericsConfig.h"

namespace NumLib
{

/*! Interface for providing shape matrices and integration point values for
 *  extrapolation,
 *
 * Local assemblers that want to have some integration point values extrapolated
 * using Extrapolator have to implement this interface.
 *
 * \tparam PropertyTag  type of the property used to query a specific kind of
 *                      integration point value, usually an enum.
 */
template <typename PropertyTag>
class Extrapolatable
{
public:
    //! Provides the shape matrix at the given integration point.
    virtual Eigen::Map<const Eigen::RowVectorXd>
    getShapeMatrix(const unsigned integration_point) const = 0;

    /*! Provides integration point values for the given property.
     *
     * \remark
     * The vector \c cache can be used to store some derived quantities that
     * should be used as integration point values. They can be written to the
     * \c cache and a reference to the \c cache can then be returned by the
     * implementation of this method.
     * Ordinary integration point values can be returned directly (if they are
     * stored in the local assembler), and the \c cache cen be left untouched.
     *
     * \returns a reference to a vector containing the integration point values.
     */
    virtual std::vector<double> const&
    getIntegrationPointValues(PropertyTag const property,
                              std::vector<double>& cache) const = 0;

    virtual ~Extrapolatable() = default;
};

/*! Interface for classes that extrapolate integration point values to nodal
 *  values.
 *
 * \tparam PropertyTag    type of the property used to query a specific kind of
 *                        integration point value, usually an enum.
 * \tparam LocalAssembler type of the local assembler being queried for
 *                        integration point values.
 */
template<typename PropertyTag, typename LocalAssembler>
class Extrapolator
{
public:
    //! Vector of local assemblers, the same as is used in the FEM process.
    using LocalAssemblers = std::vector<std::unique_ptr<LocalAssembler>>;

    //! Extrapolates the given \c property from the given local assemblers.
    virtual void extrapolate(
            LocalAssemblers const& loc_asms, PropertyTag const property) = 0;

    /*! Computes residuals from the extrapolation of the given \c property.
     *
     * The residuals are computed as element values.
     *
     * \pre extrapolate() must have been called before with the same arguments.
     */
    virtual void calculateResiduals(
            LocalAssemblers const& loc_asms, PropertyTag const property) = 0;

    //! Returns the extrapolated nodal values.
    //! \todo Maybe write directly to a MeshProperty.
    virtual GlobalVector const& getNodalValues() const = 0;

    //! Returns the extrapolation residuals.
    //! \todo Maybe write directly to a MeshProperty.
    virtual GlobalVector const& getElementResiduals() const = 0;

    virtual ~Extrapolator() = default;
};

} // namespace ProcessLib

#endif // NUMLIB_EXTRAPOLATOR_H
