/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef NUMLIB_LOCAL_LLSQ_EXTRAPOLATOR_H
#define NUMLIB_LOCAL_LLSQ_EXTRAPOLATOR_H

#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/DOF/GlobalMatrixProviders.h"
#include "Extrapolator.h"

namespace NumLib
{
/*! Extrapolator doing a linear least squares extrapolation locally in each
 *  element.
 *
 * The results of the element-wise least squares are averaged at each node.
 *
 * \note
 * \parblock
 * The described procedure is not least-squares optimal on the whole mesh! But
 * the residuals are computed from the actual interpolation result that is being
 * returned by this class.
 *
 * Currently this class only supports interpolating single component variables.
 *
 * Furthermore, the number of integration points in each element must be greater
 * than or equal to the number of nodes of that element. This restriction is due to
 * the use of the least squares which requires an exact or overdetermined equation
 * system.
 * \endparblock
 */
class LocalLinearLeastSquaresExtrapolator : public Extrapolator
{
public:
    /*! Constructs a new instance
     *
     * \note
     * The \c dof_table of \c matrix_specs must be set, and it must point to a
     * dof_table for one single component variable.
     */
    explicit LocalLinearLeastSquaresExtrapolator(
        MathLib::MatrixSpecifications const& matrix_specs,
        NumLib::LocalToGlobalIndexMap const& dof_table)
        : _nodal_values(
              NumLib::GlobalVectorProvider::provider.getVector(
                  matrix_specs))
#ifndef USE_PETSC
          ,
          _residuals(dof_table.size())
#else
          ,
          _residuals(dof_table.size(), false)
#endif
          ,
          _local_to_global(dof_table)
    {
        /* Note in case the following assertion fails.
         * If you copied the extrapolation code, for your processes from
         * somewhere, note that the code from the groundwater flow process might
         * not suit your needs: It is a special case and is therefore most
         * likely too simplistic. You better adapt the extrapolation code from
         * some more advanced process, like the TES process.
         */
        assert(dof_table.getNumberOfComponents() == 1 &&
               "The d.o.f. table passed must be for one variable that has "
               "only one component!");
    }

    void extrapolate(
            ExtrapolatableElementCollection const& extrapolatables) override;

    /*! \copydoc Extrapolator::calculateResiduals()
     *
     * The computed residuals are root-mean-square of the difference between
     * the integration point values obtained from the local assemblers and the
     * extrapolation results when interpolated back to the integration points
     * again.
     */
    void calculateResiduals(
            ExtrapolatableElementCollection const& extrapolatables) override;

    GlobalVector const& getNodalValues() const override
    {
        return _nodal_values;
    }

    GlobalVector const& getElementResiduals() const override
    {
        return _residuals;
    }

    ~LocalLinearLeastSquaresExtrapolator()
    {
        NumLib::GlobalVectorProvider::provider.releaseVector(
            _nodal_values);
    }

private:
    //! Extrapolate one element.
    void extrapolateElement(
        std::size_t const element_index,
        ExtrapolatableElementCollection const& extrapolatables,
        GlobalVector& counts);

    //! Compute the residuals for one element
    void calculateResiudalElement(
        std::size_t const element_index,
        ExtrapolatableElementCollection const& extrapolatables);

    GlobalVector& _nodal_values;  //!< extrapolated nodal values
    GlobalVector _residuals;      //!< extrapolation residuals

    //! DOF table used for writing to global vectors.
    NumLib::LocalToGlobalIndexMap const& _local_to_global;

    //! Avoids frequent reallocations.
    Eigen::MatrixXd _local_matrix_cache;

    //! Avoids frequent reallocations.
    std::vector<double> _integration_point_values_cache;
};

}  // namespace NumLib

#endif  // NUMLIB_LOCAL_LLSQ_EXTRAPOLATOR_H
