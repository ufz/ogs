/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <functional>
#include <vector>

#include "MathLib/LinAlg/GlobalMatrixVectorTypes.h"

#include "ExtrapolatableElement.h"

namespace NumLib
{
class LocalToGlobalIndexMap;

/*! Adaptor to get information needed by an Extrapolator from an "arbitrary"
 * collection of elements (e.g., local assemblers).
 *
 * This is an interface class; suitable implementations have to be provided for
 * concrete collections of extrapolatable elements.
 */
class ExtrapolatableElementCollection
{
public:
    //! Returns the shape matrix of the element with the given \c id at the
    //! given \c integration_point.
    virtual Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        std::size_t const id, unsigned const integration_point) const = 0;

    /*! Returns integration point values of some property of a specific element.
     *
     * \param id ID of the element of which the property is queried.
     * \param t The time used in the evaluation if time dependent quantities.
     * \param current_solution The current solution of a ProcessLib::Process;
     * more generally any nodal GlobalVector.
     * \param dof_table The processes d.o.f. table used to get each element's
     * local d.o.f. from \c current_solution.
     * \param cache Can be used to compute a property on-the-fly.
     *
     * \remark
     * \parblock
     * Since this method returns a reference, integration point values can not
     * be computed locally and just returned.
     *
     * Instead, the vector \c cache can be used to store some derived quantities
     * that should be used as integration point values. They can be written to
     * the \c cache and a reference to the \c cache can then be returned by the
     * implementation of this method. Ordinary integration point values can be
     * returned directly (if they are stored in the local assembler), and the
     * \c cache can be left untouched.
     *
     * For usage examples see the processes already implemented or the unit test
     * in TestExtrapolation.cpp
     * \endparblock
     */
    virtual std::vector<double> const& getIntegrationPointValues(
        std::size_t const id, const double t,
        GlobalVector const& current_solution,
        LocalToGlobalIndexMap const& dof_table,
        std::vector<double>& cache) const = 0;

    //! Returns the number of elements whose properties shall be extrapolated.
    virtual std::size_t size() const = 0;

    virtual ~ExtrapolatableElementCollection() = default;
};

template <typename LocalAssemblerCollection>
class ExtrapolatableLocalAssemblerCollection
    : public ExtrapolatableElementCollection
{
public:
    //! LocalAssemblerCollection contains many LocalAssembler's.
    using LocalAssembler =
        typename std::decay<decltype(*std::declval<LocalAssemblerCollection>()
                                         [static_cast<std::size_t>(0)])>::type;

    static_assert(std::is_base_of<ExtrapolatableElement, LocalAssembler>::value,
                  "Local assemblers used for extrapolation must be derived "
                  "from ExtrapolatableElement.");

    /*! A method providing integration point values of some property.
     *
     * The method must return reference to a vector containing the integration
     * point values.
     *
     * For further information about the \c cache parameter see
     * ExtrapolatableElementCollection::getIntegrationPointValues().
     */
    using IntegrationPointValuesMethod =
        std::function<std::vector<double> const&(
            LocalAssembler const& loc_asm, const double t,
            GlobalVector const& current_solution,
            NumLib::LocalToGlobalIndexMap const& dof_table,
            std::vector<double>& cache)>;

    /*! Constructs a new instance.
     *
     * \param local_assemblers a collection of local assemblers whose
     * integration point values shall be extrapolated.
     * \param integration_point_values_method a LocalAssembler method returning
     * integration point values of some property for that local assembler.
     */
    ExtrapolatableLocalAssemblerCollection(
        LocalAssemblerCollection const& local_assemblers,
        IntegrationPointValuesMethod const& integration_point_values_method)
        : _local_assemblers(local_assemblers),
          _integration_point_values_method{integration_point_values_method}
    {
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        std::size_t const id, unsigned const integration_point) const override
    {
        ExtrapolatableElement const& loc_asm = *_local_assemblers[id];
        return loc_asm.getShapeMatrix(integration_point);
    }

    std::vector<double> const& getIntegrationPointValues(
        std::size_t const id, const double t,
        GlobalVector const& current_solution,
        LocalToGlobalIndexMap const& dof_table,
        std::vector<double>& cache) const override
    {
        auto const& loc_asm = *_local_assemblers[id];
        return _integration_point_values_method(loc_asm, t, current_solution,
                                                dof_table, cache);
    }

    std::size_t size() const override { return _local_assemblers.size(); }
private:
    LocalAssemblerCollection const& _local_assemblers;
    IntegrationPointValuesMethod const _integration_point_values_method;
};

//! Creates an ExtrapolatableLocalAssemblerCollection, which can be used to
//! provide information to an Extrapolator.
template <typename LocalAssemblerCollection,
          typename IntegrationPointValuesMethod>
ExtrapolatableLocalAssemblerCollection<LocalAssemblerCollection>
makeExtrapolatable(LocalAssemblerCollection const& local_assemblers,
                   IntegrationPointValuesMethod integration_point_values_method)
{
    return ExtrapolatableLocalAssemblerCollection<LocalAssemblerCollection>{
        local_assemblers, integration_point_values_method};
}
}  // namespace NumLib
