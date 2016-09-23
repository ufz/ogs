/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_PROPERTYINDEXPARAMETER_H_
#define PROCESSLIB_PROPERTYINDEXPARAMETER_H_

#include "Parameter.h"

#include "MeshLib/PropertyVector.h"

namespace MeshLib
{
template <typename T>
class PropertyVector;
}  // MeshLib

namespace ProcessLib
{

/// A parameter class looking for values from index in a property vector.
/// This class can be used for material ID depependent parameters.
template <typename T, MeshLib::MeshItemType MeshItemType>
struct PropertyIndexParameter final
    : public Parameter<T>
{
    /**
     * Constructing from a property vector of index and corresponding values
     *
     * @param property    a property vector of index for mesh items
     * @param vec_values  a vector of values for each index
     */
    PropertyIndexParameter(MeshLib::PropertyVector<int> const& property,
    std::vector<std::vector<double>> const& vec_values)
        : _property_index(property), _vec_values(vec_values)
    {
    }

    bool isTimeDependent() const override { return false; }

    unsigned getNumberOfComponents() const override
    {
        return _vec_values.empty() ? 0 : _vec_values.front().size();
    }

    std::vector<T> const& operator()(double const /*t*/,
                                     SpatialPosition const& pos) const override
    {
        auto const item_id = getMeshItemID(pos, type<MeshItemType>());
        assert(item_id);
        int index = _property_index[item_id.get()];
        return _vec_values[index];
    }

private:
    template <MeshLib::MeshItemType ITEM_TYPE> struct type {};

    boost::optional<std::size_t>
    getMeshItemID(SpatialPosition const& pos, type<MeshLib::MeshItemType::Cell>) const
    {
        return pos.getElementID();
    }

    boost::optional<std::size_t>
    getMeshItemID(SpatialPosition const& pos, type<MeshLib::MeshItemType::Node>) const
    {
        return pos.getNodeID();
    }

    MeshLib::PropertyVector<int> const& _property_index;
    std::vector<std::vector<T>> const _vec_values;
};


std::unique_ptr<ParameterBase> createPropertyIndexParameter(
    BaseLib::ConfigTree const& config, MeshLib::Mesh const& mesh);

}  // ProcessLib

#endif // PROCESSLIB_PROPERTYINDEXPARAMETER_H_
