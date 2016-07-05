/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESSLIB_EXTRAPOLATORDATA_H
#define PROCESSLIB_EXTRAPOLATORDATA_H

#include <memory>
#include "NumLib/DOF/LocalToGlobalIndexMap.h"
#include "NumLib/Extrapolation/Extrapolator.h"

namespace ProcessLib
{
/*! Helper struct containing an extrapolator and a single component DOF table.
 *
 * Storage for the DOF table is managed optionally.
 *
 * \todo Later on this struct shall be moved, e.g., be merged with the process
 * output class.
 */
struct ExtrapolatorData
{
    ExtrapolatorData() = default;

    ExtrapolatorData(
        std::unique_ptr<NumLib::Extrapolator>&& extrapolator_,
        NumLib::LocalToGlobalIndexMap const* const dof_table_single_component,
        bool const manage_storage)
        : extrapolator(std::move(extrapolator_))
        , _dof_table_single_component(dof_table_single_component)
        , _manage_storage(manage_storage)
    {
    }

    ExtrapolatorData(ExtrapolatorData&& other)
        : extrapolator(std::move(other.extrapolator))
        , _dof_table_single_component(other._dof_table_single_component)
        , _manage_storage(other._manage_storage)
    {
        other._manage_storage = false;
        other._dof_table_single_component = nullptr;
    }

    ExtrapolatorData& operator=(ExtrapolatorData&& other)
    {
        cleanup();
        _manage_storage = other._manage_storage;
        _dof_table_single_component = other._dof_table_single_component;
        extrapolator = std::move(other.extrapolator);
        return *this;
    }

    ~ExtrapolatorData() { cleanup(); }
    std::unique_ptr<NumLib::Extrapolator> extrapolator;

private:
    void cleanup()
    {
        if (_manage_storage) {
            delete _dof_table_single_component;
        }
    }

    NumLib::LocalToGlobalIndexMap const* _dof_table_single_component = nullptr;
    bool _manage_storage = false;
};

} // namespace ProcessLib


#endif // PROCESSLIB_EXTRAPOLATORDATA_H
