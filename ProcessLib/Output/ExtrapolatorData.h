/**
 * \file
 * \copyright
 * Copyright (c) 2012-2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

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
class ExtrapolatorData
{
public:
    ExtrapolatorData() = default;

    /*! Constructs a new instance.
     *
     * \param extrapolator the extrapolator managed by the instance being
     * created
     * \param dof_table_single_component the d.o.f. table used by the \c
     * extrapolator
     * \param manage_storage If true the memory of \c dof_table_single_component
     * will be freed by the destructor of this class.
     */
    ExtrapolatorData(
        std::unique_ptr<NumLib::Extrapolator>&& extrapolator,
        NumLib::LocalToGlobalIndexMap const* const dof_table_single_component,
        bool const manage_storage)
        : extrapolator_(std::move(extrapolator)),
          dof_table_single_component_(dof_table_single_component),
          manage_storage_(manage_storage)
    {
    }

    ExtrapolatorData(ExtrapolatorData&& other)
        : extrapolator_(std::move(other.extrapolator_)),
          dof_table_single_component_(other.dof_table_single_component_),
          manage_storage_(other.manage_storage_)
    {
        other.manage_storage_ = false;
        other.dof_table_single_component_ = nullptr;
    }

    ExtrapolatorData& operator=(ExtrapolatorData&& other)
    {
        cleanup();
        manage_storage_ = other.manage_storage_;
        dof_table_single_component_ = other.dof_table_single_component_;
        extrapolator_ = std::move(other.extrapolator_);
        other.dof_table_single_component_ = nullptr;
        other.manage_storage_ = false;
        return *this;
    }

    NumLib::LocalToGlobalIndexMap const& getDOFTable() const
    {
        return *dof_table_single_component_;
    }
    NumLib::Extrapolator& getExtrapolator() const { return *extrapolator_; }

    ~ExtrapolatorData() { cleanup(); }

private:
    //! Deletes the d.o.f table if it is allowed to do so.
    void cleanup()
    {
        if (manage_storage_)
        {
            delete dof_table_single_component_;
            dof_table_single_component_ = nullptr;
        }
    }

    //! Extrapolator managed by the ExtrapolatorData instance.
    std::unique_ptr<NumLib::Extrapolator> extrapolator_;

    //! D.o.f. table used by the extrapolator.
    NumLib::LocalToGlobalIndexMap const* dof_table_single_component_ = nullptr;

    //! If true, free storage of the d.o.f. table in the ExtrapolatorData
    //! destructor.
    bool manage_storage_ = false;
};

}  // namespace ProcessLib
