/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>

#include "BaseLib/Logging.h"

namespace ProcessLib::Assembly
{
struct Stats
{
    std::size_t count = 0;
    std::size_t count_nonzero = 0;
    std::size_t count_global = 0;

    Stats& operator+=(Stats const& other)
    {
        count += other.count;
        count_nonzero += other.count_nonzero;
        count_global += other.count_global;

        return *this;
    }

    void print(std::string const& matrix_or_vector_name) const
    {
        DBUG("Stats [{}]: {} elements added to the matrix cache.",
             matrix_or_vector_name,
             count);
        DBUG("Stats [{}]: {} nonzero elements added to the matrix cache.",
             matrix_or_vector_name,
             count_nonzero);
        DBUG("Stats [{}]: {} elements added to the global matrix.",
             matrix_or_vector_name,
             count_global);
    }
};

struct MultiStats
{
    Stats M;
    Stats K;
    Stats b;
    Stats Jac;

    MultiStats& operator+=(MultiStats const& other)
    {
        M += other.M;
        K += other.K;
        b += other.b;
        Jac += other.Jac;

        return *this;
    }

    void print() const
    {
        M.print("M");
        K.print("K");
        b.print("b");
        Jac.print("J");
    }
};

template <typename Data>
class CumulativeStats
    : public std::enable_shared_from_this<CumulativeStats<Data>>
{
    using Base = std::enable_shared_from_this<CumulativeStats<Data>>;

public:
    Data data;

    static std::shared_ptr<CumulativeStats<Data>> create()
    {
        return std::shared_ptr<CumulativeStats<Data>>(
            new CumulativeStats<Data>());
    }

    // Could return unique_ptr, but shared_ptr is more consistent with the
    // create() method.
    std::shared_ptr<CumulativeStats<Data>> clone()
    {
        return std::make_shared<CumulativeStats<Data>>(*this);
    }

    CumulativeStats(CumulativeStats<Data> const& other) = delete;

    CumulativeStats(CumulativeStats<Data>& other)
        : Base{other},
          data{},
          parent_{other.parent_ ? other.parent_ : other.shared_from_this()},
          parent_mutex_{other.parent_mutex_}
    {
    }

    CumulativeStats(CumulativeStats<Data>&& other)
        : parent_{std::move(other.parent_)},
          parent_mutex_{std::move(other.parent_mutex_)}
    {
        std::swap(data, other.data);
    }

    ~CumulativeStats()
    {
        if (!parent_)
        {
            return;
        }

        std::lock_guard<std::mutex> const lock(*parent_mutex_);

        DBUG("Adding cumulative stats to parent.");

        parent_->data += data;
    }

    void print() const { data.print(); }

private:
    CumulativeStats() : parent_mutex_{std::make_shared<std::mutex>()} {}

    std::shared_ptr<CumulativeStats<Data>> parent_;
    std::shared_ptr<std::mutex> parent_mutex_;
};
}  // namespace ProcessLib::Assembly
