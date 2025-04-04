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

#include <cstddef>
#include <vector>

namespace BaseLib
{
/**
 * Interface class for subdivision operators
 */
class ISubdivision
{
public:
    /// Returns a vector of subdivided points
    virtual std::vector<double> operator()() const = 0;

    virtual ~ISubdivision() = default;
};

/**
 * Uniform subdivision operator
 */
class UniformSubdivision : public ISubdivision
{
public:
    /**
     * Configuration
     * @param length          total length to be subdivided
     * @param n_subdivision   the number of subdivision
     */
    UniformSubdivision(double length, std::size_t n_subdivision)
        : length_(length), n_subdivision_(n_subdivision)
    {
    }

    /// Returns a vector of subdivided points
    std::vector<double> operator()() const override
    {
        std::vector<double> x;
        x.reserve(n_subdivision_ + 1);
        const double dL = length_ / static_cast<double>(n_subdivision_);
        for (std::size_t i = 0; i < n_subdivision_ + 1; i++)
        {
            x.push_back(i * dL);
        }
        return x;
    }

private:
    const double length_;
    const std::size_t n_subdivision_;
};

/**
 * Gradual subdivision operator with a constant multiplier_
 */
class GradualSubdivision : public ISubdivision
{
public:
    /**
     * Constructor
     * @param L           total length to be subdivided
     * @param dL0         initial cell length
     * @param max_dL      maximum cell length
     * @param multiplier  multiplier to cell length
     */
    GradualSubdivision(const double L,
                       const double dL0,
                       const double max_dL,
                       const double multiplier);

    /// Returns a vector of subdivided points
    std::vector<double> operator()() const override;

private:
    const double length_;
    const double dL0_;
    const double max_dL_;
    const double multiplier_;
};

/**
 * Gradual subdivision operator with a constant multiplier_.
 *
 * In this class the number of subdivisions is known a priori.
 */
class GradualSubdivisionFixedNum : public ISubdivision
{
public:
    /**
     * Constructor
     * @param L                 total length to be subdivided
     * @param num_subdivisions  number of subdivisions to generate
     * @param multiplier        multiplier to cell length
     */
    GradualSubdivisionFixedNum(const double L,
                               const std::size_t num_subdivisions,
                               const double multiplier);

    /// Returns a vector of subdivided points
    std::vector<double> operator()() const override;

private:
    const double length_;
    const std::size_t num_subdivisions_;
    const double multiplier_;
};

}  // namespace BaseLib
