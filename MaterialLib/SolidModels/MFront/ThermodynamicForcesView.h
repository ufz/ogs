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

#include <Eigen/Core>
#include <boost/mp11.hpp>

namespace MaterialLib::Solids::MFront
{
/// Used for disambiguation with MFront's tangent operator blocks data.
struct OGSMFrontThermodynamicForcesData
{
    std::vector<double> data;
};

/**
 * Provides convenient access to the individual blocks of MFront's thermodynamic
 * forces data.
 *
 * \tparam DisplacementDim the displacement dimension
 * \tparam TDynForces a list (of types) of thermodynamic forces
 *
 * "Thermodynamic forces" is MFront nomenclature for the data an MFront
 * behaviour computes.
 *
 * \c TDynForces is a list of types. Each is expected to behave like Strain and
 * the other classes in Variable.h.
 *
 * MFront's thermodynamic forces are stored in a single vector of double values.
 * The data of forces that come earlier in \c TDynForces are stored earlier in
 * MFront's data.
 */
template <int DisplacementDim, typename TDynForces>
class OGSMFrontThermodynamicForcesView
{
    static_assert(boost::mp11::mp_is_set<TDynForces>::value,
                  "The types in the list TDynForces are not unique.");

    /// A template metafunction accessing the size of the given \c Variable.
    template <typename Variable>
    struct SizeOf
        : boost::mp11::mp_size_t<Variable::template size<DisplacementDim>()>
    {
    };

    /// Computes the offset of the given \c Force's data in MFront's
    /// thermodynamic forces data.
    template <typename Force>
    static constexpr std::size_t dataOffset()
    {
        using namespace boost::mp11;
        static_assert(mp_contains<TDynForces, Force>::value,
                      "The type Force is not in the list TDynForces.");

        using ForceIndex = mp_find<TDynForces, Force>;

        using ForcesHead = mp_take<TDynForces, ForceIndex>;

        using Sizes = mp_transform<SizeOf, ForcesHead>;

        return mp_apply<mp_plus, Sizes>::value;
    }

    /// Access a block of the given \c data as an Eigen::Map.
    template <typename Force, typename DataVector>
    static constexpr auto asEigenMap(DataVector& data)
    {
        static_assert(Force::template size<DisplacementDim>() != 1,
                      "Use asDouble for the single component case.");

        assert(data.size() == data_size_all_forces);

        constexpr std::size_t data_offset = dataOffset<Force>();

        constexpr std::size_t rows = Force::template rows<DisplacementDim>();
        constexpr std::size_t cols = Force::template cols<DisplacementDim>();
        constexpr auto order = cols == 1 ? Eigen::ColMajor : Eigen::RowMajor;

        using MatrixType =
            std::conditional_t<std::is_const_v<DataVector>,
                               const Eigen::Matrix<double, rows, cols, order>,
                               Eigen::Matrix<double, rows, cols, order>>;

        return Eigen::Map<MatrixType>(data.data() + data_offset);
    }

    /// Access a block of the given \c data as a double value.
    template <typename Force, typename DataVector>
    static constexpr auto& asDouble(DataVector& data)
    {
        static_assert(Force::template size<DisplacementDim>() == 1,
                      "Use asEigenMap for the multi component case.");

        assert(data.size() == data_size_all_forces);

        constexpr std::size_t data_offset = dataOffset<Force>();

        return data[data_offset];
    }

public:
    /// Read-only access to the data for the given thermodynamic force \c Force.
    template <typename Force>
    auto block(Force force, OGSMFrontThermodynamicForcesData const& data) const
    {
        return block(force, data.data);
    }

    /// Overload taking a std::vector.
    template <typename Force>
    auto block(Force, std::vector<double> const& data) const
    {
        constexpr std::size_t data_size =
            Force::template size<DisplacementDim>();

        if constexpr (data_size == 1)
        {
            return asDouble<Force>(data);
        }
        else
        {
            return asEigenMap<Force>(data);
        }
    }

    /// Read-write access to the data for the given thermodynamic force
    /// \c Force.
    ///
    /// This overload is chosen if \c Force has one component.
    template <typename Force,
              std::enable_if_t<SizeOf<Force>::value == 1, bool> = true>
    double& block(Force force, OGSMFrontThermodynamicForcesData& data) const
    {
        return block(force, data.data);
    }

    /// Overload taking a std::vector.
    template <typename Force,
              std::enable_if_t<SizeOf<Force>::value == 1, bool> = true>
    double& block(Force, std::vector<double>& data) const
    {
        return asDouble<Force>(data);
    }

    /// Read-write access to the data for the given thermodynamic force \c
    /// Force.
    ///
    /// This overload is chosen if \c Force has more than one component.
    template <typename Force,
              std::enable_if_t<SizeOf<Force>::value != 1, bool> = true>
    auto block(Force force, OGSMFrontThermodynamicForcesData& data) const
    {
        return block(force, data.data);
    }

    /// Overload taking a std::vector.
    template <typename Force,
              std::enable_if_t<SizeOf<Force>::value != 1, bool> = true>
    auto block(Force, std::vector<double>& data) const
    {
        return asEigenMap<Force>(data);
    }

    /// The passed data to the block() methods must have this size.
    static constexpr std::size_t data_size_all_forces = boost::mp11::mp_apply<
        boost::mp11::mp_plus,
        boost::mp11::mp_transform<SizeOf, TDynForces>>::value;
};
}  // namespace MaterialLib::Solids::MFront
