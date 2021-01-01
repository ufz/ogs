/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <memory>
#include <vector>

#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "ProcessLib/LocalAssemblerInterface.h"

#include "IntegrationPointDataNonlocalInterface.h"

namespace ProcessLib
{
namespace SmallDeformationNonlocal
{
template <int DisplacementDim>
struct SmallDeformationNonlocalLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
    virtual std::size_t setIPDataInitialConditions(
        std::string const& name, double const* values,
        int const integration_order) = 0;

    virtual void setIPDataInitialConditionsFromCellData(
        std::string const& name, std::vector<double> const& value) = 0;

    virtual void computeCrackIntegral(
        std::size_t mesh_item_id,
        NumLib::LocalToGlobalIndexMap const& dof_table, GlobalVector const& x,
        double& crack_volume) = 0;

    virtual std::vector<double> const& getIntPtEpsPV(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;
    virtual std::vector<double> const& getIntPtEpsPDXX(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;
    virtual std::vector<double> getKappaD() const = 0;
    virtual std::vector<double> const& getIntPtDamage(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtFreeEnergyDensity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> getSigma() const = 0;
    virtual std::vector<double> const& getIntPtSigma(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    virtual std::vector<double> const& getIntPtEpsilon(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const = 0;

    // TODO move to NumLib::ExtrapolatableElement
    virtual unsigned getNumberOfIntegrationPoints() const = 0;

    virtual typename MaterialLib::Solids::MechanicsBase<
        DisplacementDim>::MaterialStateVariables const&
    getMaterialStateVariablesAt(int integration_point) const = 0;

    virtual std::vector<double> const& getNodalValues(
        std::vector<double>& nodal_values) const = 0;

    virtual void nonlocal(
        std::size_t const mesh_item_id,
        std::vector<std::unique_ptr<
            SmallDeformationNonlocalLocalAssemblerInterface>> const&
            local_assemblers) = 0;

    virtual void getIntegrationPointCoordinates(
        Eigen::Vector3d const& coords,
        std::vector<double>& distances) const = 0;

    virtual IntegrationPointDataNonlocalInterface* getIPDataPtr(
        int const ip) = 0;
};

}  // namespace SmallDeformationNonlocal
}  // namespace ProcessLib
