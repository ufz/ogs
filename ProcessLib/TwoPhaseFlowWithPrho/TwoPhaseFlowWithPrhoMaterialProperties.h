/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <memory>
#include <vector>
#include "MaterialLib/Fluid/FluidPropertyHeaders.h"
#include "MaterialLib/PhysicalConstant.h"
#include "MaterialLib/PorousMedium/Porosity/Porosity.h"
#include "MaterialLib/PorousMedium/PorousPropertyHeaders.h"
#include "MaterialLib/PorousMedium/Storage/Storage.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CapillaryPressureSaturation.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/CapillaryPressure/CreateCapillaryPressureModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/CreateRelativePermeabilityModel.h"
#include "MaterialLib/PorousMedium/UnsaturatedProperty/RelativePermeability/RelativePermeability.h"
#include "MathLib/InterpolationAlgorithms/PiecewiseLinearInterpolation.h"

namespace MeshLib
{
template <typename PROP_VAL_TYPE>
class PropertyVector;
}

namespace ProcessLib
{
class SpatialPosition;
namespace TwoPhaseFlowWithPrho
{
class TwoPhaseFlowWithPrhoMaterialProperties
{
public:
    using ArrayType = MaterialLib::Fluid::FluidProperty::ArrayType;

    TwoPhaseFlowWithPrhoMaterialProperties(
        boost::optional<MeshLib::PropertyVector<int> const&> const material_ids,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>
            liquid_density,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>
            viscosity,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>
            gas_density,
        std::unique_ptr<MaterialLib::Fluid::FluidProperty>
            gas_viscosity,
        std::vector<std::unique_ptr<MaterialLib::PorousMedium::Permeability>>&&
            intrinsic_permeability_models,
        std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>&&
            porosity_models,
        std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>&&
            storage_models,
        std::vector<std::unique_ptr<
            MaterialLib::PorousMedium::CapillaryPressureSaturation>>&&
            capillary_pressure_models,
        std::vector<
            std::unique_ptr<MaterialLib::PorousMedium::RelativePermeability>>&&
            relative_permeability_models);

    int getMaterialID(const std::size_t element_id);

    Eigen::MatrixXd const& getPermeability(
        const int material_id,
        const double t,
        const ProcessLib::SpatialPosition& pos,
        const int dim) const;

    double getPorosity(const int material_id, const double t,
                       const ProcessLib::SpatialPosition& pos, const double p,
                       const double T, const double porosity_variable) const;

    double getNonwetRelativePermeability(const double t,
                                         const ProcessLib::SpatialPosition& pos,
                                         const double p, const double T,
                                         const double saturation) const;
    double getWetRelativePermeability(const double t,
                                      const ProcessLib::SpatialPosition& pos,
                                      const double p, const double T,
                                      const double saturation) const;
    double getCapillaryPressure(const int material_id, const double t,
                                const ProcessLib::SpatialPosition& pos,
                                const double p, const double T,
                                const double saturation) const;
    double getCapillaryPressureDerivative(
        const int material_id, const double t,
        const ProcessLib::SpatialPosition& pos, const double p, const double T,
        const double saturation) const;
    double getLiquidDensity(const double p, const double T) const;
    double getGasDensity(const double p, const double T) const;
    double getGasViscosity(const double p, const double T) const;
    double getLiquidViscosity(const double p, const double T) const;
    double getDerivGasDensity(double const p, double const T) const;
    bool computeConstitutiveRelation(
        double const t,
        ProcessLib::SpatialPosition const& x_position,
        const int material_id,
        double const pg,
        double const X,
        double const T,
        double& Sw,
        double& X_m,
        double& dsw_dpg,
        double& dsw_dX,
        double& dxm_dpg,
        double& dxm_dX);

protected:
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _liquid_density;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _viscosity;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _gas_density;
    std::unique_ptr<MaterialLib::Fluid::FluidProperty> _gas_viscosity;

    /** Use two phase models for different material zones.
    *  Material IDs must be given as mesh element properties.
    */
    boost::optional<MeshLib::PropertyVector<int> const&> const _material_ids;

    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Permeability>>
        _intrinsic_permeability_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Porosity>>
        _porosity_models;
    std::vector<std::unique_ptr<MaterialLib::PorousMedium::Storage>>
        _storage_models;
    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::CapillaryPressureSaturation>>
        _capillary_pressure_models;
    std::vector<
        std::unique_ptr<MaterialLib::PorousMedium::RelativePermeability>>
        _relative_permeability_models;

private:
    static int const jacobian_residual_size = 2;
    using ResidualVector = Eigen::Matrix<double, jacobian_residual_size, 1>;
    using JacobianMatrix =
        Eigen::Matrix<double, jacobian_residual_size, jacobian_residual_size,
                      Eigen::RowMajor>;
    using UnknownVector = Eigen::Matrix<double, jacobian_residual_size, 1>;

private:
    /**
    * Calculates the residual vector.
    */
    void calculateResidual(const int material_id, double const pl,
                           double const X, double const T, double Sw,
                           double X_m, ResidualVector& res);
    /**
    * Calculates the Jacobian.
    */
    void calculateJacobian(const int material_id, double const t,
                           ProcessLib::SpatialPosition const& x,
                           double const pl, double const X, double const T,
                           JacobianMatrix& Jac, double Sw, double X_m);
    /** Complementary condition 1
    * for calculating molar fraction of light component in the liquid phase
    */
    double calculateEquilibiumRhoWetLight(double const pg, double const Sw,
                                          double const rho_wet_h2) const;
    /** Complementary condition 2
    * for calculating the saturation
    */
    double calculateSaturation(double /*PL*/, double X, double Sw,
                               double rho_wet_h2, double rho_nonwet_h2,
                               double /*T*/) const;
    /**
    * Calculate the derivatives using the analytical way
    */
    double calculatedSwdP(double pl, double S, double rho_wet_h2,
                          double const T, int current_material_id) const;
    /**
    * Calculate the derivatives using the analytical way
    */
    double calculatedSwdX(double const pl, const double /*X*/, const double S,
                          const double rho_wet_h2, double const T,
                          int current_material_id) const;
    /**
    * Calculate the derivatives using the analytical way
    */
    double calculatedXmdX(double pl, double Sw, double rho_wet_h2, double dSwdX,
                          int current_material_id) const;
    /**
    * Calculate the derivatives using the analytical way
    */
    double calculatedXmdP(double pl, double Sw, double rho_wet_h2, double dSwdP,
                          int current_material_id) const;
};

}  // end of namespace
}  // end of namespace
