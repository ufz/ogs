/**
 * \copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Dense>
#include <vector>

#include "HTMaterialProperties.h"

#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Parameter/Parameter.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

#include "HTLocalAssemblerInterface.h"

namespace ProcessLib
{
namespace HT
{
template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class HTFEM : public HTLocalAssemblerInterface
{
    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;

    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using GlobalDimNodalMatrixType =
        typename ShapeMatricesType::GlobalDimNodalMatrixType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;

public:
    HTFEM(MeshLib::Element const& element,
          std::size_t const local_matrix_size,
          bool is_axially_symmetric,
          unsigned const integration_order,
          HTMaterialProperties const& material_properties,
          const unsigned dof_per_node)
        : HTLocalAssemblerInterface(),
          _element(element),
          _material_properties(material_properties),
          _integration_method(integration_order)
    {
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * dof_per_node);
        (void)local_matrix_size;

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();
        _ip_data.reserve(n_integration_points);

        auto const shape_matrices =
            initShapeMatrices<ShapeFunction, ShapeMatricesType,
                              IntegrationMethod, GlobalDim>(
                element, is_axially_symmetric, _integration_method);

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            _ip_data.emplace_back(
                shape_matrices[ip].N, shape_matrices[ip].dNdx,
                _integration_method.getWeightedPoint(ip).getWeight() *
                    shape_matrices[ip].integralMeasure *
                    shape_matrices[ip].detJ);
        }
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _ip_data[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

protected:
    MeshLib::Element const& _element;
    HTMaterialProperties const& _material_properties;

    IntegrationMethod const _integration_method;
    std::vector<
        IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType>,
        Eigen::aligned_allocator<
            IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType>>>
        _ip_data;

    double getHeatEnergyCoefficient(const double t, const SpatialPosition& pos,
                                    const double porosity,
                                    const double fluid_density,
                                    const double specific_heat_capacity_fluid)
    {
        auto const& material_properties = this->_material_properties;
        auto const specific_heat_capacity_solid =
            material_properties.specific_heat_capacity_solid(t, pos)[0];

        auto const solid_density = material_properties.density_solid(t, pos)[0];

        return solid_density * specific_heat_capacity_solid * (1 - porosity) +
               fluid_density * specific_heat_capacity_fluid * porosity;
    }

    GlobalDimMatrixType getThermalConductivityDispersivity(
        const double t, const SpatialPosition& pos, const double porosity,
        const double fluid_density, const double specific_heat_capacity_fluid,
        const GlobalDimVectorType& velocity, const GlobalDimMatrixType& I)
    {
        auto const& material_properties = this->_material_properties;

        auto const thermal_conductivity_solid =
            material_properties.thermal_conductivity_solid(t, pos)[0];
        auto const thermal_conductivity_fluid =
            material_properties.thermal_conductivity_fluid(t, pos)[0];
        double const thermal_conductivity =
            thermal_conductivity_solid * (1 - porosity) +
            thermal_conductivity_fluid * porosity;

        auto const thermal_dispersivity_longitudinal =
            material_properties.thermal_dispersivity_longitudinal(t, pos)[0];
        auto const thermal_dispersivity_transversal =
            material_properties.thermal_dispersivity_transversal(t, pos)[0];

        if (thermal_dispersivity_longitudinal == 0.0 &&
            thermal_dispersivity_transversal == 0.0)
            return thermal_conductivity * I;

        double const velocity_magnitude = velocity.norm();
        GlobalDimMatrixType const thermal_dispersivity =
            fluid_density * specific_heat_capacity_fluid *
            (thermal_dispersivity_transversal * velocity_magnitude * I +
             (thermal_dispersivity_longitudinal -
              thermal_dispersivity_transversal) /
                 velocity_magnitude * velocity * velocity.transpose());

        return thermal_conductivity * I + thermal_dispersivity;
    }

    std::vector<double> const& getIntPtDarcyVelocityLocal(
        const double t, std::vector<double> const& local_p,
        std::vector<double> const& local_T, std::vector<double>& cache) const
    {
        auto const n_integration_points =
            _integration_method.getNumberOfPoints();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<
            Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, GlobalDim, n_integration_points);

        SpatialPosition pos;
        pos.setElementID(_element.getID());

        MaterialLib::Fluid::FluidProperty::ArrayType vars;

        auto const p_nodal_values = Eigen::Map<const NodalVectorType>(
            &local_p[0], ShapeFunction::NPOINTS);

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            auto const& ip_data = _ip_data[ip];
            auto const& N = ip_data.N;
            auto const& dNdx = ip_data.dNdx;

            pos.setIntegrationPoint(ip);

            double T_int_pt = 0.0;
            double p_int_pt = 0.0;
            NumLib::shapeFunctionInterpolate(local_p, N, p_int_pt);
            NumLib::shapeFunctionInterpolate(local_T, N, T_int_pt);
            vars[static_cast<int>(
                MaterialLib::Fluid::PropertyVariableType::T)] = T_int_pt;
            vars[static_cast<int>(
                MaterialLib::Fluid::PropertyVariableType::p)] = p_int_pt;

            auto const K = _material_properties.porous_media_properties
                               .getIntrinsicPermeability(t, pos)
                               .getValue(t, pos, 0.0, T_int_pt);

            auto const mu = _material_properties.fluid_properties->getValue(
                MaterialLib::Fluid::FluidPropertyType::Viscosity, vars);
            GlobalDimMatrixType const K_over_mu = K / mu;

            cache_mat.col(ip).noalias() = -K_over_mu * dNdx * p_nodal_values;

            if (_material_properties.has_gravity)
            {
                auto const rho_w =
                    _material_properties.fluid_properties->getValue(
                        MaterialLib::Fluid::FluidPropertyType::Density, vars);
                auto const b = _material_properties.specific_body_force;
                // here it is assumed that the vector b is directed 'downwards'
                cache_mat.col(ip).noalias() += K_over_mu * rho_w * b;
            }
        }

        return cache;
    }
};

}  // namespace HT
}  // namespace ProcessLib
