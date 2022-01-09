/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RichardsComponentTransportFEM.h"

#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Property.h"
#include "MaterialLib/MPL/Utils/FormEigenTensor.h"

namespace ProcessLib
{
namespace RichardsComponentTransport
{
template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
LocalAssemblerData<ShapeFunction, IntegrationMethod, GlobalDim>::
    LocalAssemblerData(
        MeshLib::Element const& element,
        std::size_t const local_matrix_size,
        bool is_axially_symmetric,
        unsigned const integration_order,
        RichardsComponentTransportProcessData const& process_data,
        ProcessVariable const& transport_process_variable)
    : _element_id(element.getID()),
      _process_data(process_data),
      _integration_method(integration_order),
      _transport_process_variable(transport_process_variable)
{
    // This assertion is valid only if all nodal d.o.f. use the same shape
    // matrices.
    assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);
    (void)local_matrix_size;

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();
    _ip_data.reserve(n_integration_points);

    auto const shape_matrices =
        NumLib::initShapeMatrices<ShapeFunction, ShapeMatricesType, GlobalDim>(
            element, is_axially_symmetric, _integration_method);

    for (unsigned ip = 0; ip < n_integration_points; ip++)
    {
        auto const& sm = shape_matrices[ip];
        double const integration_factor = sm.integralMeasure * sm.detJ;
        _ip_data.emplace_back(
            sm.N, sm.dNdx,
            _integration_method.getWeightedPoint(ip).getWeight() *
                integration_factor,
            sm.N.transpose() * sm.N * integration_factor *
                _integration_method.getWeightedPoint(ip).getWeight());
    }
}

template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
void LocalAssemblerData<ShapeFunction, IntegrationMethod, GlobalDim>::assemble(
    double const t, double const dt, std::vector<double> const& local_x,
    std::vector<double> const& /*local_xdot*/,
    std::vector<double>& local_M_data, std::vector<double>& local_K_data,
    std::vector<double>& local_b_data)
{
    auto const local_matrix_size = local_x.size();
    // This assertion is valid only if all nodal d.o.f. use the same shape
    // matrices.
    assert(local_matrix_size == ShapeFunction::NPOINTS * NUM_NODAL_DOF);

    auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_M_data, local_matrix_size, local_matrix_size);
    auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
        local_K_data, local_matrix_size, local_matrix_size);
    auto local_b = MathLib::createZeroedVector<LocalVectorType>(
        local_b_data, local_matrix_size);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element_id);

    auto p_nodal_values = Eigen::Map<const NodalVectorType>(
        &local_x[pressure_index], pressure_size);

    auto const& b = _process_data.specific_body_force;

    MaterialPropertyLib::VariableArray vars;

    GlobalDimMatrixType const& I(
        GlobalDimMatrixType::Identity(GlobalDim, GlobalDim));

    // Get material properties
    auto const& medium = *_process_data.media_map->getMedium(_element_id);
    auto const& phase = medium.phase("AqueousLiquid");
    auto const& component =
        phase.component(_transport_process_variable.getName());

    auto KCC = local_K.template block<concentration_size, concentration_size>(
        concentration_index, concentration_index);
    auto MCC = local_M.template block<concentration_size, concentration_size>(
        concentration_index, concentration_index);
    auto Kpp = local_K.template block<pressure_size, pressure_size>(
        pressure_index, pressure_index);
    auto Mpp = local_M.template block<pressure_size, pressure_size>(
        pressure_index, pressure_index);
    auto Bp = local_b.template block<pressure_size, 1>(pressure_index, 0);

    for (std::size_t ip(0); ip < n_integration_points; ++ip)
    {
        pos.setIntegrationPoint(ip);

        auto const& ip_data = _ip_data[ip];
        auto const& N = ip_data.N;
        auto const& dNdx = ip_data.dNdx;
        auto const& w = ip_data.integration_weight;

        double C_int_pt = 0.0;
        double p_int_pt = 0.0;
        // Order matters: First C, then p!
        NumLib::shapeFunctionInterpolate(local_x, N, C_int_pt, p_int_pt);

        vars[static_cast<int>(
            MaterialPropertyLib::Variable::capillary_pressure)] = -p_int_pt;
        auto const Sw = medium[MaterialPropertyLib::PropertyType::saturation]
                            .template value<double>(vars, pos, t, dt);

        double const dSw_dpc =
            medium[MaterialPropertyLib::PropertyType::saturation]
                .template dValue<double>(
                    vars, MaterialPropertyLib::Variable::capillary_pressure,
                    pos, t, dt);

        vars[static_cast<int>(MaterialPropertyLib::Variable::concentration)] =
            C_int_pt;
        vars[static_cast<int>(MaterialPropertyLib::Variable::phase_pressure)] =
            p_int_pt;

        // \todo the argument to getValue() has to be changed for non
        // constant storage model
        auto const specific_storage =
            medium[MaterialPropertyLib::PropertyType::storage]
                .template value<double>(vars, pos, t, dt);
        // \todo the first argument has to be changed for non constant
        // porosity model
        auto const porosity =
            medium[MaterialPropertyLib::PropertyType::porosity]
                .template value<double>(vars, pos, t, dt);

        auto const retardation_factor =
            component[MaterialPropertyLib::PropertyType::retardation_factor]
                .template value<double>(vars, pos, t, dt);

        auto const solute_dispersivity_transverse =
            medium[MaterialPropertyLib::PropertyType::transversal_dispersivity]
                .template value<double>(vars, pos, t, dt);
        auto const solute_dispersivity_longitudinal =
            medium[MaterialPropertyLib::PropertyType::longitudinal_dispersivity]
                .template value<double>(vars, pos, t, dt);

        // Use the fluid density model to compute the density
        auto const density = phase[MaterialPropertyLib::PropertyType::density]
                                 .template value<double>(vars, pos, t, dt);
        auto const decay_rate =
            component[MaterialPropertyLib::PropertyType::decay_rate]
                .template value<double>(vars, pos, t, dt);
        auto const pore_diffusion_coefficient =
            MaterialPropertyLib::formEigenTensor<GlobalDim>(
                component[MaterialPropertyLib::PropertyType::pore_diffusion]
                    .value(vars, pos, t, dt));

        auto const K = MaterialPropertyLib::formEigenTensor<GlobalDim>(
            medium[MaterialPropertyLib::PropertyType::permeability].value(
                vars, pos, t, dt));
        vars[static_cast<int>(
            MaterialPropertyLib::Variable::liquid_saturation)] = Sw;
        auto const k_rel =
            medium[MaterialPropertyLib::PropertyType::relative_permeability]
                .template value<double>(vars, pos, t, dt);
        // Use the viscosity model to compute the viscosity
        auto const mu = phase[MaterialPropertyLib::PropertyType::viscosity]
                            .template value<double>(vars, pos, t, dt);
        auto const K_times_k_rel_over_mu = K * (k_rel / mu);

        GlobalDimVectorType const velocity =
            _process_data.has_gravity
                ? GlobalDimVectorType(-K_times_k_rel_over_mu *
                                      (dNdx * p_nodal_values - density * b))
                : GlobalDimVectorType(-K_times_k_rel_over_mu * dNdx *
                                      p_nodal_values);

        double const velocity_magnitude = velocity.norm();
        GlobalDimMatrixType const hydrodynamic_dispersion =
            velocity_magnitude != 0.0
                ? GlobalDimMatrixType(
                      porosity * pore_diffusion_coefficient +
                      solute_dispersivity_transverse * velocity_magnitude * I +
                      (solute_dispersivity_longitudinal -
                       solute_dispersivity_transverse) /
                          velocity_magnitude * velocity * velocity.transpose())
                : GlobalDimMatrixType(porosity * pore_diffusion_coefficient +
                                      solute_dispersivity_transverse *
                                          velocity_magnitude * I);

        // matrix assembly
        KCC.noalias() +=
            (dNdx.transpose() * hydrodynamic_dispersion * dNdx +
             N.transpose() * velocity.transpose() * dNdx +
             N.transpose() * decay_rate * porosity * retardation_factor * N) *
            w;
        MCC.noalias() += w * N.transpose() * porosity * retardation_factor * N;
        Kpp.noalias() += w * dNdx.transpose() * K_times_k_rel_over_mu * dNdx;
        // \TODO Extend to pressure dependent density.
        double const drhow_dp(0.0);
        Mpp.noalias() += (specific_storage * Sw + porosity * Sw * drhow_dp -
                          porosity * dSw_dpc) *
                         ip_data.mass_operator;

        if (_process_data.has_gravity)
        {
            Bp += w * density * dNdx.transpose() * K_times_k_rel_over_mu * b;
        }
        /* with Oberbeck-Boussing assumption density difference only exists
         * in buoyancy effects */
    }
}

template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
std::vector<double> const&
LocalAssemblerData<ShapeFunction, IntegrationMethod, GlobalDim>::
    getIntPtDarcyVelocity(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const
{
    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    constexpr int process_id = 0;  // monolithic scheme
    auto const indices =
        NumLib::getIndices(_element_id, *dof_table[process_id]);
    assert(!indices.empty());
    auto const local_x = x[process_id]->get(indices);

    cache.clear();
    auto cache_mat = MathLib::createZeroedMatrix<
        Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
        cache, GlobalDim, n_integration_points);

    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element_id);

    MaterialPropertyLib::VariableArray vars;

    // Get material properties
    auto const& medium = *_process_data.media_map->getMedium(_element_id);
    auto const& phase = medium.phase("AqueousLiquid");

    auto const p_nodal_values = Eigen::Map<const NodalVectorType>(
        &local_x[ShapeFunction::NPOINTS], ShapeFunction::NPOINTS);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        auto const& ip_data = _ip_data[ip];
        auto const& N = ip_data.N;
        auto const& dNdx = ip_data.dNdx;

        pos.setIntegrationPoint(ip);

        auto const dt = std::numeric_limits<double>::quiet_NaN();
        auto const K = MaterialPropertyLib::formEigenTensor<GlobalDim>(
            medium[MaterialPropertyLib::PropertyType::permeability].value(
                vars, pos, t, dt));
        auto const mu = phase[MaterialPropertyLib::PropertyType::viscosity]
                            .template value<double>(vars, pos, t, dt);

        double C_int_pt = 0.0;
        double p_int_pt = 0.0;
        NumLib::shapeFunctionInterpolate(local_x, N, C_int_pt, p_int_pt);

        // saturation
        vars[static_cast<int>(
            MaterialPropertyLib::Variable::capillary_pressure)] = -p_int_pt;
        auto const Sw = medium[MaterialPropertyLib::PropertyType::saturation]
                            .template value<double>(vars, pos, t, dt);

        vars[static_cast<int>(
            MaterialPropertyLib::Variable::liquid_saturation)] = Sw;
        auto const k_rel =
            medium[MaterialPropertyLib::PropertyType::relative_permeability]
                .template value<double>(vars, pos, t, dt);

        cache_mat.col(ip).noalias() = -dNdx * p_nodal_values;
        if (_process_data.has_gravity)
        {
            vars[static_cast<int>(
                MaterialPropertyLib::Variable::concentration)] = C_int_pt;
            vars[static_cast<int>(
                MaterialPropertyLib::Variable::phase_pressure)] = p_int_pt;
            auto const rho_w = phase[MaterialPropertyLib::PropertyType::density]
                                   .template value<double>(vars, pos, t, dt);
            auto const b = _process_data.specific_body_force;
            // here it is assumed that the vector b is directed 'downwards'
            cache_mat.col(ip).noalias() += rho_w * b;
        }
        cache_mat.col(ip).noalias() = k_rel / mu * (K * cache_mat.col(ip));
    }
    return cache;
}

template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
Eigen::Map<const Eigen::RowVectorXd>
LocalAssemblerData<ShapeFunction, IntegrationMethod, GlobalDim>::getShapeMatrix(
    const unsigned integration_point) const
{
    auto const& N = _ip_data[integration_point].N;

    // assumes N is stored contiguously in memory
    return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
}

template <typename ShapeFunction, typename IntegrationMethod, int GlobalDim>
std::vector<double> const&
LocalAssemblerData<ShapeFunction, IntegrationMethod, GlobalDim>::
    getIntPtSaturation(
        const double t,
        std::vector<GlobalVector*> const& x,
        std::vector<NumLib::LocalToGlobalIndexMap const*> const& dof_table,
        std::vector<double>& cache) const
{
    ParameterLib::SpatialPosition pos;
    pos.setElementID(_element_id);

    MaterialPropertyLib::VariableArray vars;

    auto const& medium = *_process_data.media_map->getMedium(_element_id);

    unsigned const n_integration_points =
        _integration_method.getNumberOfPoints();

    constexpr int process_id = 0;  // monolithic scheme
    auto const indices =
        NumLib::getIndices(_element_id, *dof_table[process_id]);
    assert(!indices.empty());
    auto const local_x = x[process_id]->get(indices);

    cache.clear();
    auto cache_vec = MathLib::createZeroedVector<LocalVectorType>(
        cache, n_integration_points);

    for (unsigned ip = 0; ip < n_integration_points; ++ip)
    {
        auto const& ip_data = _ip_data[ip];
        auto const& N = ip_data.N;

        double C_int_pt = 0.0;
        double p_int_pt = 0.0;
        NumLib::shapeFunctionInterpolate(local_x, N, C_int_pt, p_int_pt);

        // saturation
        vars[static_cast<int>(
            MaterialPropertyLib::Variable::capillary_pressure)] = -p_int_pt;
        auto const dt = std::numeric_limits<double>::quiet_NaN();
        auto const Sw = medium[MaterialPropertyLib::PropertyType::saturation]
                            .template value<double>(vars, pos, t, dt);
        cache_vec[ip] = Sw;
    }

    return cache;
}

}  // namespace RichardsComponentTransport
}  // namespace ProcessLib
