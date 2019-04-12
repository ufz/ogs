/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <Eigen/Dense>
#include <vector>

#include "ComponentTransportProcessData.h"
#include "MaterialLib/Fluid/FluidProperties/FluidProperties.h"
#include "MaterialLib/MPL/MaterialSpatialDistributionMap.h"
#include "MaterialLib/MPL/Medium.h"
#include "MaterialLib/MPL/Property.h"
#include "MathLib/LinAlg/Eigen/EigenMapTools.h"
#include "NumLib/DOF/DOFTableUtil.h"
#include "NumLib/Extrapolation/ExtrapolatableElement.h"
#include "NumLib/Fem/FiniteElement/TemplateIsoparametric.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "NumLib/Function/Interpolation.h"
#include "ParameterLib/Parameter.h"
#include "ProcessLib/CoupledSolutionsForStaggeredScheme.h"
#include "ProcessLib/LocalAssemblerInterface.h"
#include "ProcessLib/ProcessVariable.h"
#include "ProcessLib/Utils/InitShapeMatrices.h"

namespace ProcessLib
{
namespace ComponentTransport
{
template <typename NodalRowVectorType, typename GlobalDimNodalMatrixType>
struct IntegrationPointData final
{
    IntegrationPointData(NodalRowVectorType const& N_,
                         GlobalDimNodalMatrixType const& dNdx_,
                         double const& integration_weight_)
        : N(N_), dNdx(dNdx_), integration_weight(integration_weight_)
    {}

    NodalRowVectorType const N;
    GlobalDimNodalMatrixType const dNdx;
    double const integration_weight;

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
};

class ComponentTransportLocalAssemblerInterface
    : public ProcessLib::LocalAssemblerInterface,
      public NumLib::ExtrapolatableElement
{
public:
    ComponentTransportLocalAssemblerInterface() = default;

    void setStaggeredCoupledSolutions(
        std::size_t const /*mesh_item_id*/,
        CoupledSolutionsForStaggeredScheme* const coupling_term)
    {
        _coupled_solutions = coupling_term;
    }

    virtual std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        GlobalVector const& current_solution,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        std::vector<double>& cache) const = 0;

protected:
    CoupledSolutionsForStaggeredScheme* _coupled_solutions{nullptr};
};

template <typename ShapeFunction, typename IntegrationMethod,
          unsigned GlobalDim>
class LocalAssemblerData : public ComponentTransportLocalAssemblerInterface
{
    // When staggered scheme is adopted, nodal pressure and nodal concentration
    // are accessed by process id.
    static const int hydraulic_process_id = 0;
    // TODO (renchao-lu): This variable is used in the calculation of the
    // fluid's density and flux, indicating the transport process id. For now it
    // is assumed that these quantities depend on the first occurring transport
    // process only. The density and flux calculations have to be extended to
    // all processes.
    static const int first_transport_process_id = 1;

    // When monolithic scheme is adopted, nodal pressure and nodal concentration
    // are accessed by vector index.
    static const int pressure_index = 0;
    static const int first_concentration_index = ShapeFunction::NPOINTS;

    static const int pressure_size = ShapeFunction::NPOINTS;
    static const int concentration_size =
        ShapeFunction::NPOINTS;  // per component

    using ShapeMatricesType = ShapeMatrixPolicyType<ShapeFunction, GlobalDim>;
    using ShapeMatrices = typename ShapeMatricesType::ShapeMatrices;

    using LocalBlockMatrixType =
        typename ShapeMatricesType::template MatrixType<pressure_size,
                                                        pressure_size>;
    using LocalSegmentVectorType =
        typename ShapeMatricesType::template VectorType<pressure_size>;

    using LocalMatrixType =
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
    using LocalVectorType = Eigen::Matrix<double, Eigen::Dynamic, 1>;

    using NodalVectorType = typename ShapeMatricesType::NodalVectorType;
    using NodalRowVectorType = typename ShapeMatricesType::NodalRowVectorType;

    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;
    using GlobalDimNodalMatrixType =
        typename ShapeMatricesType::GlobalDimNodalMatrixType;
    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;

public:
    LocalAssemblerData(
        MeshLib::Element const& element,
        std::size_t const local_matrix_size,
        bool is_axially_symmetric,
        unsigned const integration_order,
        ComponentTransportProcessData const& process_data,
        std::vector<std::reference_wrapper<ProcessVariable>>
            transport_process_variables)
        : _element(element),
          _process_data(process_data),
          _integration_method(integration_order),
          _transport_process_variables(transport_process_variables)
    {
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

    void assemble(double const t, std::vector<double> const& local_x,
                  std::vector<double>& local_M_data,
                  std::vector<double>& local_K_data,
                  std::vector<double>& local_b_data) override
    {
        auto const local_matrix_size = local_x.size();
        // Nodal DOFs include pressure
        int const num_nodal_dof = 1 + _transport_process_variables.size();
        // This assertion is valid only if all nodal d.o.f. use the same shape
        // matrices.
        assert(local_matrix_size == ShapeFunction::NPOINTS * num_nodal_dof);

        auto local_M = MathLib::createZeroedMatrix<LocalMatrixType>(
            local_M_data, local_matrix_size, local_matrix_size);
        auto local_K = MathLib::createZeroedMatrix<LocalMatrixType>(
            local_K_data, local_matrix_size, local_matrix_size);
        auto local_b = MathLib::createZeroedVector<LocalVectorType>(
            local_b_data, local_matrix_size);

        // Get block matrices
        auto Kpp = local_K.template block<pressure_size, pressure_size>(
            pressure_index, pressure_index);
        auto Mpp = local_M.template block<pressure_size, pressure_size>(
            pressure_index, pressure_index);
        auto Bp = local_b.template segment<pressure_size>(pressure_index);

        auto local_p = Eigen::Map<const NodalVectorType>(
            &local_x[pressure_index], pressure_size);

        auto const number_of_components = num_nodal_dof - 1;
        for (int component_id = 0; component_id < number_of_components;
             ++component_id)
        {
            /*  Partitioned assembler matrix
             *  |  pp | pc1 | pc2 | pc3 |
             *  |-----|-----|-----|-----|
             *  | c1p | c1c1|  0  |  0  |
             *  |-----|-----|-----|-----|
             *  | c2p |  0  | c2c2|  0  |
             *  |-----|-----|-----|-----|
             *  | c3p |  0  |  0  | c3c3|
             */
            auto concentration_index =
                pressure_size + component_id * concentration_size;

            auto KCC =
                local_K.template block<concentration_size, concentration_size>(
                    concentration_index, concentration_index);
            auto MCC =
                local_M.template block<concentration_size, concentration_size>(
                    concentration_index, concentration_index);
            auto MCp =
                local_M.template block<concentration_size, pressure_size>(
                    concentration_index, pressure_index);
            auto MpC =
                local_M.template block<pressure_size, concentration_size>(
                    pressure_index, concentration_index);
            auto BC = local_b.template segment<concentration_size>(
                concentration_index);

            auto local_C = Eigen::Map<const NodalVectorType>(
                &local_x[concentration_index], concentration_size);

            assembleBlockMatrices(component_id, t, local_C, local_p, KCC, MCC,
                                  MCp, MpC, Kpp, Mpp, BC, Bp);
        }
    }

    void assembleBlockMatrices(
        int const component_id, double const t,
        Eigen::Ref<const NodalVectorType> const& C_nodal_values,
        Eigen::Ref<const NodalVectorType> const& p_nodal_values,
        Eigen::Ref<LocalBlockMatrixType> KCC,
        Eigen::Ref<LocalBlockMatrixType> MCC,
        Eigen::Ref<LocalBlockMatrixType> MCp,
        Eigen::Ref<LocalBlockMatrixType> MpC,
        Eigen::Ref<LocalBlockMatrixType> Kpp,
        Eigen::Ref<LocalBlockMatrixType> Mpp,
        Eigen::Ref<LocalSegmentVectorType> BC,
        Eigen::Ref<LocalSegmentVectorType> Bp)
    {
        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto const& b = _process_data.specific_body_force;

        MaterialLib::Fluid::FluidProperty::ArrayType vars;

        GlobalDimMatrixType const& I(
            GlobalDimMatrixType::Identity(GlobalDim, GlobalDim));

        // Get material properties
        auto const& medium =
            *_process_data.media_map->getMedium(_element.getID());
        // Select the only valid for component transport liquid phase.
        auto const& phase = medium.phase("AqueousLiquid");

        // Assume that the component name is the same as the process variable
        // name. Components are shifted by one because the first one is always
        // pressure.
        auto const& component = phase.component(
            _transport_process_variables[component_id].get().getName());

        for (unsigned ip(0); ip < n_integration_points; ++ip)
        {
            pos.setIntegrationPoint(ip);

            auto const& ip_data = _ip_data[ip];
            auto const& N = ip_data.N;
            auto const& dNdx = ip_data.dNdx;
            auto const& w = ip_data.integration_weight;

            double C_int_pt = 0.0;
            double p_int_pt = 0.0;

            NumLib::shapeFunctionInterpolate(C_nodal_values, N, C_int_pt);
            NumLib::shapeFunctionInterpolate(p_nodal_values, N, p_int_pt);

            // porosity model
            auto const porosity =
                _process_data.porous_media_properties.getPorosity(t, pos)
                    .getValue(t, pos, 0.0, C_int_pt);

            auto const retardation_factor =
                _process_data.retardation_factor(t, pos)[0];

            auto const& solute_dispersivity_transverse =
                medium.template value<double>(
                    MaterialPropertyLib::transversal_dispersivity);

            auto const& solute_dispersivity_longitudinal =
                medium.template value<double>(
                    MaterialPropertyLib::longitudinal_dispersivity);

            // Use the fluid density model to compute the density
            // TODO (renchao): concentration of which component as the argument
            // for calculation of fluid density
            vars[static_cast<int>(
                MaterialLib::Fluid::PropertyVariableType::C)] = C_int_pt;
            vars[static_cast<int>(
                MaterialLib::Fluid::PropertyVariableType::p)] = p_int_pt;
            auto const density = _process_data.fluid_properties->getValue(
                MaterialLib::Fluid::FluidPropertyType::Density, vars);
            auto const decay_rate = _process_data.decay_rate(t, pos)[0];

            auto const& molecular_diffusion_coefficient =
                component.template value<double>(
                    MaterialPropertyLib::molecular_diffusion);

            auto const& K =
                _process_data.porous_media_properties.getIntrinsicPermeability(
                    t, pos).getValue(t, pos, 0.0, 0.0);
            // Use the viscosity model to compute the viscosity
            auto const mu = _process_data.fluid_properties->getValue(
                MaterialLib::Fluid::FluidPropertyType::Viscosity, vars);

            GlobalDimMatrixType const K_over_mu = K / mu;
            GlobalDimVectorType const velocity =
                _process_data.has_gravity
                    ? GlobalDimVectorType(-K_over_mu *
                                          (dNdx * p_nodal_values - density * b))
                    : GlobalDimVectorType(-K_over_mu * dNdx * p_nodal_values);

            const double drho_dp = _process_data.fluid_properties->getdValue(
                MaterialLib::Fluid::FluidPropertyType::Density,
                vars,
                MaterialLib::Fluid::PropertyVariableType::p);

            const double drho_dC = _process_data.fluid_properties->getdValue(
                MaterialLib::Fluid::FluidPropertyType::Density,
                vars,
                MaterialLib::Fluid::PropertyVariableType::C);
            double const velocity_magnitude = velocity.norm();
            GlobalDimMatrixType const hydrodynamic_dispersion =
                velocity_magnitude != 0.0
                    ? GlobalDimMatrixType(
                          (porosity * molecular_diffusion_coefficient +
                           solute_dispersivity_transverse *
                               velocity_magnitude) *
                              I +
                          (solute_dispersivity_longitudinal -
                           solute_dispersivity_transverse) /
                              velocity_magnitude * velocity *
                              velocity.transpose())
                    : GlobalDimMatrixType(
                          (porosity * molecular_diffusion_coefficient +
                           solute_dispersivity_transverse *
                               velocity_magnitude) *
                          I);
            const double R_times_phi(retardation_factor * porosity);
            GlobalDimVectorType const mass_density_flow = velocity * density;
            if (_process_data.non_advective_form)
            {
                const double C_int_pt(N.dot(C_nodal_values));
                MCp.noalias() +=
                    N.transpose() * N * (C_int_pt * R_times_phi * drho_dp * w);
                MCC.noalias() +=
                    N.transpose() * N * (R_times_phi * C_int_pt * drho_dC * w);
                KCC.noalias() -= dNdx.transpose() * mass_density_flow * N * w;
            }
            else
            {
                BC.noalias() -=
                    N.transpose() *
                    (mass_density_flow.dot(dNdx * C_nodal_values) * w);
            }
            MCC.noalias() += N.transpose() * N * (R_times_phi * density * w);
            KCC.noalias() +=
                dNdx.transpose() * hydrodynamic_dispersion * dNdx *
                    (density * w) +
                N.transpose() * N * (decay_rate * R_times_phi * density * w);

            MpC.noalias() += N.transpose() * N * (porosity * drho_dC * w);

            // Calculate Mpp, Kpp, and bp in the first loop over components
            if (component_id == 0)
            {
                Mpp.noalias() += N.transpose() * N * (porosity * drho_dp * w);
                Kpp.noalias() +=
                    dNdx.transpose() * K_over_mu * dNdx * (density * w);

                if (_process_data.has_gravity)
                {
                    Bp.noalias() += dNdx.transpose() * K_over_mu * b *
                                    (density * density * w);
                }
            }
        }
    }

    void assembleForStaggeredScheme(
        double const t,
        std::vector<double>& local_M_data,
        std::vector<double>& local_K_data,
        std::vector<double>& local_b_data,
        LocalCoupledSolutions const& coupled_xs) override
    {
        if (coupled_xs.process_id == hydraulic_process_id)
        {
            assembleHydraulicEquation(t, local_M_data, local_K_data,
                                      local_b_data, coupled_xs);
        }
        else
        {
            // Go for assembling in an order of transport process id.
            assembleComponentTransportEquation(t, local_M_data, local_K_data,
                                               local_b_data, coupled_xs);
        }
    }

    void assembleHydraulicEquation(double const t,
                                   std::vector<double>& local_M_data,
                                   std::vector<double>& local_K_data,
                                   std::vector<double>& local_b_data,
                                   LocalCoupledSolutions const& coupled_xs)
    {
        auto local_p = Eigen::Map<const NodalVectorType>(
            coupled_xs.local_coupled_xs[hydraulic_process_id].data(),
            pressure_size);
        auto local_C = Eigen::Map<const NodalVectorType>(
            coupled_xs.local_coupled_xs[first_transport_process_id].data(),
            concentration_size);
        auto local_C0 = Eigen::Map<const NodalVectorType>(
            coupled_xs.local_coupled_xs0[first_transport_process_id].data(),
            concentration_size);

        auto const dt = coupled_xs.dt;

        auto local_M = MathLib::createZeroedMatrix<LocalBlockMatrixType>(
            local_M_data, pressure_size, pressure_size);
        auto local_K = MathLib::createZeroedMatrix<LocalBlockMatrixType>(
            local_K_data, pressure_size, pressure_size);
        auto local_b = MathLib::createZeroedVector<LocalSegmentVectorType>(
            local_b_data, pressure_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto const& b = _process_data.specific_body_force;

        MaterialLib::Fluid::FluidProperty::ArrayType vars;

        for (unsigned ip(0); ip < n_integration_points; ++ip)
        {
            pos.setIntegrationPoint(ip);

            auto const& ip_data = _ip_data[ip];
            auto const& N = ip_data.N;
            auto const& dNdx = ip_data.dNdx;
            auto const& w = ip_data.integration_weight;

            double C_int_pt = 0.0;
            double p_int_pt = 0.0;

            NumLib::shapeFunctionInterpolate(local_C, N, C_int_pt);
            NumLib::shapeFunctionInterpolate(local_p, N, p_int_pt);

            // porosity model
            auto const porosity =
                _process_data.porous_media_properties.getPorosity(t, pos)
                    .getValue(t, pos, 0.0, C_int_pt);

            // Use the fluid density model to compute the density
            // TODO: Concentration of which component as one of arguments for
            // calculation of fluid density
            vars[static_cast<int>(
                MaterialLib::Fluid::PropertyVariableType::C)] = C_int_pt;
            vars[static_cast<int>(
                MaterialLib::Fluid::PropertyVariableType::p)] = p_int_pt;
            auto const density = _process_data.fluid_properties->getValue(
                MaterialLib::Fluid::FluidPropertyType::Density, vars);

            auto const& K =
                _process_data.porous_media_properties.getIntrinsicPermeability(
                    t, pos).getValue(t, pos, 0.0, 0.0);
            // Use the viscosity model to compute the viscosity
            auto const mu = _process_data.fluid_properties->getValue(
                MaterialLib::Fluid::FluidPropertyType::Viscosity, vars);

            GlobalDimMatrixType const K_over_mu = K / mu;

            const double drho_dp = _process_data.fluid_properties->getdValue(
                MaterialLib::Fluid::FluidPropertyType::Density,
                vars,
                MaterialLib::Fluid::PropertyVariableType::p);
            const double drho_dC = _process_data.fluid_properties->getdValue(
                MaterialLib::Fluid::FluidPropertyType::Density,
                vars,
                MaterialLib::Fluid::PropertyVariableType::C);

            // matrix assembly
            local_M.noalias() += w * N.transpose() * porosity * drho_dp * N;
            local_K.noalias() +=
                w * dNdx.transpose() * density * K_over_mu * dNdx;

            if (_process_data.has_gravity)
                local_b.noalias() +=
                    w * density * density * dNdx.transpose() * K_over_mu * b;

            // coupling term
            {
                double C0_int_pt = 0.0;
                NumLib::shapeFunctionInterpolate(local_C0, N, C0_int_pt);

                local_b.noalias() -= w * N.transpose() * porosity * drho_dC *
                                     (C_int_pt - C0_int_pt) / dt;
            }
        }
    }

    void assembleComponentTransportEquation(
        double const t, std::vector<double>& local_M_data,
        std::vector<double>& local_K_data, std::vector<double>& local_b_data,
        LocalCoupledSolutions const& coupled_xs)
    {
        auto const& transport_process_id = coupled_xs.process_id;
        auto local_C = Eigen::Map<const NodalVectorType>(
            coupled_xs.local_coupled_xs[transport_process_id].data(),
            concentration_size);
        auto local_p = Eigen::Map<const NodalVectorType>(
            coupled_xs.local_coupled_xs[hydraulic_process_id].data(),
            pressure_size);
        auto local_p0 = Eigen::Map<const NodalVectorType>(
            coupled_xs.local_coupled_xs0[hydraulic_process_id].data(),
            pressure_size);

        auto local_M = MathLib::createZeroedMatrix<LocalBlockMatrixType>(
            local_M_data, concentration_size, concentration_size);
        auto local_K = MathLib::createZeroedMatrix<LocalBlockMatrixType>(
            local_K_data, concentration_size, concentration_size);

        auto const dt = coupled_xs.dt;
        auto local_b = MathLib::createZeroedVector<LocalSegmentVectorType>(
            local_b_data, concentration_size);

        unsigned const n_integration_points =
            _integration_method.getNumberOfPoints();

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        auto const& b = _process_data.specific_body_force;

        MaterialLib::Fluid::FluidProperty::ArrayType vars;

        GlobalDimMatrixType const& I(
            GlobalDimMatrixType::Identity(GlobalDim, GlobalDim));

        auto const& medium =
            *_process_data.media_map->getMedium(_element.getID());
        auto const& phase = medium.phase("AqueousLiquid");
        // Hydraulic process id is 0 and thus transport process id starts
        // from 1.
        auto const component_id = transport_process_id - 1;
        auto const& component = phase.component(
            _transport_process_variables[component_id].get().getName());

        for (unsigned ip(0); ip < n_integration_points; ++ip)
        {
            pos.setIntegrationPoint(ip);

            auto const& ip_data = _ip_data[ip];
            auto const& N = ip_data.N;
            auto const& dNdx = ip_data.dNdx;
            auto const& w = ip_data.integration_weight;

            double C_int_pt = 0.0;
            double p_int_pt = 0.0;

            NumLib::shapeFunctionInterpolate(local_C, N, C_int_pt);
            NumLib::shapeFunctionInterpolate(local_p, N, p_int_pt);

            // porosity model
            auto const porosity =
                _process_data.porous_media_properties.getPorosity(t, pos)
                    .getValue(t, pos, 0.0, C_int_pt);

            auto const retardation_factor =
                _process_data.retardation_factor(t, pos)[0];

            auto const& solute_dispersivity_transverse =
                medium.template value<double>(
                    MaterialPropertyLib::transversal_dispersivity);
            auto const& solute_dispersivity_longitudinal =
                medium.template value<double>(
                    MaterialPropertyLib::longitudinal_dispersivity);

            // Use the fluid density model to compute the density
            vars[static_cast<int>(
                MaterialLib::Fluid::PropertyVariableType::C)] = C_int_pt;
            vars[static_cast<int>(
                MaterialLib::Fluid::PropertyVariableType::p)] = p_int_pt;
            auto const density = _process_data.fluid_properties->getValue(
                MaterialLib::Fluid::FluidPropertyType::Density, vars);
            auto const decay_rate = _process_data.decay_rate(t, pos)[0];

            auto const& molecular_diffusion_coefficient =
                component.template value<double>(
                    MaterialPropertyLib::molecular_diffusion);

            auto const& K =
                _process_data.porous_media_properties.getIntrinsicPermeability(
                    t, pos).getValue(t, pos, 0.0, 0.0);
            // Use the viscosity model to compute the viscosity
            auto const mu = _process_data.fluid_properties->getValue(
                MaterialLib::Fluid::FluidPropertyType::Viscosity, vars);

            GlobalDimMatrixType const K_over_mu = K / mu;
            GlobalDimVectorType const velocity =
                _process_data.has_gravity
                    ? GlobalDimVectorType(-K_over_mu *
                                          (dNdx * local_p - density * b))
                    : GlobalDimVectorType(-K_over_mu * dNdx * local_p);


            double const velocity_magnitude = velocity.norm();
            GlobalDimMatrixType const hydrodynamic_dispersion =
                velocity_magnitude != 0.0
                    ? GlobalDimMatrixType(
                          (porosity * molecular_diffusion_coefficient +
                           solute_dispersivity_transverse *
                               velocity_magnitude) *
                              I +
                          (solute_dispersivity_longitudinal -
                           solute_dispersivity_transverse) /
                              velocity_magnitude * velocity *
                              velocity.transpose())
                    : GlobalDimMatrixType(
                          (porosity * molecular_diffusion_coefficient +
                           solute_dispersivity_transverse *
                               velocity_magnitude) *
                          I);

            double const R_times_phi = retardation_factor * porosity;

            if (_process_data.non_advective_form)
            {
                const double drho_dC =
                    _process_data.fluid_properties->getdValue(
                        MaterialLib::Fluid::FluidPropertyType::Density,
                        vars,
                        MaterialLib::Fluid::PropertyVariableType::C);
                local_M.noalias() +=
                    N.transpose() * N * (R_times_phi * C_int_pt * drho_dC * w);
            }
            else
            {
                local_b.noalias() -= N.transpose() *
                                     velocity.dot(dNdx * local_C) *
                                     (density * w);
            }
            local_M.noalias() +=
                N.transpose() * N * (R_times_phi * density * w);

            // coupling term

            if (_process_data.non_advective_form)
            {
                double p0_int_pt = 0.0;

                NumLib::shapeFunctionInterpolate(local_p0, N, p0_int_pt);
                const double drho_dp =
                    _process_data.fluid_properties->getdValue(
                        MaterialLib::Fluid::FluidPropertyType::Density,
                        vars,
                        MaterialLib::Fluid::PropertyVariableType::p);
                local_K.noalias() +=
                    N.transpose() * N *
                        ((R_times_phi * drho_dp * (p_int_pt - p0_int_pt) / dt) *
                         w) -
                    dNdx.transpose() * velocity * N * (density * w);
            }
            local_K.noalias() +=
                dNdx.transpose() * hydrodynamic_dispersion * dNdx *
                    (density * w) +
                N.transpose() * N * (decay_rate * R_times_phi * density * w);
        }
    }

    std::vector<double> const& getIntPtDarcyVelocity(
        const double t,
        GlobalVector const& current_solution,
        NumLib::LocalToGlobalIndexMap const& dof_table,
        std::vector<double>& cache) const override
    {
        auto const indices = NumLib::getIndices(_element.getID(), dof_table);
        assert(!indices.empty());

        if (_coupled_solutions == nullptr)  // monolithic scheme
        {
            auto const local_x = current_solution.get(indices);

            // Assuming that fluid density always depends on the concentration
            // of the component which is placed at the uppermost.
            auto const local_p = Eigen::Map<const NodalVectorType>(
                &local_x[pressure_index], pressure_size);
            auto const local_C = Eigen::Map<const NodalVectorType>(
                &local_x[first_concentration_index], concentration_size);

            return calculateIntPtDarcyVelocity(t, local_p, local_C, cache);
        }
        else  // staggered scheme
        {
            std::vector<std::vector<GlobalIndexType>>
                indices_of_all_coupled_processes(
                    _coupled_solutions->coupled_xs.size(), indices);

            auto const local_xs = getCurrentLocalSolutions(
                *(this->_coupled_solutions), indices_of_all_coupled_processes);

            auto const local_p = MathLib::toVector(local_xs[hydraulic_process_id]);
            auto const local_C =
                MathLib::toVector(local_xs[first_transport_process_id]);

            return calculateIntPtDarcyVelocity(t, local_p, local_C, cache);
        }
    }

    std::vector<double> const& calculateIntPtDarcyVelocity(
        const double t,
        Eigen::Ref<const NodalVectorType> const& p_nodal_values,
        Eigen::Ref<const NodalVectorType> const& C_nodal_values,
        std::vector<double>& cache) const
    {
        auto const n_integration_points =
            _integration_method.getNumberOfPoints();

        cache.clear();
        auto cache_mat = MathLib::createZeroedMatrix<
            Eigen::Matrix<double, GlobalDim, Eigen::Dynamic, Eigen::RowMajor>>(
            cache, GlobalDim, n_integration_points);

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        MaterialLib::Fluid::FluidProperty::ArrayType vars;

        for (unsigned ip = 0; ip < n_integration_points; ++ip)
        {
            auto const& ip_data = _ip_data[ip];
            auto const& N = ip_data.N;
            auto const& dNdx = ip_data.dNdx;

            pos.setIntegrationPoint(ip);

            auto const& K =
                _process_data.porous_media_properties.getIntrinsicPermeability(
                    t, pos).getValue(t, pos, 0.0, 0.0);
            auto const mu = _process_data.fluid_properties->getValue(
                MaterialLib::Fluid::FluidPropertyType::Viscosity, vars);
            GlobalDimMatrixType const K_over_mu = K / mu;

            cache_mat.col(ip).noalias() = -K_over_mu * dNdx * p_nodal_values;
            if (_process_data.has_gravity)
            {
                double C_int_pt = 0.0;
                double p_int_pt = 0.0;

                NumLib::shapeFunctionInterpolate(C_nodal_values, N, C_int_pt);
                NumLib::shapeFunctionInterpolate(p_nodal_values, N, p_int_pt);

                vars[static_cast<int>(
                    MaterialLib::Fluid::PropertyVariableType::C)] = C_int_pt;
                vars[static_cast<int>(
                    MaterialLib::Fluid::PropertyVariableType::p)] = p_int_pt;

                auto const rho_w = _process_data.fluid_properties->getValue(
                    MaterialLib::Fluid::FluidPropertyType::Density, vars);
                auto const b = _process_data.specific_body_force;
                // here it is assumed that the vector b is directed 'downwards'
                cache_mat.col(ip).noalias() += K_over_mu * rho_w * b;
            }
        }

        return cache;
    }

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N = _ip_data[integration_point].N;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N.data(), N.size());
    }

    Eigen::Vector3d getFlux(MathLib::Point3d const& pnt_local_coords,
                            double const t,
                            std::vector<double> const& local_x) const override
    {
        auto const local_p = Eigen::Map<const NodalVectorType>(
            &local_x[pressure_index], pressure_size);
        auto const local_C = Eigen::Map<const NodalVectorType>(
            &local_x[first_concentration_index], concentration_size);

        return calculateFlux(pnt_local_coords, t, local_p, local_C);
    }

    Eigen::Vector3d getFlux(
        MathLib::Point3d const& pnt_local_coords,
        double const t,
        std::vector<std::vector<double>> const& local_xs) const override
    {
        auto const local_p = MathLib::toVector(local_xs[hydraulic_process_id]);
        auto const local_C =
            MathLib::toVector(local_xs[first_transport_process_id]);

        return calculateFlux(pnt_local_coords, t, local_p, local_C);
    }

    Eigen::Vector3d calculateFlux(
        MathLib::Point3d const& pnt_local_coords,
        double const t,
        Eigen::Ref<const NodalVectorType> const& p_nodal_values,
        Eigen::Ref<const NodalVectorType> const& C_nodal_values) const
    {
        // eval shape matrices at given point
        auto const shape_matrices = [&]() {
            using FemType =
                NumLib::TemplateIsoparametric<ShapeFunction, ShapeMatricesType>;

            FemType fe(*static_cast<const typename ShapeFunction::MeshElement*>(
                &_element));

            typename ShapeMatricesType::ShapeMatrices shape_matrices(
                ShapeFunction::DIM, GlobalDim, ShapeFunction::NPOINTS);

            // Note: Axial symmetry is set to false here, because we only need
            // dNdx here, which is not affected by axial symmetry.
            fe.template computeShapeFunctions<NumLib::ShapeMatrixType::DNDX>(
                pnt_local_coords.getCoords(), shape_matrices, GlobalDim, false);

            return shape_matrices;
        }();

        ParameterLib::SpatialPosition pos;
        pos.setElementID(_element.getID());

        MaterialLib::Fluid::FluidProperty::ArrayType vars;

        // local_x contains the local concentration and pressure values
        NumLib::shapeFunctionInterpolate(
            C_nodal_values, shape_matrices.N,
            vars[static_cast<int>(
                MaterialLib::Fluid::PropertyVariableType::C)]);
        NumLib::shapeFunctionInterpolate(
            p_nodal_values, shape_matrices.N,
            vars[static_cast<int>(
                MaterialLib::Fluid::PropertyVariableType::p)]);

        auto const K =
            _process_data.porous_media_properties
                .getIntrinsicPermeability(t, pos)
                .getValue(t, pos, 0.0,
                          vars[static_cast<int>(
                              MaterialLib::Fluid::PropertyVariableType::C)]);

        auto const mu = _process_data.fluid_properties->getValue(
            MaterialLib::Fluid::FluidPropertyType::Viscosity, vars);
        GlobalDimMatrixType const K_over_mu = K / mu;

        GlobalDimVectorType q =
            -K_over_mu * shape_matrices.dNdx * p_nodal_values;

        if (_process_data.has_gravity)
        {
            auto const rho_w =
                _process_data.fluid_properties->getValue(
                    MaterialLib::Fluid::FluidPropertyType::Density, vars);
            auto const b = _process_data.specific_body_force;
            q += K_over_mu * rho_w * b;
        }
        Eigen::Vector3d flux(0.0, 0.0, 0.0);
        flux.head<GlobalDim>() = q;
        return flux;
    }

private:
    MeshLib::Element const& _element;
    ComponentTransportProcessData const& _process_data;

    IntegrationMethod const _integration_method;
    std::vector<std::reference_wrapper<ProcessVariable>> const
        _transport_process_variables;

    std::vector<
        IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType>,
        Eigen::aligned_allocator<
            IntegrationPointData<NodalRowVectorType, GlobalDimNodalMatrixType>>>
        _ip_data;
};

}  // namespace ComponentTransport
}  // namespace ProcessLib
