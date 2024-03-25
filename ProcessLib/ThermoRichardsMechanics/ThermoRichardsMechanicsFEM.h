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
#include <vector>

#include "IntegrationPointData.h"
#include "LocalAssemblerInterface.h"
#include "MathLib/EigenBlockMatrixView.h"
#include "MathLib/KelvinVector.h"
#include "NumLib/DOF/LocalDOF.h"
#include "NumLib/Fem/InitShapeMatrices.h"
#include "NumLib/Fem/Integration/GenericIntegrationMethod.h"
#include "NumLib/Fem/ShapeMatrixPolicy.h"
#include "ProcessLib/Deformation/BMatrixPolicy.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/Porosity.h"
#include "ProcessLib/ThermoRichardsMechanics/ConstitutiveCommon/TransportPorosity.h"
#include "ThermoRichardsMechanicsProcessData.h"

namespace ProcessLib
{
namespace ThermoRichardsMechanics
{
namespace MPL = MaterialPropertyLib;

template <typename ShapeFunctionDisplacement, typename ShapeFunction,
          int DisplacementDim, typename ConstitutiveTraits>
class ThermoRichardsMechanicsLocalAssembler
    : public LocalAssemblerInterface<DisplacementDim, ConstitutiveTraits>
{
    static constexpr int temperature_index = 0;
    static constexpr int temperature_size = ShapeFunction::NPOINTS;
    static constexpr int pressure_index = temperature_size;
    static constexpr int pressure_size = ShapeFunction::NPOINTS;
    static constexpr int displacement_index = 2 * ShapeFunction::NPOINTS;
    static constexpr int displacement_size =
        ShapeFunctionDisplacement::NPOINTS * DisplacementDim;

public:
    using ShapeMatricesTypeDisplacement =
        ShapeMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;
    // Note: temperature variable uses the same shape functions as that are used
    // by pressure variable.
    using ShapeMatricesType =
        ShapeMatrixPolicyType<ShapeFunction, DisplacementDim>;

    using GlobalDimMatrixType = typename ShapeMatricesType::GlobalDimMatrixType;
    using GlobalDimVectorType = typename ShapeMatricesType::GlobalDimVectorType;

    using BMatricesType =
        BMatrixPolicyType<ShapeFunctionDisplacement, DisplacementDim>;
    using KelvinVectorType = typename BMatricesType::KelvinVectorType;

    using IpData = IntegrationPointData<ShapeMatricesTypeDisplacement,
                                        ShapeMatricesType, DisplacementDim,
                                        ShapeFunctionDisplacement::NPOINTS>;

    static int const KelvinVectorSize =
        MathLib::KelvinVector::kelvin_vector_dimensions(DisplacementDim);
    using Invariants = MathLib::KelvinVector::Invariants<KelvinVectorSize>;

    using SymmetricTensor = Eigen::Matrix<double, KelvinVectorSize, 1>;

    static constexpr auto& N_u_op = MathLib::eigenBlockMatrixView<
        DisplacementDim,
        typename ShapeMatricesTypeDisplacement::NodalRowVectorType>;

    ThermoRichardsMechanicsLocalAssembler(
        ThermoRichardsMechanicsLocalAssembler const&) = delete;
    ThermoRichardsMechanicsLocalAssembler(
        ThermoRichardsMechanicsLocalAssembler&&) = delete;

    ThermoRichardsMechanicsLocalAssembler(
        MeshLib::Element const& e,
        std::size_t const /*local_matrix_size*/,
        NumLib::GenericIntegrationMethod const& integration_method,
        bool const is_axially_symmetric,
        ThermoRichardsMechanicsProcessData<DisplacementDim, ConstitutiveTraits>&
            process_data);

    void setInitialConditionsConcrete(Eigen::VectorXd const local_x,
                                      double const t,
                                      int const process_id) override;

    class LocalMatrices
    {
        using NodalMatrix = typename ShapeMatricesType::NodalMatrixType;

        static auto constexpr local_matrix_dim =
            displacement_size + pressure_size + temperature_size;

        template <Eigen::Index rows, Eigen::Index cols>
        using Mat =
            typename ShapeMatricesTypeDisplacement::template MatrixType<rows,
                                                                        cols>;
        using Vec = typename ShapeMatricesTypeDisplacement::template VectorType<
            local_matrix_dim>;

    public:
        void setZero()
        {
            M_TT = NodalMatrix::Zero(temperature_size, temperature_size);
            M_Tp = NodalMatrix::Zero(temperature_size, pressure_size);
            K_TT = NodalMatrix::Zero(temperature_size, temperature_size);
            K_Tp = NodalMatrix::Zero(temperature_size, pressure_size);
            dK_TT_dp = NodalMatrix::Zero(temperature_size, pressure_size);

            M_pu = Mat<pressure_size, displacement_size>::Zero(
                pressure_size, displacement_size);

            M_pT = NodalMatrix::Zero(pressure_size, temperature_size);

            K_pp = NodalMatrix::Zero(pressure_size, pressure_size);
            K_pT = NodalMatrix::Zero(pressure_size, temperature_size);

            storage_p_a_p = NodalMatrix::Zero(pressure_size, pressure_size);
            storage_p_a_S_Jpp = NodalMatrix::Zero(pressure_size, pressure_size);
            storage_p_a_S = NodalMatrix::Zero(pressure_size, pressure_size);

            Jac = Mat<local_matrix_dim, local_matrix_dim>::Zero(
                local_matrix_dim, local_matrix_dim);
            res = Vec::Zero(local_matrix_dim);
        }

        LocalMatrices& operator+=(LocalMatrices const& other)
        {
            M_TT += other.M_TT;
            M_Tp += other.M_Tp;
            K_TT += other.K_TT;
            K_Tp += other.K_Tp;
            dK_TT_dp += other.dK_TT_dp;

            M_pu += other.M_pu;

            M_pT += other.M_pT;

            K_pp += other.K_pp;
            K_pT += other.K_pT;

            storage_p_a_p += other.storage_p_a_p;
            storage_p_a_S_Jpp += other.storage_p_a_S_Jpp;
            storage_p_a_S += other.storage_p_a_S;

            Jac += other.Jac;
            res += other.res;

            return *this;
        }

        LocalMatrices& operator*=(double const a)
        {
            M_TT *= a;
            M_Tp *= a;
            K_TT *= a;
            K_Tp *= a;
            dK_TT_dp *= a;

            M_pu *= a;

            M_pT *= a;

            K_pp *= a;
            K_pT *= a;

            storage_p_a_p *= a;
            storage_p_a_S_Jpp *= a;
            storage_p_a_S *= a;

            Jac *= a;
            res *= a;

            return *this;
        }

        template <typename OStream>
        friend OStream& operator<<(OStream& os, LocalMatrices const& loc_mat)
        {
            os << "- M_TT:\n"
               << loc_mat.M_TT  //
               << "\n- M_Tp:\n"
               << loc_mat.M_Tp  //
               << "\n- K_TT:\n"
               << loc_mat.K_TT  //
               << "\n- K_Tp:\n"
               << loc_mat.K_Tp  //
               << "\n- dK_TT_dp:\n"
               << loc_mat.dK_TT_dp  //
               << "\n- M_pu:\n"
               << loc_mat.M_pu  //
               << "\n- M_pT:\n"
               << loc_mat.M_pT  //
               << "\n- K_pp:\n"
               << loc_mat.K_pp  //
               << "\n- K_pT:\n"
               << loc_mat.K_pT  //
               << "\n- storage_p_a_p:\n"
               << loc_mat.storage_p_a_p  //
               << "\n- storage_p_a_S_Jpp:\n"
               << loc_mat.storage_p_a_S_Jpp  //
               << "\n- storage_p_a_S:\n"
               << loc_mat.storage_p_a_S  //
               << "\n- Jac:\n"
               << loc_mat.Jac  //
               << "\n- res:\n"
               << loc_mat.res;
            return os;
        }

        NodalMatrix M_TT;
        NodalMatrix M_Tp;
        NodalMatrix K_TT;
        NodalMatrix K_Tp;
        NodalMatrix dK_TT_dp;

        Mat<pressure_size, displacement_size> M_pu;

        NodalMatrix M_pT;

        NodalMatrix K_pp;
        NodalMatrix K_pT;

        NodalMatrix storage_p_a_p;
        NodalMatrix storage_p_a_S_Jpp;
        NodalMatrix storage_p_a_S;

        //! "Direct" contributions to the Jacobian, without those from K and M
        //! matrices.
        Mat<local_matrix_dim, local_matrix_dim> Jac;

        //! "Direct" contributions to the residual, without those from K and M
        //! matrices.
        Vec res;
    };

    void assembleWithJacobian(double const t, double const dt,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_x_prev,
                              std::vector<double>& /*local_M_data*/,
                              std::vector<double>& /*local_K_data*/,
                              std::vector<double>& local_rhs_data,
                              std::vector<double>& local_Jac_data) override;

private:
    void assembleWithJacobianSingleIP(
        double const t, double const dt,
        ParameterLib::SpatialPosition const& x_position,
        std::vector<double> const& local_x,
        std::vector<double> const& local_x_prev, IpData const& ip_data,
        typename ConstitutiveTraits::ConstitutiveSetting& CS,
        MaterialPropertyLib::Medium& medium, LocalMatrices& out,
        typename ConstitutiveTraits::StatefulData& current_state,
        typename ConstitutiveTraits::StatefulDataPrev const& prev_state,
        MaterialStateData<DisplacementDim>& mat_state,
        typename ConstitutiveTraits::OutputData& output_data) const;

    void addToLocalMatrixData(double const dt,
                              std::vector<double> const& local_x,
                              std::vector<double> const& local_x_prev,
                              LocalMatrices const& loc_mat,
                              std::vector<double>& local_rhs_data,
                              std::vector<double>& local_Jac_data) const;

    void massLumping(LocalMatrices& loc_mat) const;

public:
    void initializeConcrete() override
    {
        unsigned const n_integration_points =
            this->integration_method_.getNumberOfPoints();
        auto const time_independent = std::numeric_limits<double>::quiet_NaN();
        auto const& medium =
            *this->process_data_.media_map.getMedium(this->element_.getID());

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            ParameterLib::SpatialPosition const x_position{
                std::nullopt, this->element_.getID(), ip,
                MathLib::Point3d(NumLib::interpolateCoordinates<
                                 ShapeFunctionDisplacement,
                                 ShapeMatricesTypeDisplacement>(
                    this->element_, this->ip_data_[ip].N_u))};
            auto& current_state = this->current_states_[ip];

            // Set initial stress from parameter.
            if (this->process_data_.initial_stress.value)
            {
                ConstitutiveTraits::ConstitutiveSetting::statefulStress(
                    current_state) =
                    MathLib::KelvinVector::symmetricTensorToKelvinVector<
                        DisplacementDim>(
                        (*this->process_data_.initial_stress.value)(
                            time_independent, x_position));
            }

            if (*this->process_data_.initialize_porosity_from_medium_property)
            {
                // Initial porosity. Could be read from integration point data
                // or mesh.
                std::get<PorosityData>(current_state).phi =
                    medium.property(MPL::porosity)
                        .template initialValue<double>(x_position,
                                                       time_independent);

                if (medium.hasProperty(MPL::PropertyType::transport_porosity))
                {
                    std::get<TransportPorosityData>(current_state).phi =
                        medium.property(MPL::transport_porosity)
                            .template initialValue<double>(x_position,
                                                           time_independent);
                }
                else
                {
                    std::get<TransportPorosityData>(current_state).phi =
                        std::get<PorosityData>(current_state).phi;
                }
            }

            double const t = 0;  // TODO (naumov) pass t from top
            this->solid_material_.initializeInternalStateVariables(
                t, x_position,
                *this->material_states_[ip].material_state_variables);
        }

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            this->material_states_[ip].pushBackState();
        }

        for (unsigned ip = 0; ip < n_integration_points; ip++)
        {
            this->prev_states_[ip] = this->current_states_[ip];
        }
    }

    void computeSecondaryVariableConcrete(
        double const t, double const dt, Eigen::VectorXd const& local_x,
        Eigen::VectorXd const& local_x_prev) override;

    Eigen::Map<const Eigen::RowVectorXd> getShapeMatrix(
        const unsigned integration_point) const override
    {
        auto const& N_u = ip_data_[integration_point].N_u;

        // assumes N is stored contiguously in memory
        return Eigen::Map<const Eigen::RowVectorXd>(N_u.data(), N_u.size());
    }

private:
    std::vector<IpData> ip_data_;

    static constexpr auto localDOF(std::vector<double> const& x)
    {
        return NumLib::localDOF<
            ShapeFunction, ShapeFunction,
            NumLib::Vectorial<ShapeFunctionDisplacement, DisplacementDim>>(x);
    }

    static auto block_uu(auto& mat)
    {
        return mat.template block<displacement_size, displacement_size>(
            displacement_index, displacement_index);
    }
    static auto block_up(auto& mat)
    {
        return mat.template block<displacement_size, pressure_size>(
            displacement_index, pressure_index);
    }
    static auto block_uT(auto& mat)
    {
        return mat.template block<displacement_size, temperature_size>(
            displacement_index, temperature_index);
    }
    static auto block_pu(auto& mat)
    {
        return mat.template block<pressure_size, displacement_size>(
            pressure_index, displacement_index);
    }
    static auto block_pp(auto& mat)
    {
        return mat.template block<pressure_size, pressure_size>(pressure_index,
                                                                pressure_index);
    }
    static auto block_pT(auto& mat)
    {
        return mat.template block<pressure_size, temperature_size>(
            pressure_index, temperature_index);
    }
    static auto block_Tp(auto& mat)
    {
        return mat.template block<temperature_size, pressure_size>(
            temperature_index, pressure_index);
    }
    static auto block_TT(auto& mat)
    {
        return mat.template block<temperature_size, temperature_size>(
            temperature_index, temperature_index);
    }

    static auto block_u(auto& vec)
    {
        return vec.template segment<displacement_size>(displacement_index);
    }
    static auto block_p(auto& vec)
    {
        return vec.template segment<pressure_size>(pressure_index);
    }
    static auto block_T(auto& vec)
    {
        return vec.template segment<temperature_size>(temperature_index);
    }

    /// This function resets the initial stress type according to the input
    /// initial stress type, either total or effective.
    /// If subtype = `StressSaturation_StrainPressureTemperature` is bing used
    /// in the process setting, the initial effective stress is converted to
    /// total stress. Otherwise, the initial total stress is converted to
    /// effective stress.
    void convertInitialStressType(
        unsigned const ip, double const t,
        ParameterLib::SpatialPosition const x_position,
        MaterialPropertyLib::Medium const& medium,
        MPL::VariableArray const& variables, double const p_at_ip);
};

}  // namespace ThermoRichardsMechanics
}  // namespace ProcessLib

#include "ThermoRichardsMechanicsFEM-impl.h"
