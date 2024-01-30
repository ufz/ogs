/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include "Base.h"
#include "MaterialLib/SolidModels/MFront/MFrontGeneric.h"
#include "MaterialLib/SolidModels/MFront/Variable.h"
#include "MaterialState.h"
#include "ProcessLib/ConstitutiveRelations/StrainData.h"
#include "ProcessLib/ConstitutiveRelations/StressData.h"

namespace ProcessLib::LargeDeformation
{
namespace ConstitutiveRelations
{
template <int DisplacementDim>
struct SolidMechanicsDataStateless
{
    KelvinMatrix<DisplacementDim> stiffness_tensor =
        KV::KMnan<DisplacementDim>();
    KelvinVector<DisplacementDim> J_uT_BT_K_N = KV::KVnan<DisplacementDim>();
    KelvinVector<DisplacementDim> J_up_BT_K_N = KV::KVnan<DisplacementDim>();
};

namespace MSM = MaterialLib::Solids::MFront;

template <int DisplacementDim>
using SolidConstitutiveRelation =
    MSM::MFrontGeneric<DisplacementDim,
                       boost::mp11::mp_list<MSM::DeformationGradient>,
                       boost::mp11::mp_list<MSM::SecondPiolaKirchhoffStress>,
                       boost::mp11::mp_list<MSM::Temperature>>;

template <int DisplacementDim>
struct SolidMechanicsModel
{
    explicit SolidMechanicsModel(
        SolidConstitutiveRelation<DisplacementDim> const& solid_material)
        : solid_material_(solid_material),
          tangent_operator_blocks_view_{
              solid_material.template createTangentOperatorBlocksView<
                  typename MSM::ForcesGradsCombinations<
                      boost::mp11::mp_list<MSM::GreenLagrangeStrain>,
                      boost::mp11::mp_list<MSM::SecondPiolaKirchhoffStress>,
                      boost::mp11::mp_list<MSM::Temperature>>::type>()}
    {
    }

    void eval(
        SpaceTimeData const& x_t,
        Temperature const& temperature,
        DeformationGradientData<DisplacementDim> const&
            deformation_gradient_data,
        PrevState<DeformationGradientData<DisplacementDim>> const&
            deformation_gradient_data_prev,
        MaterialStateData<DisplacementDim>& mat_state,
        PrevState<StressData<DisplacementDim>> const& stress_data_prev,
        StressData<DisplacementDim>& stress_data,
        SolidMechanicsDataStateless<DisplacementDim>& current_stateless) const;

    auto getInternalVariables() const
    {
        return solid_material_.getInternalVariables();
    }

private:
    SolidConstitutiveRelation<DisplacementDim> const& solid_material_;

    MSM::OGSMFrontTangentOperatorBlocksView<
        DisplacementDim,
        MSM::ForcesGradsCombinations<
            boost::mp11::mp_list<MSM::GreenLagrangeStrain>,
            boost::mp11::mp_list<MSM::SecondPiolaKirchhoffStress>,
            boost::mp11::mp_list<MSM::Temperature>>::type>
        tangent_operator_blocks_view_;
};

extern template struct SolidMechanicsModel<2>;
extern template struct SolidMechanicsModel<3>;
}  // namespace ConstitutiveRelations
}  // namespace ProcessLib::LargeDeformation
