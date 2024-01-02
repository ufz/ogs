/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "CreateThermoRichardsMechanicsLocalAssemblers.h"

#include "ConstitutiveStress_StrainTemperature/Traits.h"
#ifdef OGS_USE_MFRONT
#include "ConstitutiveStressSaturation_StrainPressureTemperature/Traits.h"
#endif

#include "ProcessLib/Utils/CreateLocalAssemblersTaylorHood.h"
#include "ThermoRichardsMechanicsFEM.h"

namespace ProcessLib::ThermoRichardsMechanics
{
template <typename ConstitutiveTraits>
struct ThermoRichardsMechanicsLocalAssembler3Args
{
    template <typename ShapeFunctionDisplacement, typename ShapeFunction,
              int DisplacementDim>
    using LocalAssemblerImplementation =
        ThermoRichardsMechanicsLocalAssembler<ShapeFunctionDisplacement,
                                              ShapeFunction, DisplacementDim,
                                              ConstitutiveTraits>;
};

template <int DisplacementDim, typename ConstitutiveTraits>
void createLocalAssemblers(
    std::vector<MeshLib::Element*> const& mesh_elements,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::unique_ptr<
        LocalAssemblerInterface<DisplacementDim, ConstitutiveTraits>>>&
        local_assemblers,
    NumLib::IntegrationOrder const integration_order,
    bool const is_axially_symmetric,
    ThermoRichardsMechanicsProcessData<DisplacementDim, ConstitutiveTraits>&
        process_data)
{
    createLocalAssemblersHM<
        DisplacementDim,
        ThermoRichardsMechanicsLocalAssembler3Args<
            ConstitutiveTraits>::template LocalAssemblerImplementation>(
        mesh_elements, dof_table, local_assemblers, integration_order,
        is_axially_symmetric, process_data);
}

template void createLocalAssemblers<
    2, ConstitutiveStress_StrainTemperature::ConstitutiveTraits<2>>(
    std::vector<MeshLib::Element*> const& mesh_elements,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::unique_ptr<LocalAssemblerInterface<
        2, ConstitutiveStress_StrainTemperature::ConstitutiveTraits<2>>>>&
        local_assemblers,
    NumLib::IntegrationOrder const integration_order,
    bool const is_axially_symmetric,
    ThermoRichardsMechanicsProcessData<
        2, ConstitutiveStress_StrainTemperature::ConstitutiveTraits<2>>&
        process_data);

template void createLocalAssemblers<
    3, ConstitutiveStress_StrainTemperature::ConstitutiveTraits<3>>(
    std::vector<MeshLib::Element*> const& mesh_elements,
    NumLib::LocalToGlobalIndexMap const& dof_table,
    std::vector<std::unique_ptr<LocalAssemblerInterface<
        3, ConstitutiveStress_StrainTemperature::ConstitutiveTraits<3>>>>&
        local_assemblers,
    NumLib::IntegrationOrder const integration_order,
    bool const is_axially_symmetric,
    ThermoRichardsMechanicsProcessData<
        3, ConstitutiveStress_StrainTemperature::ConstitutiveTraits<3>>&
        process_data);

#ifdef OGS_USE_MFRONT
template void createLocalAssemblers<
    2,
    ConstitutiveStressSaturation_StrainPressureTemperature::ConstitutiveTraits<
        2>>(std::vector<MeshLib::Element*> const& mesh_elements,
            NumLib::LocalToGlobalIndexMap const& dof_table,
            std::vector<std::unique_ptr<LocalAssemblerInterface<
                2, ConstitutiveStressSaturation_StrainPressureTemperature::
                       ConstitutiveTraits<2>>>>& local_assemblers,
            NumLib::IntegrationOrder const integration_order,
            bool const is_axially_symmetric,
            ThermoRichardsMechanicsProcessData<
                2, ConstitutiveStressSaturation_StrainPressureTemperature::
                       ConstitutiveTraits<2>>& process_data);

template void createLocalAssemblers<
    3,
    ConstitutiveStressSaturation_StrainPressureTemperature::ConstitutiveTraits<
        3>>(std::vector<MeshLib::Element*> const& mesh_elements,
            NumLib::LocalToGlobalIndexMap const& dof_table,
            std::vector<std::unique_ptr<LocalAssemblerInterface<
                3, ConstitutiveStressSaturation_StrainPressureTemperature::
                       ConstitutiveTraits<3>>>>& local_assemblers,
            NumLib::IntegrationOrder const integration_order,
            bool const is_axially_symmetric,
            ThermoRichardsMechanicsProcessData<
                3, ConstitutiveStressSaturation_StrainPressureTemperature::
                       ConstitutiveTraits<3>>& process_data);
#endif
}  // namespace ProcessLib::ThermoRichardsMechanics
