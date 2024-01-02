/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on September 29, 2022, 3:22 PM
 */

#include "IntegrationPointDataTools.h"

#include "MeshLib/Elements/Elements.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Properties.h"
#include "MeshLib/PropertyVector.h"
#include "MeshLib/Utils/IntegrationPointWriter.h"
#include "NumLib/Fem/Integration/GaussLegendreIntegrationPolicy.h"

namespace MeshToolsLib
{
template <typename ElementType>
int getNumberOfElementIntegrationPointsGeneral(
    MeshLib::IntegrationPointMetaData const& ip_meta_data)
{
    using IntegrationPolicy =
        NumLib::GaussLegendreIntegrationPolicy<ElementType>;
    using NonGenericIntegrationMethod =
        typename IntegrationPolicy::IntegrationMethod;
    NonGenericIntegrationMethod int_met{
        (unsigned int)ip_meta_data.integration_order};
    return int_met.getNumberOfPoints();
}

int getNumberOfElementIntegrationPoints(
    MeshLib::IntegrationPointMetaData const& ip_meta_data,
    MeshLib::Element const& e)
{
    switch (e.getCellType())
    {
        case MeshLib::CellType::LINE2:
            return getNumberOfElementIntegrationPointsGeneral<MeshLib::Line>(
                ip_meta_data);
        case MeshLib::CellType::LINE3:
            return getNumberOfElementIntegrationPointsGeneral<MeshLib::Line3>(
                ip_meta_data);
        case MeshLib::CellType::TRI3:
            return getNumberOfElementIntegrationPointsGeneral<MeshLib::Tri>(
                ip_meta_data);
        case MeshLib::CellType::TRI6:
            return getNumberOfElementIntegrationPointsGeneral<MeshLib::Tri6>(
                ip_meta_data);
        case MeshLib::CellType::QUAD4:
            return getNumberOfElementIntegrationPointsGeneral<MeshLib::Quad>(
                ip_meta_data);
        case MeshLib::CellType::QUAD8:
            return getNumberOfElementIntegrationPointsGeneral<MeshLib::Quad8>(
                ip_meta_data);
        case MeshLib::CellType::QUAD9:
            return getNumberOfElementIntegrationPointsGeneral<MeshLib::Quad9>(
                ip_meta_data);
        case MeshLib::CellType::TET4:
            return getNumberOfElementIntegrationPointsGeneral<MeshLib::Tet>(
                ip_meta_data);
        case MeshLib::CellType::TET10:
            return getNumberOfElementIntegrationPointsGeneral<MeshLib::Tet10>(
                ip_meta_data);
        case MeshLib::CellType::HEX8:
            return getNumberOfElementIntegrationPointsGeneral<MeshLib::Hex>(
                ip_meta_data);
        case MeshLib::CellType::HEX20:
            return getNumberOfElementIntegrationPointsGeneral<MeshLib::Hex20>(
                ip_meta_data);
        case MeshLib::CellType::PRISM6:
            return getNumberOfElementIntegrationPointsGeneral<MeshLib::Prism>(
                ip_meta_data);
        case MeshLib::CellType::PRISM15:
            return getNumberOfElementIntegrationPointsGeneral<MeshLib::Prism15>(
                ip_meta_data);
        case MeshLib::CellType::PYRAMID5:
            return getNumberOfElementIntegrationPointsGeneral<MeshLib::Pyramid>(
                ip_meta_data);
        case MeshLib::CellType::PYRAMID13:
            return getNumberOfElementIntegrationPointsGeneral<
                MeshLib::Pyramid13>(ip_meta_data);
        case MeshLib::CellType::PRISM18:
        case MeshLib::CellType::HEX27:
            OGS_FATAL("Mesh element type {:s} is not supported",
                      MeshLib::CellType2String(e.getCellType()));
        case MeshLib::CellType::INVALID:
        default:
            OGS_FATAL("Invalid Element Type");
    }
    return 0;
}

std::vector<std::size_t> getIntegrationPointDataOffsetsOfMeshElements(
    std::vector<MeshLib::Element*> const& mesh_elements,
    MeshLib::PropertyVectorBase const& pv,
    MeshLib::Properties const& properties)
{
    // For special field data such as OGS_VERSION, IntegrationPointMetaData,
    // etc., which are not "real" integration points:
    if (pv.getPropertyName().find("_ip") == std::string::npos)
    {
        return {};
    }

    auto const n_components = pv.getNumberOfGlobalComponents();

    std::vector<std::size_t> element_ip_data_offsets(mesh_elements.size() + 1);
    std::size_t counter = 0;
    auto const ip_meta_data =
        MeshLib::getIntegrationPointMetaData(properties, pv.getPropertyName());
    for (std::size_t i = 0; i < mesh_elements.size(); i++)
    {
        auto const* const element = mesh_elements[i];

        // Assuming that the order of elements in mesh_elements is not touched.
        element_ip_data_offsets[i] = counter;
        counter += getNumberOfElementIntegrationPoints(ip_meta_data, *element) *
                   n_components;
    }
    element_ip_data_offsets[mesh_elements.size()] = counter;

    return element_ip_data_offsets;
}

}  // namespace MeshToolsLib
