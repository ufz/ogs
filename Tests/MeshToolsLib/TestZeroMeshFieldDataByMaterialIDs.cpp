/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on July 10, 2023, 1:32 PM
 */
#include <gtest/gtest.h>

#include <algorithm>
#include <string>
#include <vector>

#include "InfoLib/TestInfo.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/IO/VtkIO/VtuInterface.h"
#include "MeshLib/IO/readMeshFromFile.h"
#include "MeshLib/Mesh.h"
#include "MeshToolsLib/ZeroMeshFieldDataByMaterialIDs.h"

// Note: the two meshes are the same except the values of double type field
// data.
bool compareMeshDoubleTypeIpData(MeshLib::Mesh const& mesh_a,
                                 MeshLib::Mesh const& mesh_b)
{
    MeshLib::Properties const& properties_a = mesh_a.getProperties();
    MeshLib::Properties const& properties_b = mesh_b.getProperties();

    for (auto const& [name, property] : properties_a)
    {
        auto const item_type = property->getMeshItemType();

        if (item_type != MeshLib::MeshItemType::IntegrationPoint)
        {
            continue;
        }

        // For special field data such as OGS_VERSION,
        // IntegrationPointMetaData,
        // etc., which are not "real" integration points:
        if (property->getPropertyName().find("_ip") == std::string::npos)
        {
            continue;
        }

        if (properties_a.template hasPropertyVector<double>(
                name, MeshLib::MeshItemType::IntegrationPoint))
        {
            auto const& pv_a =
                *properties_a.template getPropertyVector<double>(name);
            auto const& pv_b =
                *properties_b.template getPropertyVector<double>(name);
            if (!std::equal(pv_a.begin(), pv_a.end(), pv_b.begin()))
            {
                return false;
            }
        }
    }

    return true;
}

void compareMeshes(std::string const& mesh_file_name,
                   std::string const& ref_mesh_file_name,
                   std::vector<int> const& material_ids)
{
    auto mesh = std::unique_ptr<MeshLib::Mesh>(
        MeshLib::IO::readMeshFromFile(mesh_file_name));
    MeshToolsLib::zeroMeshFieldDataByMaterialIDs(*mesh, material_ids);

    auto const ref_mesh = std::unique_ptr<MeshLib::Mesh>(
        MeshLib::IO::readMeshFromFile(ref_mesh_file_name));

    ASSERT_TRUE(compareMeshDoubleTypeIpData(*mesh, *ref_mesh));
}

TEST(MeshLib, ZeroMeshFieldDataByMaterialIDs)
{
    // 2D
    {
        std::string const mesh_file_name =
            TestInfoLib::TestInfo::data_path +
            "/ThermoMechanics/CreepBGRa/CreepAfterExcavation/"
            "CreepAfterExcavation_ts_61_t_4320000.000000.vtu";
        std::string ref_mesh_file_name =
            TestInfoLib::TestInfo::data_path +
            "/MeshLib/2D_mesh_with_reset_sigma_ip.vtu";

        std::vector<int> const material_ids = {1};
        compareMeshes(mesh_file_name, ref_mesh_file_name, material_ids);
    }

    // 3D
    {
        std::string const mesh_file_name =
            TestInfoLib::TestInfo::data_path +
            "/NodePartitionedMesh/WithIntegrationPointStress/MixedElements/"
            "mesh_with_3D_different_elements_sigma_ip.vtu";

        // reset the elements with materialIDs 0, 2:
        {
            std::vector<int> const material_ids = {0, 2};
            std::string ref_mesh_file_name =
                TestInfoLib::TestInfo::data_path +
                "/MeshLib/3D_mesh_with_reset_sigma_ip_1.vtu";

            compareMeshes(mesh_file_name, ref_mesh_file_name, material_ids);
        }

        // reset the elements with materialIDs 1, 3:
        {
            std::vector<int> const material_ids = {1, 3};
            std::string ref_mesh_file_name =
                TestInfoLib::TestInfo::data_path +
                "/MeshLib/3D_mesh_with_reset_sigma_ip_2.vtu";

            compareMeshes(mesh_file_name, ref_mesh_file_name, material_ids);
        }
    }
}
