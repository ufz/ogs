// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "Applications/FileIO/GocadIO/GocadAsciiReader.h"
#include "BaseLib/StringTools.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"
#include "MeshLib/Properties.h"
#include "MeshLib/PropertyVector.h"

/// A GoCAD TSurf data set may consist of several TFACE sections that are all
/// parsed into a single mesh. Writes such a file with three TFACE sections,
/// each carrying one node property column. The last section holds two elements
/// so that a section contributing more than one element is covered.
struct GocadAsciiReaderSurfacesTest : public ::testing::Test
{
    GocadAsciiReaderSurfacesTest()
    {
        std::ofstream out(file_name);
        out << "GOCAD TSurf 1\n";
        out << "HEADER {\n";
        out << "name:Surfaces\n";
        out << "}\n";
        out << "PROPERTIES Temperature\n";
        out << "ESIZES 1\n";
        out << "TFACE\n";
        out << "PVRTX 1 0 0 0 10\n";
        out << "PVRTX 2 1 0 0 20\n";
        out << "PVRTX 3 0 1 0 30\n";
        out << "TRGL 1 2 3\n";
        out << "TFACE\n";
        out << "PVRTX 4 2 0 0 40\n";
        out << "PVRTX 5 3 0 0 50\n";
        out << "PVRTX 6 2 1 0 60\n";
        out << "TRGL 4 5 6\n";
        out << "TFACE\n";
        out << "PVRTX 7 4 0 0 70\n";
        out << "PVRTX 8 5 0 0 80\n";
        out << "PVRTX 9 4 1 0 90\n";
        out << "PVRTX 10 5 1 0 100\n";
        out << "TRGL 7 8 9\n";
        out << "TRGL 8 10 9\n";
        out << "END\n";
        out.close();
    }

    ~GocadAsciiReaderSurfacesTest() override { std::remove(file_name.c_str()); }

    std::string const file_name = (std::filesystem::temp_directory_path() /=
                                   BaseLib::randomString(32) + ".ts")
                                      .string();
};

/// A GoCAD PLine data set may consist of several ILINE sections that are all
/// parsed into a single mesh. Writes such a file with two ILINE sections, the
/// first one holding two segments, the second one a single segment.
struct GocadAsciiReaderLinesTest : public ::testing::Test
{
    GocadAsciiReaderLinesTest()
    {
        std::ofstream out(file_name);
        out << "GOCAD PLine 1\n";
        out << "HEADER {\n";
        out << "name:Lines\n";
        out << "}\n";
        out << "PROPERTIES Temperature\n";
        out << "ESIZES 1\n";
        out << "ILINE\n";
        out << "PVRTX 1 0 0 0 10\n";
        out << "PVRTX 2 1 0 0 20\n";
        out << "PVRTX 3 2 0 0 30\n";
        out << "SEG 1 2\n";
        out << "SEG 2 3\n";
        out << "ILINE\n";
        out << "PVRTX 4 3 0 0 40\n";
        out << "PVRTX 5 4 0 0 50\n";
        out << "SEG 4 5\n";
        out << "END\n";
        out.close();
    }

    ~GocadAsciiReaderLinesTest() override { std::remove(file_name.c_str()); }

    std::string const file_name = (std::filesystem::temp_directory_path() /=
                                   BaseLib::randomString(32) + ".pl")
                                      .string();
};

// The segments of all ILINE sections of one data set are parsed into a single
// mesh, and the node properties and material IDs are handled the same way as
// for surfaces.
TEST_F(GocadAsciiReaderLinesTest, ParsesSegmentsOfAllLines)
{
    std::vector<std::unique_ptr<MeshLib::Mesh>> meshes;
    ASSERT_TRUE(FileIO::Gocad::GocadAsciiReader::readFile(file_name, meshes));
    ASSERT_EQ(1u, meshes.size());

    MeshLib::Mesh const& mesh = *meshes[0];
    ASSERT_EQ(5u, mesh.getNumberOfNodes());
    ASSERT_EQ(3u, mesh.getNumberOfElements());
    for (auto const* const element : mesh.getElements())
    {
        EXPECT_EQ(MeshLib::MeshElemType::LINE, element->getGeomType());
    }

    auto const* const temperature =
        mesh.getProperties().getPropertyVector<double>(
            "Temperature", MeshLib::MeshItemType::Node, 1);
    ASSERT_NE(nullptr, temperature);
    EXPECT_EQ((std::vector<double>{10, 20, 30, 40, 50}),
              std::vector<double>(temperature->begin(), temperature->end()));

    auto const* const mat_ids = mesh.getProperties().getPropertyVector<int>(
        "MaterialIDs", MeshLib::MeshItemType::Cell, 1);
    ASSERT_NE(nullptr, mat_ids);
    EXPECT_EQ((std::vector<int>{0, 0, 1}),
              std::vector<int>(mat_ids->begin(), mat_ids->end()));
}

/// A data set with invalid contents, named after the defect it carries.
struct InvalidGocadFile
{
    std::string name;
    std::string content;
};

/// Writes a data set with invalid contents to a temporary file that is removed
/// again when the test ends.
struct GocadAsciiReaderInvalidFileTest
    : public ::testing::TestWithParam<InvalidGocadFile>
{
    ~GocadAsciiReaderInvalidFileTest() override
    {
        std::remove(file_name.c_str());
    }

    void writeFile(std::string const& content) const
    {
        std::ofstream out(file_name);
        out << content;
        out.close();
    }

    std::string const file_name = (std::filesystem::temp_directory_path() /=
                                   BaseLib::randomString(32) + ".gcd")
                                      .string();
};

// Reading a malformed data set must abort the parsing and leave no mesh
// behind, instead of silently constructing a mesh from incomplete input.
TEST_P(GocadAsciiReaderInvalidFileTest, RejectsMalformedFile)
{
    writeFile(GetParam().content);

    std::vector<std::unique_ptr<MeshLib::Mesh>> meshes;
    ASSERT_FALSE(FileIO::Gocad::GocadAsciiReader::readFile(file_name, meshes));
    EXPECT_TRUE(meshes.empty());
}

INSTANTIATE_TEST_SUITE_P(
    GocadAsciiReader,
    GocadAsciiReaderInvalidFileTest,
    ::testing::Values(
        // An element referencing an undefined node ID would create an element
        // with a node that does not exist.
        InvalidGocadFile{"SurfaceWithUnknownNodeId",
                         "GOCAD TSurf 1\n"
                         "HEADER {\n"
                         "name:Surface\n"
                         "}\n"
                         "TFACE\n"
                         "VRTX 1 0 0 0\n"
                         "VRTX 2 1 0 0\n"
                         "VRTX 3 0 1 0\n"
                         "TRGL 1 2 99\n"
                         "END\n"},
        InvalidGocadFile{"LineWithUnknownNodeId",
                         "GOCAD PLine 1\n"
                         "HEADER {\n"
                         "name:Line\n"
                         "}\n"
                         "ILINE\n"
                         "VRTX 1 0 0 0\n"
                         "VRTX 2 1 0 0\n"
                         "SEG 1 99\n"
                         "END\n"},
        // A TRGL line must carry one node ID per triangle corner. A missing ID
        // would create an element referencing node 0.
        InvalidGocadFile{"SurfaceWithTruncatedTriangle",
                         "GOCAD TSurf 1\n"
                         "HEADER {\n"
                         "name:Surface\n"
                         "}\n"
                         "TFACE\n"
                         "VRTX 1 0 0 0\n"
                         "VRTX 2 1 0 0\n"
                         "VRTX 3 0 1 0\n"
                         "TRGL 1 2\n"
                         "END\n"},
        // A PVRTX line must carry one value per declared property. A missing
        // value would silently append a zero, because a failed extraction
        // value-initialises its argument.
        InvalidGocadFile{"PVRTXWithMissingPropertyValue",
                         "GOCAD TSurf 1\n"
                         "HEADER {\n"
                         "name:Surface\n"
                         "}\n"
                         "PROPERTIES Temperature\n"
                         "ESIZES 1\n"
                         "TFACE\n"
                         "PVRTX 1 0 0 0 10\n"
                         "PVRTX 2 1 0 0\n"
                         "PVRTX 3 0 1 0 30\n"
                         "TRGL 1 2 3\n"
                         "END\n"}),
    [](::testing::TestParamInfo<InvalidGocadFile> const& info)
    { return info.param.name; });

// The node properties of all TFACE sections of one data set are collected in
// the same property vector. Every section must append its values, so that the
// vector holds one value per node of the whole mesh and the values parsed for
// earlier sections survive.
TEST_F(GocadAsciiReaderSurfacesTest, AppendsNodePropertiesAcrossSurfaces)
{
    std::vector<std::unique_ptr<MeshLib::Mesh>> meshes;
    ASSERT_TRUE(FileIO::Gocad::GocadAsciiReader::readFile(file_name, meshes));
    ASSERT_EQ(1u, meshes.size());

    MeshLib::Mesh const& mesh = *meshes[0];
    ASSERT_EQ(10u, mesh.getNumberOfNodes());
    ASSERT_EQ(4u, mesh.getNumberOfElements());

    auto const* const temperature =
        mesh.getProperties().getPropertyVector<double>(
            "Temperature", MeshLib::MeshItemType::Node, 1);
    ASSERT_NE(nullptr, temperature);
    EXPECT_EQ((std::vector<double>{10, 20, 30, 40, 50, 60, 70, 80, 90, 100}),
              std::vector<double>(temperature->begin(), temperature->end()));
}

// Each TFACE section appends its elements to the shared MaterialIDs vector and
// must receive a fresh material ID one larger than the current maximum, while
// the IDs already assigned to earlier sections stay untouched. The former
// in-place `(*max_element)++` reused the old maximum and mutated an existing
// element.
TEST_F(GocadAsciiReaderSurfacesTest, AssignsFreshMaterialIdPerSurface)
{
    std::vector<std::unique_ptr<MeshLib::Mesh>> meshes;
    ASSERT_TRUE(FileIO::Gocad::GocadAsciiReader::readFile(file_name, meshes));
    ASSERT_EQ(1u, meshes.size());

    auto const* const mat_ids =
        meshes[0]->getProperties().getPropertyVector<int>(
            "MaterialIDs", MeshLib::MeshItemType::Cell, 1);
    ASSERT_NE(nullptr, mat_ids);
    EXPECT_EQ((std::vector<int>{0, 1, 2, 2}),
              std::vector<int>(mat_ids->begin(), mat_ids->end()));
}
