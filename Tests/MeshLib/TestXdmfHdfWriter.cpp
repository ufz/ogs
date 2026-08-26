// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>
#include <hdf5.h>
#include <spdlog/fmt/fmt.h>

#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include <vector>

#include "BaseLib/StringTools.h"
#include "MeshLib/IO/XDMF/XdmfHdfWriter.h"
#include "MeshLib/Mesh.h"
#include "MeshLib/MeshEnums.h"
#include "MeshToolsLib/MeshGenerators/MeshGenerator.h"

// ---------------------------------------------------------------------------
// XdmfHdfWriter integration tests (file I/O)
//
// Mesh: 2×2 quads  →  9 nodes, 4 cells  (name "mesh")
// Properties:
//   "pressure"    – variable, float64, per-node
//   "MaterialIDs" – constant, int32,   per-cell
// ---------------------------------------------------------------------------

class XdmfHdfWriterTest : public ::testing::Test
{
public:
    XdmfHdfWriterTest()
        : tmp_dir(std::filesystem::temp_directory_path() /
                  BaseLib::randomString(12)),
          h5_path(tmp_dir / "out.h5"),
          static_h5_path(tmp_dir / "out_static.h5"),
          // 2×2 grid: (2+1)^2 = 9 nodes, 4 cells, mesh name = "mesh"
          mesh(MeshToolsLib::MeshGenerator::generateRegularQuadMesh(1.0, 2))
    {
        std::filesystem::create_directories(tmp_dir);

        auto& props = mesh->getProperties();

        auto* pressure = props.createNewPropertyVector<double>(
            "pressure", MeshLib::MeshItemType::Node, 1);
        pressure->resize(mesh->getNumberOfNodes(), 0.0);

        auto* static_data = props.createNewPropertyVector<double>(
            "example static attribute", MeshLib::MeshItemType::Node, 1);
        static_data->resize(mesh->getNumberOfNodes(), 0.0);

        auto* mat_ids = props.createNewPropertyVector<int>(
            "MaterialIDs", MeshLib::MeshItemType::Cell, 1);
        mat_ids->resize(mesh->getNumberOfElements(), 0);
    }

    ~XdmfHdfWriterTest() override { std::filesystem::remove_all(tmp_dir); }

protected:
    std::filesystem::path const tmp_dir;
    std::filesystem::path const h5_path;
    std::filesystem::path const static_h5_path;
    std::unique_ptr<MeshLib::Mesh> const mesh;
    static constexpr auto path_format_string = "/meshes/{}/{}";

    // Write mesh with "pressure" and "MaterialIDs" where pressure is a
    // variable attribute and "MaterialIDs" is a static attribute, output stem
    // = tmp_dir/"out".
    void createExampleFile(bool const store_static_data_separately) const
    {
        MeshLib::IO::XdmfHdfWriter writer(
            {std::cref(*mesh)}, tmp_dir / "out", 0, 0.0,
            {"MaterialIDs", "pressure"},
            {"example static attribute"} /* static_attribute_names */, false, 1,
            1024, store_static_data_separately);
    }

    // Return true when the named dataset exists under /meshes/<mesh>/ in
    // h5_path.
    bool datasetExists(std::filesystem::path const& h5_path,
                       std::string const& name) const
    {
        hid_t const f =
            H5Fopen(h5_path.string().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
        if (f < 0)
        {
            return false;
        }
        std::string const path =
            fmt::format(path_format_string, mesh->getName(), name);
        htri_t const exists = H5Lexists(f, path.c_str(), H5P_DEFAULT);
        H5Fclose(f);
        return exists > 0;
    }

    // Open a named dataset from out_static.h5 / out.h5. Returns
    // {file_id, dataset_id}; both negative on failure. Caller owns the handles
    // and must close them with H5Dclose / H5Fclose.
    std::pair<hid_t, hid_t> openDataset(std::string const& name,
                                        bool const open_static) const
    {
        hid_t const file = open_static
                               ? H5Fopen(static_h5_path.string().c_str(),
                                         H5F_ACC_RDONLY, H5P_DEFAULT)
                               : H5Fopen(h5_path.string().c_str(),
                                         H5F_ACC_RDONLY, H5P_DEFAULT);
        if (file < 0)
        {
            return {-1, -1};
        }
        std::string const path =
            fmt::format(path_format_string, mesh->getName(), name);
        return {file, H5Dopen2(file, path.c_str(), H5P_DEFAULT)};
    }
};

TEST_F(XdmfHdfWriterTest, SeparateStaticModeMakesTwoHDF5Files)
{
    {
        MeshLib::IO::XdmfHdfWriter writer(
            {std::cref(*mesh)}, tmp_dir / "out", 0, 0.0, {"pressure"},
            {} /* static_attribute_names */, false, 1, 1024,
            /*store_static_data_separately=*/true);
        writer.writeStep(1.0);
    }
    EXPECT_TRUE(std::filesystem::exists(h5_path));
    EXPECT_TRUE(std::filesystem::exists(static_h5_path));
}

TEST_F(XdmfHdfWriterTest, SingleFileModeDoesNotCreateStaticFile)
{
    {
        MeshLib::IO::XdmfHdfWriter writer(
            {std::cref(*mesh)}, tmp_dir / "out", 0, 0.0, {"pressure"},
            {} /* static_attribute_names */, false, 1, 1024,
            /*store_static_data_separately=*/false);
        writer.writeStep(1.0);
    }
    EXPECT_TRUE(std::filesystem::exists(h5_path));
    EXPECT_FALSE(std::filesystem::exists(static_h5_path));
}

TEST_F(XdmfHdfWriterTest, XdmfFileReferencesStaticFilenameForGeometry)
{
    createExampleFile(/*store_static_data_separately=*/true);

    auto const xdmf_path = tmp_dir / ("out_" + mesh->getName() + ".xdmf");
    ASSERT_TRUE(std::filesystem::exists(xdmf_path));

    std::ifstream f(xdmf_path);
    std::string const content((std::istreambuf_iterator<char>(f)),
                              std::istreambuf_iterator<char>());

    EXPECT_NE(content.find("out_static.h5"), std::string::npos)
        << "XDMF must reference the static HDF5 file for geometry/topology";
    EXPECT_NE(content.find("out.h5"), std::string::npos)
        << "XDMF must reference the dynamic HDF5 file for variable data";
}

TEST_F(XdmfHdfWriterTest, StaticFileGeometryHasTimeDimension)
{
    // Static file geometry: shape (1, n_nodes, 3) — rank 3, one time step
    // prepended
    createExampleFile(/*store_static_data_separately=*/true);

    auto const [file, ds] = openDataset("geometry", /*open_static=*/true);
    ASSERT_GE(file, 0);
    ASSERT_GE(ds, 0) << "geometry dataset must exist in static file";

    hid_t const sp = H5Dget_space(ds);
    EXPECT_EQ(H5Sget_simple_extent_ndims(sp), 3)
        << "static geometry must be rank-3 (time, n_nodes, 3)";

    hsize_t dims[3] = {0, 0, 0};
    H5Sget_simple_extent_dims(sp, dims, nullptr);
    EXPECT_EQ(dims[0], 1u) << "time dimension must be 1";
    EXPECT_EQ(dims[1], mesh->getNumberOfNodes());
    EXPECT_EQ(dims[2], 3u);

    H5Sclose(sp);
    H5Dclose(ds);
    H5Fclose(file);
}

TEST_F(XdmfHdfWriterTest, DynamicFileVariableAttributeHasTimeDimension)
{
    // Dynamic file pressure: shape (n_steps, n_nodes) — rank 2, time as axis 0
    constexpr unsigned n_steps = 3;  // constructor + 2 writeStep calls
    {
        MeshLib::IO::XdmfHdfWriter writer(
            {std::cref(*mesh)}, tmp_dir / "out", 0, 0.0, {"pressure"},
            {} /* static_attribute_names */, false, 1, 1024, true);
        for (unsigned i = 1; i < n_steps; ++i)
        {
            writer.writeStep(static_cast<double>(i));
        }
    }

    auto const [file, ds] = openDataset("pressure", false);
    ASSERT_GE(file, 0);
    ASSERT_GE(ds, 0) << "pressure dataset must exist in dynamic file";

    hid_t const sp = H5Dget_space(ds);
    EXPECT_EQ(H5Sget_simple_extent_ndims(sp), 2)
        << "dynamic scalar attribute must be rank-2 (n_steps x n_nodes)";

    std::vector<hsize_t> dims(2);
    H5Sget_simple_extent_dims(sp, dims.data(), nullptr);
    EXPECT_EQ(dims[0], static_cast<hsize_t>(n_steps));
    EXPECT_EQ(dims[1], mesh->getNumberOfNodes());

    H5Sclose(sp);
    H5Dclose(ds);
    H5Fclose(file);
}

TEST_F(XdmfHdfWriterTest, StaticFileIsNotModifiedAfterConstruction)
{
    // The static file is fully written in HdfWriter's constructor and closed.
    // Subsequent writeStep calls must not alter it.
    MeshLib::IO::XdmfHdfWriter writer(
        {std::cref(*mesh)}, tmp_dir / "out", 0, 0.0, {"pressure"},
        {} /* static_attribute_names */, false, 1, 1024, true);

    auto const mtime = std::filesystem::last_write_time(static_h5_path);

    writer.writeStep(1.0);
    writer.writeStep(2.0);

    // comparison via EXPECT_EQ results in compile error on Mac, error msg.:
    // 'to_chars' has been explicitly marked unavailable her
    // EXPECT_EQ(std::filesystem::last_write_time(static_h5_path), mtime)
    //     << "static file mtime must not change after writeStep calls";

    EXPECT_TRUE(std::filesystem::last_write_time(static_h5_path) == mtime)
        << "static file mtime must not change after writeStep calls";
}

TEST_F(XdmfHdfWriterTest, StaticFileAbsentFromDynamicFile)
{
    // In separate-file mode, constant attributes (geometry, topology,
    // MaterialIDs) must not appear in the dynamic file, and the variable
    // attribute (pressure) must not appear in the static file.
    createExampleFile(/*store_static_data_separately=*/true);

    EXPECT_FALSE(datasetExists(h5_path, "geometry"))
        << "geometry must NOT be in the dynamic file";
    EXPECT_FALSE(datasetExists(h5_path, "MaterialIDs"))
        << "MaterialIDs (constant) must NOT be in the dynamic file";
    EXPECT_FALSE(datasetExists(static_h5_path, "pressure"))
        << "pressure (variable) must NOT be in the static file";
    EXPECT_TRUE(datasetExists(static_h5_path, "geometry"))
        << "geometry must be in the static file";
}

TEST_F(XdmfHdfWriterTest, UserDeclaredStaticAttributeGoesToStaticFile)
{
    {
        MeshLib::IO::XdmfHdfWriter writer(
            {std::cref(*mesh)}, tmp_dir / "out",
            /*time_step=*/0,
            /*initial_time=*/0.0,
            /*output_variable_names=*/
            {"MaterialIDs", "pressure", "example static attribute"},
            /*static_attribute_names=*/{"example static attribute"},
            /*use_compression=*/false,
            /*n_files=*/1,
            /* chunk_size_bytes=*/1024,
            /*store_static_data_separately=*/true);
        writer.writeStep(1.0);
    }

    EXPECT_TRUE(datasetExists(static_h5_path, "example static attribute"))
        << "user-declared static attribute must be in the static file";
    EXPECT_FALSE(datasetExists(h5_path, "example static attribute"))
        << "user-declared static attribute must NOT be in the dynamic file";
}

TEST_F(XdmfHdfWriterTest, GlobalIdsAreAlwaysStatic)
{
    auto& props = mesh->getProperties();
    auto* g_nodes = props.createNewPropertyVector<std::size_t>(
        "global_node_ids", MeshLib::MeshItemType::Node, 1);
    g_nodes->resize(mesh->getNumberOfNodes(), 0);
    auto* g_elems = props.createNewPropertyVector<std::size_t>(
        "global_element_ids", MeshLib::MeshItemType::Cell, 1);
    g_elems->resize(mesh->getNumberOfElements(), 0);

    {
        MeshLib::IO::XdmfHdfWriter writer(
            {std::cref(*mesh)}, tmp_dir / "out", 0, 0.0,
            {"pressure", "global_node_ids", "global_element_ids"},
            {} /* static_attribute_names */, false, 1, 1024,
            /*store_static_data_separately=*/true);
        writer.writeStep(1.0);
    }

    EXPECT_TRUE(datasetExists(static_h5_path, "global_node_ids"));
    EXPECT_TRUE(datasetExists(static_h5_path, "global_element_ids"));
    EXPECT_FALSE(datasetExists(h5_path, "global_node_ids"));
    EXPECT_FALSE(datasetExists(h5_path, "global_element_ids"));
}

TEST_F(XdmfHdfWriterTest, MaterialIDsInStaticFileHasCorrectShape)
{
    // In separate-file mode, MaterialIDs in the static hdf5 file must be
    // rank-2 with shape (1, n_elements) — one time step prepended.
    createExampleFile(/*store_static_data_separately=*/true);

    auto const [file, ds] = openDataset("MaterialIDs", /*open_static=*/true);
    ASSERT_GE(file, 0);
    ASSERT_GE(ds, 0) << "MaterialIDs must exist in static file";

    hid_t const sp = H5Dget_space(ds);
    EXPECT_EQ(H5Sget_simple_extent_ndims(sp), 2)
        << "MaterialIDs must be rank-2 (time, n_cells)";

    hsize_t dims[2] = {0, 0};
    H5Sget_simple_extent_dims(sp, dims, nullptr);
    EXPECT_EQ(dims[0], 1u) << "time dimension must be 1";
    EXPECT_EQ(dims[1], mesh->getNumberOfElements());

    H5Sclose(sp);
    H5Dclose(ds);
    H5Fclose(file);
}

TEST_F(XdmfHdfWriterTest, MaterialIDsInStaticFileHasCorrectValues)
{
    auto* mat_ids = mesh->getProperties().getPropertyVector<int>("MaterialIDs");
    ASSERT_NE(mat_ids, nullptr);
    std::vector<int> const expected = {10, 20, 30, 40};
    for (std::size_t i = 0; i < expected.size(); ++i)
    {
        (*mat_ids)[i] = expected[i];
    }

    createExampleFile(/*store_static_data_separately=*/true);

    auto const [file, ds] = openDataset("MaterialIDs", /*open_static=*/true);
    ASSERT_GE(file, 0);
    ASSERT_GE(ds, 0);

    std::vector<int> actual(expected.size());
    herr_t const status = H5Dread(ds, H5T_NATIVE_INT, H5S_ALL, H5S_ALL,
                                  H5P_DEFAULT, actual.data());
    ASSERT_GE(status, 0);
    EXPECT_EQ(actual, expected);

    H5Dclose(ds);
    H5Fclose(file);
}
