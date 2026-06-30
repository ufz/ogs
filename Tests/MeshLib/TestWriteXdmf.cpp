// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "MeshLib/IO/XDMF/MeshPropertyDataType.h"
#include "MeshLib/IO/XDMF/XdmfData.h"
#include "MeshLib/IO/XDMF/writeXdmf.h"
#include "MeshLib/MeshEnums.h"

// ---------------------------------------------------------------------------
// write_xdmf string-generation tests  (no file I/O)
//
// All datasets (static and dynamic) carry a leading time dimension.
//
// Geometry: 4 nodes × 3 coords  → starts=[0,0], strides=[1,1]
//   static  hyperslab: |0 0 0:1 1 1:1 4 3:1 4 3
//   dynamic hyperslab: |<step> 0 0:1 1 1:1 4 3:<max_step> 4 3
//
// Scalar attribute: 4 nodes × 1  → starts=[0], strides=[1]
//   static  hyperslab: |0 0:1 1:1 4:1 4
//   dynamic hyperslab: |<step> 0:1 1:1 4:<max_step> 4
// ---------------------------------------------------------------------------

namespace
{
// Minimal XdmfData helpers for 4-node single-quad mesh (serial, n_files=1).
MeshLib::IO::XdmfData makeGeometry()
{
    return {
        4, 3,           MeshPropertyDataType::float64, "geometry", std::nullopt,
        1, std::nullopt};
}

MeshLib::IO::XdmfData makeTopology()
{
    return {4,
            1,
            MeshPropertyDataType::uint64,
            "topology",
            std::nullopt,
            1,
            ParentDataType::QUADRILATERAL};
}

MeshLib::IO::XdmfData makeScalarNodeAttribute(std::string const& name)
{
    return {4,
            1,
            MeshPropertyDataType::float64,
            name,
            MeshLib::MeshItemType::Node,
            1,
            std::nullopt};
}
}  // namespace

TEST(WriteXdmf, StaticDataItemReferencesStaticFile)
{
    auto const fn = MeshLib::IO::write_xdmf(
        makeGeometry(), makeTopology(), {}, {makeScalarNodeAttribute("p")},
        "out_static.h5", "out.h5", "1.0", "m");

    std::string const xdmf = fn({0.0});

    EXPECT_NE(xdmf.find("out_static.h5:/meshes/m/geometry"), std::string::npos)
        << "geometry must reference static file";
    EXPECT_NE(xdmf.find("out_static.h5:/meshes/m/topology"), std::string::npos)
        << "topology must reference static file";
    EXPECT_NE(xdmf.find("out.h5:/meshes/m/p"), std::string::npos)
        << "variable attribute must reference dynamic file";

    EXPECT_EQ(xdmf.find("out.h5:/meshes/m/geometry"), std::string::npos)
        << "geometry must NOT reference dynamic file";
    EXPECT_EQ(xdmf.find("out_static.h5:/meshes/m/p"), std::string::npos)
        << "variable attribute must NOT reference static file";
}

TEST(WriteXdmf, BackwardCompatModeAllDataInSingleFile)
{
    // When store_static_data_separately=false the caller passes the dynamic
    // filename for both static_h5filename and dynamic_h5filename.
    auto const fn = MeshLib::IO::write_xdmf(makeGeometry(), makeTopology(), {},
                                            {makeScalarNodeAttribute("p")},
                                            "out.h5", "out.h5", "1.0", "m");

    std::string const xdmf = fn({0.0});

    EXPECT_NE(xdmf.find("out.h5:/meshes/m/geometry"), std::string::npos);
    EXPECT_NE(xdmf.find("out.h5:/meshes/m/p"), std::string::npos);
    EXPECT_EQ(xdmf.find("_static.h5"), std::string::npos)
        << "static filename must not appear in backward-compat XDMF";
}

TEST(WriteXdmf, StaticGeometryHyperslabHasLeadingTimeStep)
{
    // Static geometry (2D, starts=[0,0]) is now written with a time dimension.
    // Its DataItem always uses time_step=0, max_step=1:
    //   hyperslab = |0 0 0:1 1 1:1 4 3:1 4 3
    auto const fn =
        MeshLib::IO::write_xdmf(makeGeometry(), makeTopology(), {}, {},
                                "static.h5", "dynamic.h5", "1.0", "m");

    std::string const xdmf = fn({0.0, 1.0});

    EXPECT_NE(xdmf.find("static.h5:/meshes/m/geometry|0 0 0:1"),
              std::string::npos)
        << "static geometry hyperslab must carry leading time step 0";
    EXPECT_EQ(xdmf.find("static.h5:/meshes/m/geometry|0 0:1"),
              std::string::npos)
        << "static geometry must not use the no-time-dim hyperslab format";
}

TEST(WriteXdmf, DynamicHyperslabCarriesTimeStepIndex)
{
    // Scalar attribute (1D, starts=[0]): hyperslab = |<step> 0:1 1:1 4:4
    auto const fn = MeshLib::IO::write_xdmf(makeGeometry(), makeTopology(), {},
                                            {makeScalarNodeAttribute("p")},
                                            "s.h5", "d.h5", "1.0", "m");

    std::string const xdmf = fn({0.0, 1.0});

    EXPECT_NE(xdmf.find("d.h5:/meshes/m/p|0 0:1"), std::string::npos)
        << "step 0: time index 0 before spatial starts";
    EXPECT_NE(xdmf.find("d.h5:/meshes/m/p|1 0:1"), std::string::npos)
        << "step 1: time index 1 before spatial starts";
}

TEST(WriteXdmf, ConstantAttributeUsesStaticFile)
{
    // Constant attributes (non-variable) should also go to the static file.
    MeshLib::IO::XdmfData const mat_ids{4,
                                        1,
                                        MeshPropertyDataType::int32,
                                        "MaterialIDs",
                                        MeshLib::MeshItemType::Cell,
                                        1,
                                        std::nullopt};

    auto const fn =
        MeshLib::IO::write_xdmf(makeGeometry(), makeTopology(), {mat_ids}, {},
                                "st.h5", "dyn.h5", "1.0", "m");

    std::string const xdmf = fn({0.0});

    EXPECT_NE(xdmf.find("st.h5:/meshes/m/MaterialIDs"), std::string::npos)
        << "constant attribute must reference static file";
    EXPECT_EQ(xdmf.find("dyn.h5:/meshes/m/MaterialIDs"), std::string::npos)
        << "constant attribute must NOT reference dynamic file";
}
