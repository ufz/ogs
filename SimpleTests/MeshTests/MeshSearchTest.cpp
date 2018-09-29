/**
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 *
 * @file MeshSearchTest.cpp
 *
 *  Created on  Aug 30, 2012 by Thomas Fischer
 */

#include <tclap/CmdLine.h>
#include <logog/include/logog.hpp>

#include "BaseLib/LogogSimpleFormatter.h"
#include "BaseLib/MemWatch.h"
#include "BaseLib/RunTime.h"

#include "MeshLib/IO/Legacy/MeshIO.h"

#include "GeoLib/Grid.h"

#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"
#include "MeshLib/Mesh.h"

void testMeshGridAlgorithm(MeshLib::Mesh const*const mesh,
                std::vector<GeoLib::Point*>& pnts_for_search,
                std::vector<std::size_t> &idx_found_nodes, bool contiguous)
{
    // constructing Grid
    INFO ("[MeshGridAlgorithm] constructing mesh grid object ...");

    if (contiguous) {
        std::vector<MeshLib::Node> mesh_nodes;
        std::size_t n_nodes(mesh->getNodes().size());
        mesh_nodes.reserve(n_nodes);
        for (std::size_t k(0); k<n_nodes; k++) {
            mesh_nodes.emplace_back(*(mesh->getNodes()[k]));
        }
#ifndef WIN32
        BaseLib::MemWatch mem_watch;
        unsigned long mem_without_mesh (mem_watch.getVirtMemUsage());
#endif
        clock_t start_grid_construction = clock();
        GeoLib::Grid<MeshLib::Node> mesh_grid(mesh_nodes.begin(), mesh_nodes.end(), 511);
        clock_t end_grid_construction = clock();
#ifndef WIN32
        unsigned long mem_with_mesh (mem_watch.getVirtMemUsage());
#endif
        INFO("\tdone, construction time: %f seconds", (end_grid_construction-start_grid_construction)/(double)(CLOCKS_PER_SEC));
    #ifndef WIN32
        INFO ("[MeshGridAlgorithm] mem for mesh grid: %i MB", (mem_with_mesh - mem_without_mesh)/(1024*1024));
    #endif
        const std::size_t n_pnts_for_search(pnts_for_search.size());
        INFO ("[MeshGridAlgorithm] searching %d points ...", pnts_for_search.size());
        clock_t start = clock();
        for (std::size_t k(0); k<n_pnts_for_search; k++) {
            MeshLib::Node const* node(mesh_grid.getNearestPoint(*(pnts_for_search[k])));
            idx_found_nodes.push_back(node->getID());
        }
        clock_t stop = clock();
        INFO("\tdone, search time: %f seconds", (stop-start)/(double)(CLOCKS_PER_SEC));
    } else {
#ifndef WIN32
        BaseLib::MemWatch mem_watch;
        unsigned long mem_without_mesh (mem_watch.getVirtMemUsage());
#endif
        clock_t start_grid_construction = clock();
        GeoLib::Grid<MeshLib::Node> mesh_grid(mesh->getNodes().begin(), mesh->getNodes().end(), 511);
        clock_t end_grid_construction = clock();
#ifndef WIN32
        unsigned long mem_with_mesh (mem_watch.getVirtMemUsage());
#endif
        INFO("\tdone, construction time: %f seconds", (end_grid_construction-start_grid_construction)/(double)(CLOCKS_PER_SEC));
#ifndef WIN32
        INFO ("[MeshGridAlgorithm] mem for mesh grid: %i MB", (mem_with_mesh - mem_without_mesh)/(1024*1024));
#endif
        const std::size_t n_pnts_for_search(pnts_for_search.size());
        INFO ("[MeshGridAlgorithm] searching %d points ...", pnts_for_search.size());
        clock_t start = clock();
        for (std::size_t k(0); k<n_pnts_for_search; k++) {
            MeshLib::Node const* node(mesh_grid.getNearestPoint(*(pnts_for_search[k])));
            idx_found_nodes.push_back(node->getID());
        }
        clock_t stop = clock();
        INFO("\tdone, search time: %f seconds", (stop-start)/(double)(CLOCKS_PER_SEC));
    }
}

int main(int argc, char *argv[])
{
    LOGOG_INITIALIZE();
    auto* logog_cout(new logog::Cout);
    auto* custom_format(new BaseLib::LogogSimpleFormatter);
    logog_cout->SetFormatter(*custom_format);

    TCLAP::CmdLine cmd("Simple mesh search test", ' ', "0.1");

    // Define a value argument and add it to the command line.
    // A value arg defines a flag and a type of value that it expects,
    // such as "-m meshfile".
    TCLAP::ValueArg<std::string> mesh_arg("m","mesh","input mesh file",true,"test.msh","string");

    // Add the argument mesh_arg to the CmdLine object. The CmdLine object
    // uses this Arg to parse the command line.
    cmd.add( mesh_arg );

    TCLAP::ValueArg<unsigned> number_arg("n","number-of-test-points","the number of test points",true,10000,"positive number");
    cmd.add( number_arg );

    TCLAP::ValueArg<bool> contiguous_arg("c","use-contiguous-memory","use a contiguous memory for the test",false,true,"yes or no | 1 or 0");
    cmd.add( contiguous_arg );

    cmd.parse( argc, argv );

    std::string fname (mesh_arg.getValue());

    MeshLib::IO::Legacy::MeshIO mesh_io;
#ifndef WIN32
    BaseLib::MemWatch mem_watch;
    unsigned long mem_without_mesh (mem_watch.getVirtMemUsage());
#endif
    BaseLib::RunTime run_time;
    run_time.start();
    MeshLib::Mesh* mesh (mesh_io.loadMeshFromFile(fname));
#ifndef WIN32
    unsigned long mem_with_mesh (mem_watch.getVirtMemUsage());
    INFO ("mem for mesh: %i MB", (mem_with_mesh - mem_without_mesh)/(1024*1024));
#endif
    INFO ("time for reading: %f s", run_time.elapsed());

    // *** preparing test data
    std::vector<MeshLib::Node*> const& nodes(mesh->getNodes());
    std::vector<GeoLib::Point*> pnts_for_search;
    unsigned n(std::min(static_cast<unsigned>(nodes.size()), number_arg.getValue()));
    for (std::size_t k(0); k<n; k++) {
        pnts_for_search.push_back(new GeoLib::Point(*(nodes[k]), k));
    }

    std::vector<std::size_t> idx_found_nodes;
    testMeshGridAlgorithm(mesh, pnts_for_search, idx_found_nodes, contiguous_arg.getValue());

    for (std::size_t k(0); k<n; k++) {
        delete pnts_for_search[k];
    }

    delete mesh;
    delete custom_format;
    delete logog_cout;
    LOGOG_SHUTDOWN();
}
