/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_PROCESS_VARIABLE_H_
#define PROCESS_LIB_PROCESS_VARIABLE_H_

#include <boost/property_tree/ptree.hpp>

#include "GeoLib/GEOObjects.h"
#include "MeshLib/Mesh.h"

#include "MeshGeoToolsLib/MeshNodeSearcher.h"
#include "MeshGeoToolsLib/BoundaryElementsSearcher.h"

#include "NeumannBcConfig.h"
#include "NeumannBC.h"

namespace ProcessLib
{
    class InitialCondition;
    class UniformDirichletBoundaryCondition;
}

namespace ProcessLib
{

/// A named process variable. Its properties includes the mesh, and the initial
/// and boundary conditions.
class ProcessVariable
{
    using ConfigTree = boost::property_tree::ptree;
public:
    ProcessVariable(ConfigTree const& config, MeshLib::Mesh const& mesh,
            GeoLib::GEOObjects const& geometries);

    ~ProcessVariable();

    std::string const& getName() const;

    /// Returns a mesh on which the process variable is defined.
    MeshLib::Mesh const& getMesh() const;

    void initializeDirichletBCs(MeshGeoToolsLib::MeshNodeSearcher& searcher,
            std::vector<std::size_t>& global_ids, std::vector<double>& values);

    template <typename OutputIterator, typename GlobalSetup, typename ...Args>
    void createNeumannBCs(OutputIterator bcs,
        MeshGeoToolsLib::BoundaryElementsSearcher& searcher,
        GlobalSetup const& global_setup,
        Args&&... args)
    {
        for (NeumannBcConfig* config : _neumann_bc_configs)
        {
            config->initialize(searcher);
            bcs = new NeumannBC<GlobalSetup>(config, std::forward<Args>(args)...);
        }
    }

private:
    std::string const _name;
    MeshLib::Mesh const& _mesh;
    InitialCondition* _initial_condition;
    std::vector<UniformDirichletBoundaryCondition*> _dirichlet_bcs;
    std::vector<NeumannBcConfig*> _neumann_bc_configs;
};

}   // namespace ProcessLib

#endif  // PROCESS_LIB_PROCESS_VARIABLE_H_
