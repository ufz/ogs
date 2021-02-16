/**
 * \file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "RandomFieldMeshElementParameter.h"

#include "BaseLib/ConfigTree.h"
#include "MeshLib/Mesh.h"
#include <random>
#include <functional>

namespace ParameterLib
{
std::unique_ptr<ParameterBase> createRandomFieldMeshElementParameter(
    std::string const& name, BaseLib::ConfigTree const& config,
    MeshLib::Mesh& mesh)
{
    //! \ogs_file_param{prj__parameters__parameter__type}
    config.checkConfigParameter("type", "RandomFieldMeshElement");
    auto const field_name =
        //! \ogs_file_param{prj__parameters__parameter__RandomFieldMeshElementParameter__field_name}
        config.getConfigParameter<std::string>("field_name");
    auto const range =
        //! \ogs_file_param{prj__parameters__parameter__RandomFieldMeshElementParameter__range}
        config.getConfigParameter<std::vector<double>>("range");
    if (range.size() != 2)
    {
        OGS_FATAL(
            "The range needs to have two components, but {:d} were given.",
            range.size());
    }
    auto const seed =
        //! \ogs_file_param{prj__parameters__parameter__RandomFieldMeshElementParameter__seed}
        config.getConfigParameter<int>("seed");
    DBUG("Generating field {:s} with range {:g} to {:g} and seed {:d}.",
         field_name, range[0], range[1], seed);

    std::vector<double> values(mesh.getElements().size());

    std::mt19937 generator(seed);
    std::uniform_real_distribution<> distr(range[0], range[1]);
    auto gen = [&distr, &generator]() { return distr(generator); };
    generate(begin(values), end(values), gen);

    MeshLib::addPropertyToMesh(mesh, field_name, MeshLib::MeshItemType::Cell, 1,
                               values);

    auto const& property =
        mesh.getProperties().getPropertyVector<double>(field_name);

    if (property->getMeshItemType() != MeshLib::MeshItemType::Cell)
    {
        OGS_FATAL("The mesh property `{:s}' is not an element property.",
                  field_name);
    }

    return std::make_unique<RandomFieldMeshElementParameter<double>>(name, mesh,
                                                                     *property);
}

}  // namespace ParameterLib
