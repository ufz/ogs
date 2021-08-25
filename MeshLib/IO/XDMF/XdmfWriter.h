/**
 * \file
 * \author Tobias Meisel
 * \date   2021-06-30
 * \brief  Collects and holds all metadata for writing XDMF file
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <functional>
#include <string>

namespace MeshLib::IO
{
class XdmfWriter final
{
public:
    /**
     * \brief Writes xdmf string into file on class destruction
     * @param xdmf_filename absolute or relative filepath to the xdmf file
     * @param xdmf_writer_fn function that generates xdmf string
     */
    XdmfWriter(std::string xdmf_filename,
               std::function<std::string(std::vector<double>)>
                   xdmf_writer_fn);
    XdmfWriter(XdmfWriter&&) = default;
    XdmfWriter& operator=(XdmfWriter&&) = default;
    XdmfWriter(XdmfWriter const&) = delete;
    XdmfWriter& operator=(XdmfWriter const&) = delete;
    ~XdmfWriter();

    /**
     * \brief Adds data for lazy (xdmf) writing algorithm
     * @param time_step time value of the current time_step
     */
    void addTimeStep(double const& time_step);

private:
    std::string filename;
    std::vector<double> times;
    std::function<std::string(std::vector<double>)> xdmf_writer;
};
}  // namespace MeshLib::IO
