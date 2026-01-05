// SPDX-FileCopyrightText: Copyright (c) OpenGeoSys Community (opengeosys.org)
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <functional>
#include <string>
#include <vector>

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
