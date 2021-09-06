#include "XdmfWriter.h"

#include <fstream>

#include "BaseLib/Logging.h"
#include "BaseLib/RunTime.h"
#include "writeXdmf.h"

namespace MeshLib::IO
{
XdmfWriter::XdmfWriter(std::string xdmf_filename,
                       std::function<std::string(std::vector<double>)>
                           xdmf_writer_fn)
    : filename(std::move(xdmf_filename)), xdmf_writer(xdmf_writer_fn)
{
}

XdmfWriter::~XdmfWriter()
{
    BaseLib::RunTime time_output;
    time_output.start();
    std::ofstream fout;
    fout.open(filename);
    fout << xdmf_writer(times);

    INFO("[time] Output of XDMF to {:s} took {:g} s.",
         filename,
         time_output.elapsed());
}

void XdmfWriter::addTimeStep(double const& time_step)
{
    times.push_back(time_step);
}

}  // namespace MeshLib::IO