#include "XdmfWriter.h"

#include <fstream>

#include "BaseLib/Logging.h"
#include "BaseLib/RunTime.h"
#include "writeXdmf.h"

namespace MeshLib::IO
{
XdmfWriter::XdmfWriter(std::string const& xdmf_filename,
                       std::function<std::string(std::vector<double>)>
                           xdmf_writer_fn)
    : filename(xdmf_filename), xdmf_writer(xdmf_writer_fn)
{
}

XdmfWriter::~XdmfWriter()
{
    BaseLib::RunTime time_output;
    time_output.start();
    std::ofstream fout;
    fout.open(this->filename);
    fout << xdmf_writer(times);

    INFO("[time] Output of XDMF took {:g} s.", time_output.elapsed());
}

void XdmfWriter::addTimeStep(double const& time_step)
{
    times.push_back(time_step);
}

}  // namespace MeshLib::IO