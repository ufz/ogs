/**
 * \file
 * \copyright
 * Copyright (c) 2012-2024, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "Histogram.h"

#include <cmath>
#include <fstream>

#include "BaseLib/Logging.h"

namespace BaseLib
{

template <typename T>
int Histogram<T>::write(std::string const& file_name,
                        std::string const& data_set_name,
                        std::string const& param_name) const
{
    if (file_name.empty())
    {
        ERR("No file name specified.");
        return 1;
    }

    std::ofstream out(file_name);
    if (!out)
    {
        ERR("Error writing histogram: Could not open file.");
        return 1;
    }

    out << "# Histogram for parameter " << param_name << " of data set "
        << data_set_name << "\n";
    std::size_t const n_bins = this->getNumberOfBins();
    std::vector<std::size_t> const& bin_cnts(this->getBinCounts());
    double const min(this->getMinimum());
    double const bin_width(this->getBinWidth());

    for (std::size_t k(0); k < n_bins; k++)
    {
        out << min + k * bin_width << " " << bin_cnts[k] << "\n";
    }
    out.close();
    return 0;
}

template <typename T>
void Histogram<T>::prettyPrint(std::ostream& os,
                               const unsigned int line_width) const
{
    const std::size_t count_max =
        *std::max_element(histogram_.begin(), histogram_.end());
    for (unsigned int bin = 0; bin < nr_bins_; ++bin)
    {
        os << "[" << min_ + bin * bin_width_ << ", "
           << min_ + (bin + 1) * bin_width_ << ")\t";
        os << histogram_[bin] << "\t";

        const int n_stars = static_cast<int>(
            std::ceil(line_width * ((double)histogram_[bin] / count_max)));
        for (int star = 0; star < n_stars; star++)
        {
            os << "*";
        }
        os << "\n";
    }
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const Histogram<T>& h)
{
    os << h.getNumberOfBins() << " " << h.getMinimum() << " " << h.getMaximum()
       << " ";
    std::copy(h.getBinCounts().begin(), h.getBinCounts().end(),
              std::ostream_iterator<std::size_t>(os, " "));
    return os << std::endl;
}

template class Histogram<double>;
template std::ostream& operator<<(std::ostream& os, const Histogram<double>& h);
}  // namespace BaseLib
