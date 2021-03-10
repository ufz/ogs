/**
 * \file
 * \author Dmitrij Naumov
 * \brief Implementation of Histogram class.
 *
 * \copyright
 * Copyright (c) 2012-2021, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#pragma once

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iosfwd>
#include <iterator>
#include <utility>
#include <vector>

#include "Logging.h"

namespace BaseLib
{
/** Basic Histogram implementation.
 *
 * Creates histogram from input data of type \c T.
 */
template <typename T>
class Histogram
{
public:
    using Data =
        typename std::vector<double>;  /// Underlying input data vector type.

public:
    /** Creates histogram of the given element in the range \c [first, last).
     *
     * Input data is copied into \c std::vector.
     *
     * \param first Range of elements to create histogram from.
     * \param last Range of elements to create histogram from.
     * \param nr_bins Number of bins in histogram.
     * \param computeHistogram Compute histogram if set. If not set user must
     *        call \c update() before accessing data.
     */
    template <typename InputIterator>
    Histogram(InputIterator first, InputIterator last, const int nr_bins = 16,
              const bool computeHistogram = true )
        : data_(first, last), nr_bins_(nr_bins), dirty_(true)
    {
        init(computeHistogram);
    }

    /** Creates histogram from \c std::vector.
     * \param data Input vector.
     * \param nr_bins Number of bins in histogram.
     * \param computeHistogram Compute histogram if set. If not set user must call
     * \c update() before accessing data.
     */
    explicit Histogram(std::vector<T> data, const unsigned int nr_bins = 16,
                       const bool computeHistogram = true)
        : data_(std::move(data)), nr_bins_(nr_bins), dirty_(true)
    {
        init(computeHistogram);
    }

    /** Updates histogram using sorted \c data_ vector.
     *
     * Start histogram creation with first element. Then find first element in
     * the next histogram bin. Number of elements in the bin is the difference
     * between these two iterators.
     * \verbatim
       [0.1, 0.2, ..., 0.7 , ..., 0.7+binWidth = 0.9,  1.0  , ..., last]
                       it             itEnd - 1      itEnd
       \endverbatim
     */
    void update()
    {
        if (!dirty_)
        {
            return;
        }

        bin_width_ = (max_ - min_) / nr_bins_;

        auto it = data_.begin();
        for (unsigned int bin = 0; bin < nr_bins_; bin++)
        {
            auto itEnd = std::upper_bound(it, data_.end(),
                                          min_ + (bin + 1) * bin_width_);
            histogram_[bin] = std::distance(it, itEnd);
            it = itEnd;
        }
        dirty_ = false;
    }

    void setMinimum(const T& minimum) { min_ = minimum; dirty_ = true; }
    void setMaximum(const T& maximum) { max_ = maximum; dirty_ = true; }

    const Data& getSortedData() const { return data_; }
    const std::vector<std::size_t>& getBinCounts() const { return histogram_; }
    const unsigned int& getNumberOfBins() const { return nr_bins_; }
    const T& getMinimum() const { return min_; }
    const T& getMaximum() const { return max_; }
    const T& getBinWidth() const { return bin_width_; }

    void
    prettyPrint(std::ostream& os, const unsigned int line_width = 16) const
    {
        const std::size_t count_max =
                *std::max_element(histogram_.begin(), histogram_.end());
        for (unsigned int bin = 0; bin < nr_bins_; ++bin)
        {
            os << "[" << min_ + bin * bin_width_ << ", " << min_ +
                    (bin + 1) * bin_width_ << ")\t";
            os << histogram_[bin] << "\t";

            const int n_stars =
                    std::ceil(line_width * ((double)histogram_[bin] / count_max));
            for (int star = 0; star < n_stars; star++)
            {
                os << "*";
            }
            os << "\n";
        }
    }

    int write(std::string const& file_name, std::string const& data_set_name, std::string const& param_name) const
    {
        if (file_name.empty())
        {
            ERR ("No file name specified.");
            return 1;
        }

        std::ofstream out (file_name);
        if (!out)
        {
            ERR("Error writing histogram: Could not open file.");
            return 1;
        }

        out << "# Histogram for parameter " << param_name << " of data set " << data_set_name << "\n";
        std::size_t const n_bins = this->getNumberOfBins();
        std::vector<std::size_t> const& bin_cnts(this->getBinCounts());
        double const min (this->getMinimum());
        double const bin_width (this->getBinWidth());

        for (std::size_t k(0); k < n_bins; k++)
        {
            out << min + k * bin_width << " " << bin_cnts[k] << "\n";
        }
        out.close ();
        return 0;
    }

protected:
    /** Initialize class members after constructor call.
     */
    void init(const bool computeHistogram = true)
    {
        std::sort(data_.begin(), data_.end());
        histogram_.resize(nr_bins_);
        min_ = data_.front();
        max_ = data_.back();
        bin_width_ = (max_ - min_) / nr_bins_;

        dirty_ = true;
        if (computeHistogram)
        {
            update();
        }
    }

    Data data_;
    const unsigned int nr_bins_;
    std::vector<std::size_t> histogram_;
    T min_, max_; ///< Minimum and maximum input data values.
    T bin_width_;

private:
    bool dirty_; ///< When set \c update() will recompute histogram.
};

/** Writes histogram to output stream.
 *
 * Writes histogram properties in this order:
 * number of bins, minimum, maximum, bin0 count, ..., binN-1 count.
 */
template <typename T>
std::ostream&
operator<<(std::ostream& os, const Histogram<T>& h)
{
    os << h.getNumberOfBins() << " "
       << h.getMinimum() << " "
       << h.getMaximum() << " ";
    std::copy(h.getBinCounts().begin(), h.getBinCounts().end(),
              std::ostream_iterator<T>(os, " "));
    return os << std::endl;
}
}   // namespace BaseLib
