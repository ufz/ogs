/**
 * @date 2014-09-19
 * @brief Base class for different search length strategies.
 *
 * @copyright
 * Copyright (c) 2012-2017, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/LICENSE.txt
 */

#ifndef SEARCHLENGTH_H_
#define SEARCHLENGTH_H_

namespace MeshGeoToolsLib
{

/// Base class for all SearchLength strategy implementations.
/// The default implementation is mesh independent and provides a strong
/// criterion for searching mesh nodes near a geometry. The algorithm can be
/// used for meshes that have nearly equi-sized elements.
class SearchLength
{
public:
    /// Constructor for SearchLength object with a default search length
    /// of 10 angstrom (\f$10^{-9}\f$ m)
    explicit SearchLength(double search_length = 1e-9)
        : _search_length(search_length) {}

    SearchLength(SearchLength const&) = default;
    SearchLength& operator=(SearchLength const&) = default;

    virtual ~SearchLength() = default;

    virtual double getSearchLength() const
    {
        return _search_length;
    }

protected:
    double _search_length;
};

} // end namespace MeshGeoToolsLib

#endif

