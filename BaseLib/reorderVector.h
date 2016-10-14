/**
 * \brief Reorder vector elements by given indices
 *
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file   reorderVector.h
 * Created on October 13, 2016, 5:37 PM
 */

#ifndef OGS_BASELIB_REORDERVECTOR_H
#define OGS_BASELIB_REORDERVECTOR_H

namespace BaseLib
{
/**
 * Reorder a vector by a given index vector.
 *  From
 *  <a href="reorderV">http://stackoverflow.com/questions/838384/reorder-vector-using-a-vector-of-indices</a>
 *
 *  Note: The simple version on that page is taken, which is good enough in performance
 *        for medium size vectors.
 */
template <typename ValueType, typename IndexType>
void reorderVector(std::vector<ValueType>& v,
                   std::vector<IndexType> const& order)
{
    for (std::size_t s = 1, d; s < order.size(); ++s)
    {
        for (d = order[s]; d < s; d = order[d]);
        if (d == s)
            while (d = order[d], d != s)
                std::swap(v[s], v[d]);
    }
}

} // end of namespace
#endif /* OGS_BASELIB_REORDERVECTOR_H */

