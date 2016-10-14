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
 *  Reorder a vector by a given index vector.
 *
 *  Note: It is good enough in performance for medium size vectors.
 */
template <typename ValueType, typename IndexType>
void reorderVector(std::vector<ValueType>& v,
                   std::vector<IndexType> const& order)
{
    std::vector<ValueType> temp_v = v;

    for (std::size_t i=0; i<order.size(); i++)
    {
        v[i] =  temp_v[order[i]];
    }
}

} // end of namespace
#endif /* OGS_BASELIB_REORDERVECTOR_H */

