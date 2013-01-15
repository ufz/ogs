/**
 * \file
 * \author Norihiro Watanabe
 * \date   2012-02-17
 * \brief Helper macros.
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef CODINGTOOLS_H
#define CODINGTOOLS_H

#define DISALLOW_COPY_AND_ASSIGN(TypeName) \
    TypeName(const TypeName&);   \
    TypeName &operator=(const TypeName&)

namespace BaseLib
{

template <typename T>
static void releaseObjectsInStdVector(T &object)
{
    if (object.size()>0) {
        const std::size_t vec_size(object.size());
        for (std::size_t i=0; i<vec_size; i++)
            if (object[i]!=nullptr) delete object[i];
        object.clear();
    }
};

template <typename T>
static void releaseObjectsInStdMap(T &object)
{
    if (object.size()>0) {
        for (typename T::iterator itr=object.begin(); itr!=object.end(); ++itr)
            if (itr->second!=nullptr) delete itr->second;
        object.clear();
    }
};

template <typename T>
static void releaseObjectsInStdQueue(T &container)
{
    while (!container.empty()) {
        if (container.front() !=nullptr) delete container.front();
        container.pop();
    }
};

} //BaseLib

#endif
