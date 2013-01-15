/**
 * \file   DiscreteObjectWithID.h
 * \author Norihiro Watanabe
 * \date   2012-08-03
 * \brief  Discrete object
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */
#pragma once

#ifndef DISCRETEOBJECTWITHID_H_
#define DISCRETEOBJECTWITHID_H_

namespace DiscreteLib
{

/**
 * \brief Object in discrete systems
 */
class DiscreteObjectWithID
{
public:
    ///
    DiscreteObjectWithID() : _obj_id(0) {};

    ///
    virtual ~DiscreteObjectWithID() {};

    /// return an object id
    std::size_t getObjectID() const {return _obj_id;};

    ///  set an object id
    void setObjectID(std::size_t i) {_obj_id = i;};
    
private:
    std::size_t _obj_id;
};

} // end

#endif //DISCRETEOBJECTWITHID_H_
