/**
 * \file
 * \author Norihiro Watanabe
 * \date   2013-04-16
 * \brief
 *
 * \copyright
 * Copyright (c) 2013, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */


#ifndef MESHSUBSET_H_
#define MESHSUBSET_H_

#include <vector>
#include <algorithm>

#include "MeshLib/Mesh.h"
#include "MeshLib/Node.h"
#include "MeshLib/Elements/Element.h"


namespace VecMatOnMeshLib
{

/**
 * A subset of nodes or elements on a single mesh.
 */
class MeshSubset
{
public:
    /// construct from nodes
    MeshSubset(const MeshLib::Mesh& msh, std::vector<MeshLib::Node*> const& vec_items)
    : _msh(msh), _nodes(&vec_items), _eles(nullptr)
    {}

    /// construct from elements
    MeshSubset(const MeshLib::Mesh& msh, std::vector<MeshLib::Element*> const& vec_items)
    : _msh(msh), _nodes(nullptr), _eles(&vec_items)
    {}

    /// construct from both nodes and elements
    MeshSubset(const MeshLib::Mesh& msh,
                std::vector<MeshLib::Node*> const& vec_nodes, std::vector<MeshLib::Element*> const& vec_eles)
    : _msh(msh), _nodes(&vec_nodes), _eles(&vec_eles)
    {}

    ~MeshSubset() {};

    /// return the total number of mesh items
    std::size_t getNTotalItems() const { return getNNodes() + getNElements(); }

    /// return this mesh ID
    std::size_t getMeshID() const { return _msh.getID(); }

    /// return the number of registered nodes
    std::size_t getNNodes() const { return (_nodes==nullptr) ? 0 : _nodes->size(); }

    /// return the number of registered elements
    std::size_t getNElements() const { return (_eles==nullptr) ? 0 : _eles->size(); }


//    const MeshLib::Node* getNode(std::size_t index) const { return (_nodes==nullptr) ? nullptr : (*_nodes)[index]; }

//    const MeshLib::Element* getElement(std::size_t index) const { return (_eles==nullptr) ? nullptr : (*_eles)[index]; }

//    bool has(const MeshLib::Node &item) const { return (std::count(_nodes->begin(), _nodes->end(), &item)>0); }

//    bool has(const MeshLib::Element &item) const { return (std::count(_eles->begin(), _eles->end(), &item)>0); }


private:
    const MeshLib::Mesh& _msh;
    std::vector<MeshLib::Node*> const* _nodes;
    std::vector<MeshLib::Element*> const* _eles;

};

}

#endif /* MESHSUBSET_H_ */
