/**
 * \copyright
 * Copyright (c) 2012-2015, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#ifndef PROCESS_LIB_NEUMANN_BC_CONFIG_H_
#define PROCESS_LIB_NEUMANN_BC_CONFIG_H_

#include <boost/property_tree/ptree.hpp>
extern template class boost::property_tree::basic_ptree<std::basic_string<char>,
      std::basic_string<char>, std::less<std::basic_string<char> > >;

#include <logog/include/logog.hpp>

#include "MathLib/ConstantFunction.h"
#include "MeshGeoToolsLib/BoundaryElementsSearcher.h"
#include "MeshLib/Elements/Element.h"

namespace GeoLib
{
    class GeoObject;
}

namespace ProcessLib
{

class BoundaryConditionConfig
{
public:
    BoundaryConditionConfig(GeoLib::GeoObject const* const geometry)
        : _geometry(geometry)
    { }

    virtual ~BoundaryConditionConfig() = default;

protected:
    GeoLib::GeoObject const* const _geometry;
};


/// Configuration of a Neumann type boundary condition read from input file.
class NeumannBcConfig : public BoundaryConditionConfig
{
    using ConfigTree = boost::property_tree::ptree;
public:
    NeumannBcConfig(GeoLib::GeoObject const* const geometry,
            ConfigTree const& config)
        : BoundaryConditionConfig(geometry)
    {
        DBUG("Constructing NeumannBcConfig from config.");

        double const value = config.get<double>("value", 0);
        DBUG("Using value %g", value);

        _function = new MathLib::ConstantFunction<double>(value);
    }

    ~NeumannBcConfig()
    {
        for (auto e : _elements)
            delete e;

        delete _function;
    }

    /// Initialize Neumann type boundary conditions.
    /// Fills in elements of the particular geometry of the boundary condition.
    /// The elements are appended to the \c elements vector.
    void initialize(MeshGeoToolsLib::BoundaryElementsSearcher& searcher)
    {
        std::vector<MeshLib::Element*> elems = searcher.getBoundaryElements(*_geometry);

        // deep copy because the searcher destroys the elements.
        std::transform(elems.cbegin(), elems.cend(),
                std::back_inserter(_elements),
                std::mem_fn(&MeshLib::Element::clone));
    }

    /// Iterator over elements of the boundary condition.
    std::vector<MeshLib::Element*>::const_iterator
    elementsBegin() const
    {
        return _elements.cbegin();
    }

    /// Past the end iterator over elements of the boundary condition.
    /// \sa elementsBegin().
    std::vector<MeshLib::Element*>::const_iterator
    elementsEnd() const
    {
        return _elements.cend();
    }

    /// Returns the right-hand-side function of the boundary condition.
    MathLib::ConstantFunction<double>*
    getFunction() const
    {
        return _function;
    }

private:
    /// Domain of the boundary condition.
    std::vector<MeshLib::Element*> _elements;

    /// A function given on the domain of the boundary condition.
    MathLib::ConstantFunction<double>* _function = nullptr;
};

}   // namespace ProcessLib

#endif  // PROCESS_LIB_NEUMANN_BC_CONFIG_H_
