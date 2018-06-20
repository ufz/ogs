/**
 * \file   BoostXmlGmlInterface.cpp
 * \author Karsten Rink
 * \date   2014-01-31
 * \brief  Implementation of the BoostXmlGmlInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "BoostXmlGmlInterface.h"

#include <boost/property_tree/xml_parser.hpp>
#include <logog/include/logog.hpp>

#include "BaseLib/ConfigTreeUtil.h"
#include "BaseLib/uniqueInsert.h"
#include "GeoLib/GEOObjects.h"
#include "GeoLib/Point.h"
#include "GeoLib/PointVec.h"
#include "GeoLib/Polyline.h"
#include "GeoLib/Surface.h"
#include "GeoLib/Triangle.h"

namespace GeoLib
{
namespace IO
{

BoostXmlGmlInterface::BoostXmlGmlInterface(GeoLib::GEOObjects& geo_objs) :
        _geo_objects(geo_objs)
{}

bool BoostXmlGmlInterface::readFile(const std::string &fname)
{
    //! \todo Reading geometries is always strict.
    auto doc = BaseLib::makeConfigTree(fname, true, "OpenGeoSysGLI");

    // ignore attributes related to XML schema
    doc->ignoreConfigAttribute("xmlns:xsi");
    doc->ignoreConfigAttribute("xsi:noNamespaceSchemaLocation");
    doc->ignoreConfigAttribute("xmlns:ogs");

    auto points = std::make_unique<std::vector<GeoLib::Point*>>();
    auto polylines = std::make_unique<std::vector<GeoLib::Polyline*>>();
    auto surfaces = std::make_unique<std::vector<GeoLib::Surface*>>();

    using MapNameId = std::map<std::string, std::size_t>;
    auto pnt_names = std::make_unique<MapNameId>();
    auto ply_names = std::make_unique<MapNameId>();
    auto sfc_names = std::make_unique<MapNameId>();

    //! \ogs_file_param{gml__name}
    auto geo_name = doc->getConfigParameter<std::string>("name");
    if (geo_name.empty())
    {
        OGS_FATAL("BoostXmlGmlInterface::readFile(): <name> tag is empty.");
    }

    //! \ogs_file_param{gml__points}
    for (auto st : doc->getConfigSubtreeList("points"))
    {
        readPoints(st, *points, *pnt_names);
        _geo_objects.addPointVec(std::move(points), geo_name,
                                 std::move(pnt_names));
    }

    //! \ogs_file_param{gml__polylines}
    for (auto st : doc->getConfigSubtreeList("polylines"))
    {
        readPolylines(st,
                      *polylines,
                      *_geo_objects.getPointVec(geo_name),
                      _geo_objects.getPointVecObj(geo_name)->getIDMap(),
                      *ply_names);
    }

    //! \ogs_file_param{gml__surfaces}
    for (auto st : doc->getConfigSubtreeList("surfaces"))
    {
        readSurfaces(st,
                     *surfaces,
                     *_geo_objects.getPointVec(geo_name),
                     _geo_objects.getPointVecObj(geo_name)->getIDMap(),
                     *sfc_names);
    }

    if (!polylines->empty()) {
        _geo_objects.addPolylineVec(std::move(polylines), geo_name,
                                    std::move(ply_names));
    }

    if (!surfaces->empty()) {
        _geo_objects.addSurfaceVec(std::move(surfaces), geo_name,
                                   std::move(sfc_names));
    }

    return true;
}

void BoostXmlGmlInterface::readPoints(BaseLib::ConfigTree const& pointsRoot,
                                      std::vector<GeoLib::Point*>& points,
                                      std::map<std::string, std::size_t>& pnt_names )
{
    //! \ogs_file_param{gml__points__point}
    for (auto const pt : pointsRoot.getConfigParameterList("point"))
    {
        //! \ogs_file_attr{gml__points__point__id}
        auto const p_id = pt.getConfigAttribute<std::size_t>("id");
        //! \ogs_file_attr{gml__points__point__x}
        auto const p_x  = pt.getConfigAttribute<double>("x");
        //! \ogs_file_attr{gml__points__point__y}
        auto const p_y  = pt.getConfigAttribute<double>("y");
        //! \ogs_file_attr{gml__points__point__z}
        auto const p_z  = pt.getConfigAttribute<double>("z");

        auto const p_size = points.size();
        BaseLib::insertIfKeyUniqueElseError(_idx_map, p_id, p_size,
            "The point id is not unique.");
        points.push_back(new GeoLib::Point(p_x, p_y, p_z, p_id));

        //! \ogs_file_attr{gml__points__point__name}
        if (auto const p_name = pt.getConfigAttributeOptional<std::string>("name"))
        {
            if (p_name->empty()) {
                OGS_FATAL("Empty point name found in geometry file.");
            }

            BaseLib::insertIfKeyUniqueElseError(pnt_names, *p_name, p_size,
                "The point name is not unique.");
        }
    }
}

void BoostXmlGmlInterface::readPolylines(
    BaseLib::ConfigTree const& polylinesRoot,
    std::vector<GeoLib::Polyline*>& polylines,
    std::vector<GeoLib::Point*> const& points,
    std::vector<std::size_t> const& pnt_id_map,
    std::map<std::string, std::size_t>& ply_names)
{
    //! \ogs_file_param{gml__polylines__polyline}
    for (auto const pl : polylinesRoot.getConfigSubtreeList("polyline"))
    {
        //! \ogs_file_attr{gml__polylines__polyline__id}
        auto const id = pl.getConfigAttribute<std::size_t>("id");
        // The id is not used but must be present in the GML file.
        // That's why pl.ignore...() cannot be used.
        (void) id;

        polylines.push_back(new GeoLib::Polyline(points));

        //! \ogs_file_attr{gml__polylines__polyline__name}
        if (auto const p_name = pl.getConfigAttributeOptional<std::string>("name"))
        {
            if (p_name->empty()) {
                OGS_FATAL("Empty polyline name found in geometry file.");
            }

            BaseLib::insertIfKeyUniqueElseError(ply_names, *p_name, polylines.size()-1,
                "The polyline name is not unique.");

            //! \ogs_file_param{gml__polylines__polyline__pnt}
            for (auto const pt : pl.getConfigParameterList<std::size_t>("pnt")) {
                polylines.back()->addPoint(pnt_id_map[_idx_map[pt]]);
            }
        }
        else
        {
            // polyline has no name, ignore it.
            pl.ignoreConfigParameterAll("pnt");
        }
    }
}

void BoostXmlGmlInterface::readSurfaces(
    BaseLib::ConfigTree const&  surfacesRoot,
    std::vector<GeoLib::Surface*>& surfaces,
    std::vector<GeoLib::Point*> const& points,
    const std::vector<std::size_t>& pnt_id_map,
    std::map<std::string, std::size_t>& sfc_names)
{
    //! \ogs_file_param{gml__surfaces__surface}
    for (auto const& sfc : surfacesRoot.getConfigSubtreeList("surface"))
    {
        //! \ogs_file_attr{gml__surfaces__surface__id}
        auto const id = sfc.getConfigAttribute<std::size_t>("id");
        // The id is not used but must be present in the GML file.
        // That's why sfc.ignore...() cannot be used.
        (void) id;
        surfaces.push_back(new GeoLib::Surface(points));

        //! \ogs_file_attr{gml__surfaces__surface__name}
        if (auto const s_name = sfc.getConfigAttributeOptional<std::string>("name"))
        {
            if (s_name->empty()) {
                OGS_FATAL("Empty surface name found in geometry file.");
            }

            BaseLib::insertIfKeyUniqueElseError(sfc_names, *s_name, surfaces.size()-1,
                "The surface name is not unique.");

            //! \ogs_file_param{gml__surfaces__surface__element}
            for (auto const& element : sfc.getConfigParameterList("element")) {
                //! \ogs_file_attr{gml__surfaces__surface__element__p1}
                auto const p1_attr = element.getConfigAttribute<std::size_t>("p1");
                //! \ogs_file_attr{gml__surfaces__surface__element__p2}
                auto const p2_attr = element.getConfigAttribute<std::size_t>("p2");
                //! \ogs_file_attr{gml__surfaces__surface__element__p3}
                auto const p3_attr = element.getConfigAttribute<std::size_t>("p3");

                auto const p1 = pnt_id_map[_idx_map[p1_attr]];
                auto const p2 = pnt_id_map[_idx_map[p2_attr]];
                auto const p3 = pnt_id_map[_idx_map[p3_attr]];
                surfaces.back()->addTriangle(p1,p2,p3);
            }
        }
        else
        {
            // surface has no name, ignore it.
            sfc.ignoreConfigParameterAll("element");
        }
    }
}

bool BoostXmlGmlInterface::write()
{
    if (this->_exportName.empty()) {
        ERR("BoostXmlGmlInterface::write(): No geometry specified.");
        return false;
    }

    GeoLib::PointVec const*const pnt_vec(_geo_objects.getPointVecObj(_exportName));
    if (! pnt_vec) {
        ERR("BoostXmlGmlInterface::write(): No PointVec within the geometry \"%s\".",
            _exportName.c_str());
        return false;
    }

    std::vector<GeoLib::Point*> const*const pnts(pnt_vec->getVector());
    if (! pnts) {
        ERR("BoostXmlGmlInterface::write(): No vector of points within the geometry \"%s\".",
            _exportName.c_str());
        return false;
    }
    if (pnts->empty()) {
        ERR("BoostXmlGmlInterface::write(): No points within the geometry \"%s\".",
            _exportName.c_str());
        return false;
    }

    // create a property tree for writing it to file
    boost::property_tree::ptree pt;

    // put header in property tree
    pt.put("<xmlattr>.xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance");
    pt.put("<xmlattr>.xsi:noNamespaceSchemaLocation",
        "https://www.opengeosys.org/images/xsd/OpenGeoSysGLI.xsd");
    pt.put("<xmlattr>.xmlns:ogs", "http://www.opengeosys.net");
    auto& geometry_set = pt.add("OpenGeoSysGLI", "");

    geometry_set.add("name", _exportName);
    auto& pnts_tag = geometry_set.add("points", "");
    for (std::size_t k(0); k<pnts->size(); k++) {
        auto& pnt_tag = pnts_tag.add("point", "");
        pnt_tag.put("<xmlattr>.id", k);
        pnt_tag.put("<xmlattr>.x", (*((*pnts)[k]))[0]);
        pnt_tag.put("<xmlattr>.y", (*((*pnts)[k]))[1]);
        pnt_tag.put("<xmlattr>.z", (*((*pnts)[k]))[2]);
        std::string const& point_name(pnt_vec->getItemNameByID(k));
        if (!point_name.empty())
            pnt_tag.put("<xmlattr>.name", point_name);
    }

    addPolylinesToPropertyTree(geometry_set);
    addSurfacesToPropertyTree(geometry_set);

    boost::property_tree::xml_writer_settings<std::string> settings('\t', 1);
    setPrecision(std::numeric_limits<double>::digits10);
    write_xml(_out, pt, settings);
    return true;
}

void BoostXmlGmlInterface::addSurfacesToPropertyTree(
    boost::property_tree::ptree & geometry_set)
{
    GeoLib::SurfaceVec const*const sfc_vec(_geo_objects.getSurfaceVecObj(_exportName));
    if (!sfc_vec) {
        INFO("BoostXmlGmlInterface::addSurfacesToPropertyTree(): "
            "No surfaces within the geometry \"%s\".", _exportName.c_str());
        return;
    }

    std::vector<GeoLib::Surface*> const*const surfaces(sfc_vec->getVector());
    if (!surfaces || surfaces->empty())
    {
        INFO(
            "BoostXmlGmlInterface::addSurfacesToPropertyTree(): "
            "No surfaces within the geometry \"%s\".",
            _exportName.c_str());
        return;
    }

    auto& surfaces_tag = geometry_set.add("surfaces", "");
    for (std::size_t i=0; i<surfaces->size(); ++i) {
        GeoLib::Surface const*const surface((*surfaces)[i]);
        std::string sfc_name("");
        sfc_vec->getNameOfElement(surface, sfc_name);
        auto& surface_tag = surfaces_tag.add("surface", "");
        surface_tag.put("<xmlattr>.id", i);
        if (!sfc_name.empty())
            surface_tag.put("<xmlattr>.name", sfc_name);
        for (std::size_t j=0; j<surface->getNumberOfTriangles(); ++j) {
            auto& element_tag = surface_tag.add("element", "");
            element_tag.put("<xmlattr>.p1", (*(*surface)[j])[0]);
            element_tag.put("<xmlattr>.p2", (*(*surface)[j])[1]);
            element_tag.put("<xmlattr>.p3", (*(*surface)[j])[2]);
        }
    }
}

void BoostXmlGmlInterface::addPolylinesToPropertyTree(
    boost::property_tree::ptree & geometry_set)
{
    GeoLib::PolylineVec const*const vec(_geo_objects.getPolylineVecObj(_exportName));
    if (!vec) {
        INFO("BoostXmlGmlInterface::addPolylinesToPropertyTree(): "
            "No polylines within the geometry \"%s\".", _exportName.c_str());
        return;
    }

    std::vector<GeoLib::Polyline*> const*const polylines(vec->getVector());
    if (!polylines || polylines->empty())
    {
        INFO(
            "BoostXmlGmlInterface::addPolylinesToPropertyTree(): "
            "No polylines within the geometry \"%s\".",
            _exportName.c_str());
        return;
    }

    auto& polylines_tag = geometry_set.add("polylines", "");
    for (std::size_t i=0; i<polylines->size(); ++i) {
        GeoLib::Polyline const*const polyline((*polylines)[i]);
        std::string ply_name("");
        vec->getNameOfElement(polyline, ply_name);
        auto& polyline_tag = polylines_tag.add("polyline", "");
        polyline_tag.put("<xmlattr>.id", i);
        if (!ply_name.empty())
            polyline_tag.put("<xmlattr>.name", ply_name);
        for (std::size_t j=0; j<polyline->getNumberOfPoints(); ++j) {
            polyline_tag.add("pnt", polyline->getPointID(j));
        }
    }
}

} // end namespace IO
} // end namespace GeoLib
