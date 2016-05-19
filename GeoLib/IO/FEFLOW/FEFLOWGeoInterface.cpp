/**
 * \copyright
 * Copyright (c) 2012-2016, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "FEFLOWGeoInterface.h"

#include <cctype>
#include <QtXml>

#include <logog/include/logog.hpp>

#include "BaseLib/FileTools.h"
#include "BaseLib/StringTools.h"

#include "GeoLib/GEOObjects.h"
#include "GeoLib/Point.h"
#include "GeoLib/Polygon.h"


namespace GeoLib
{
namespace IO
{

void FEFLOWGeoInterface::readFEFLOWFile(const std::string &filename, GeoLib::GEOObjects& geoObjects)
{
    std::ifstream in(filename.c_str());
    if (!in)
    {
        ERR("FEFLOWInterface::readFEFLOWFile(): Could not open file %s.", filename.c_str());
    }

    FEM_CLASS fem_class;
    std::vector<GeoLib::Point*>* points = NULL;
    std::vector<GeoLib::Polyline*>* lines = NULL;

    bool isXZplane = false;

    std::string line_string;
    std::stringstream line_stream;
    while (!in.eof())
    {
        getline(in, line_string);
        //....................................................................
        // CLASS
        if (line_string.find("CLASS") != std::string::npos) // the version number follows afterward, e.g. CLASS (v.5.313)
        {
            getline(in, line_string);
            line_stream.str(line_string);
            // problem class, time mode, problem orientation, dimension, nr. layers for 3D, saturation switch, precision of results, precision of coordinates
            line_stream >> fem_class.problem_class >> fem_class.time_mode >> fem_class.orientation >> fem_class.dimension >> fem_class.n_layers3d;
            line_stream.clear();
        }
        //....................................................................
        // GRAVITY
        else if (line_string.compare("GRAVITY") == 0)
        {
            getline(in, line_string);
            line_stream.str(line_string);
            double vec[3] = { };
            line_stream >> vec[0] >> vec[1] >> vec[2];
            if (vec[0] == 0.0 && vec[1] == -1.0 && vec[2] == 0.0)
                // x-z plane
                isXZplane = true;
            line_stream.clear();
        }
        //....................................................................
        // SUPERMESH
        else if (line_string.compare("SUPERMESH") == 0)
        {
            readSuperMesh(in, fem_class, &points, &lines);
        }
        //....................................................................
    }
    in.close();

    std::string project_name(
        BaseLib::extractBaseNameWithoutExtension(filename));
    if (points)
        geoObjects.addPointVec(
            std::unique_ptr<std::vector<GeoLib::Point*>>(points), project_name);
    if (lines)
        geoObjects.addPolylineVec(
            std::unique_ptr<std::vector<GeoLib::Polyline*>>(lines),
            project_name);

    if (isXZplane && points)
    {
        for (auto* pt : *points)
        {
            (*pt)[2] = (*pt)[1];
            (*pt)[1] = .0;
        }
    }

}


void FEFLOWGeoInterface::readPoints(QDomElement &nodesEle, const std::string &tag, int dim, std::vector<GeoLib::Point*> &points)
{
    QDomElement xmlEle = nodesEle.firstChildElement(QString::fromStdString(tag));
    if (xmlEle.isNull())
        return;
    QString str_pt_list1 = xmlEle.text();
    std::istringstream ss(str_pt_list1.toStdString());
    std::string line_str;
    while (!ss.eof())
    {
        std::getline(ss, line_str);
        BaseLib::trim(line_str, ' ');
        if (line_str.empty()) continue;
        std::istringstream line_ss(line_str);
        std::size_t pt_id = 0;
        std::array<double,3> pt_xyz;
        line_ss >> pt_id;
        for (int i = 0; i < dim; i++)
            line_ss >> pt_xyz[i];
        points[pt_id - 1] = new GeoLib::Point(pt_xyz, pt_id);
    }
}


void FEFLOWGeoInterface::readSuperMesh(std::ifstream &in, const FEM_CLASS &fem_class, std::vector<GeoLib::Point*>** p_points, std::vector<GeoLib::Polyline*>** p_lines)
{
    // get XML strings
    std::ostringstream oss;
    std::string line_string;
    while (true)
    {
        getline(in, line_string);
        BaseLib::trim(line_string);
        oss << line_string << "\n";
        if (line_string.find("</supermesh>") != std::string::npos)
            break;
    }
    const QString strXML(oss.str().c_str());

    // convert string to XML
    QDomDocument doc;
    if (!doc.setContent(strXML))
    {
        ERR("FEFLOWInterface::readSuperMesh(): Illegal XML format error");
        return;
    }

    // get geometry data from XML
    QDomElement docElem = doc.documentElement(); // #supermesh
    // #nodes
    *p_points = new std::vector<GeoLib::Point*>();
    std::vector<GeoLib::Point*>* points = *p_points;
    QDomElement nodesEle = docElem.firstChildElement("nodes");
    if (nodesEle.isNull())
        return;

    {
        const QString str = nodesEle.attribute("count");
        const long n_points = str.toLong();
        points->resize(n_points);
        //fixed
        readPoints(nodesEle, "fixed", fem_class.dimension, *points);
        readPoints(nodesEle, "linear", fem_class.dimension, *points);
        readPoints(nodesEle, "parabolic", fem_class.dimension, *points);
    }

    // #polygons
    *p_lines = new std::vector<GeoLib::Polyline*>();
    std::vector<GeoLib::Polyline*>* lines = *p_lines;
    QDomElement polygonsEle = docElem.firstChildElement("polygons");
    if (polygonsEle.isNull())
        return;

    {
        QDomNode child = polygonsEle.firstChild();
        while (!child.isNull())
        {
            if (child.nodeName() != "polygon")
            {
                child = child.nextSibling();
                continue;
            }
            QDomElement xmlEle = child.firstChildElement("nodes");
            if (xmlEle.isNull())
                continue;
            const QString str = xmlEle.attribute("count");
            const std::size_t n_points = str.toLong();
            QString str_ptId_list = xmlEle.text().simplified();
            {
                GeoLib::Polyline* line = new GeoLib::Polyline(*points);
                lines->push_back(line);
                std::istringstream ss(str_ptId_list.toStdString());
                for (std::size_t i = 0; i < n_points; i++)
                {
                    int pt_id = 0;
                    ss >> pt_id;
                    line->addPoint(pt_id - 1);
                }
                line->addPoint(line->getPointID(0));
            }
            child = child.nextSibling();
        }
    }
}

} // IO
} // GeoLib
