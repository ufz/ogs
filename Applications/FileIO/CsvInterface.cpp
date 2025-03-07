/**
 * \file
 * \author Karsten Rink
 * \date 2015-03-25
 * \brief Implementation of the CsvInterface class.
 *
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 */

#include "CsvInterface.h"

#include <algorithm>
#include <numeric>

#include "GeoLib/Point.h"

namespace FileIO
{
CsvInterface::CsvInterface() = default;

std::vector<std::string> CsvInterface::getColumnNames(std::string const& fname,
                                                      char const delim)
{
    std::ifstream in(fname.c_str());

    if (!in.is_open())
    {
        ERR("CsvInterface::getColumnNames(): Could not open file {:s}.", fname);
        return std::vector<std::string>();
    }
    std::string line;
    if (!std::getline(in, line))
    {
        ERR("CsvInterface::getColumnNames(): Could not read line from file "
            "{:s}. Is it empty?",
            fname);
        return {};
    }

    if (delim == '\n')
    {
        return {};
    }

    std::list<std::string> fields = BaseLib::splitString(line, delim);
    if (fields.size() < 2)
    {
        for (char const d : {'\t', ';', ','})
        {
            fields = BaseLib::splitString(line, d);
            if (fields.size() > 1)
            {
                break;
            }
        }
    }
    return {begin(fields), end(fields)};
}

int CsvInterface::readPoints(std::string const& fname, char delim,
                             std::vector<GeoLib::Point*>& points)
{
    std::ifstream in(fname.c_str());

    if (!in.is_open())
    {
        ERR("CsvInterface::readPoints(): Could not open file {:s}.", fname);
        return -1;
    }

    std::string line;
    std::getline(in, line);

    std::size_t line_count(0);
    std::size_t error_count(0);
    std::list<std::string>::const_iterator it;
    while (std::getline(in, line))
    {
        line_count++;
        std::list<std::string> const fields = BaseLib::splitString(line, delim);

        if (fields.size() < 3)
        {
            ERR("Line {:d} contains not enough columns of data. Skipping "
                "line...",
                line_count);
            error_count++;
            continue;
        }
        it = fields.begin();
        try
        {
            std::array<double, 3> point{};
            point[0] = std::stod(*it);
            point[1] = std::stod(*(++it));
            point[2] = std::stod(*(++it));
            points.push_back(new GeoLib::Point(point[0], point[1], point[2]));
        }
        catch (const std::invalid_argument&)
        {
            ERR("Error converting data to coordinates in line {:d}.",
                line_count);
        }
    }
    return error_count;
}

int CsvInterface::readPoints(std::string const& fname, char delim,
                             std::vector<GeoLib::Point*>& points,
                             std::string const& x_column_name,
                             std::string const& y_column_name,
                             std::string const& z_column_name)
{
    std::ifstream in(fname.c_str());
    std::array<std::string, 3> const column_names = {
        {x_column_name, y_column_name, z_column_name}};

    if (!in.is_open())
    {
        ERR("CsvInterface::readPoints(): Could not open file {:s}.", fname);
        return -1;
    }

    std::string line;
    std::getline(in, line);
    std::array<std::size_t, 3> const column_idx = {
        {CsvInterface::findColumn(line, delim, x_column_name),
         CsvInterface::findColumn(line, delim, y_column_name),
         (z_column_name.empty())
             ? CsvInterface::findColumn(line, delim, y_column_name)
             : CsvInterface::findColumn(line, delim, z_column_name)}};

    for (std::size_t i = 0; i < 3; ++i)
    {
        if (column_idx[i] == std::numeric_limits<std::size_t>::max())
        {
            ERR("Column '{:s}' not found in file header.", column_names[i]);
            return -1;
        }
    }

    return readPoints(in, delim, points, column_idx);
}

int CsvInterface::readPoints(std::string const& fname, char delim,
                             std::vector<GeoLib::Point*>& points,
                             std::size_t x_column_idx, std::size_t y_column_idx,
                             std::size_t z_column_idx)
{
    std::ifstream in(fname.c_str());

    if (!in.is_open())
    {
        ERR("CsvInterface::readPoints(): Could not open file {:s}.", fname);
        return -1;
    }

    if (z_column_idx == std::numeric_limits<std::size_t>::max())
    {
        z_column_idx = y_column_idx;
    }
    std::array<std::size_t, 3> const column_idx = {
        {x_column_idx, y_column_idx, z_column_idx}};

    return readPoints(in, delim, points, column_idx);
}

int CsvInterface::readPoints(std::ifstream& in, char delim,
                             std::vector<GeoLib::Point*>& points,
                             std::array<std::size_t, 3> const& column_idx)
{
    std::array<std::size_t, 3> order = {{0, 1, 2}};
    std::sort(order.begin(), order.end(),
              [&column_idx](std::size_t idx1, std::size_t idx2)
              { return column_idx[idx1] < column_idx[idx2]; });
    std::array<std::size_t, 3> const column_advance = {
        {column_idx[order[0]], column_idx[order[1]] - column_idx[order[0]],
         column_idx[order[2]] - column_idx[order[1]]}};

    std::string line;
    std::size_t line_count(0);
    std::size_t error_count(0);
    std::list<std::string>::const_iterator it;

    while (std::getline(in, line))
    {
        line_count++;
        std::list<std::string> const fields = BaseLib::splitString(line, delim);

        if (fields.size() < column_idx[order[2]] + 1)
        {
            ERR("Line {:d} contains not enough columns of data. Skipping "
                "line...",
                line_count);
            error_count++;
            continue;
        }

        it = fields.begin();
        try
        {
            std::advance(it, column_advance[0]);
            std::array<double, 3> point{};
            point[order[0]] = std::stod(*it);
            std::advance(it, column_advance[1]);
            point[order[1]] = std::stod(*it);
            std::advance(it, column_advance[2]);
            point[order[2]] =
                (column_idx[1] == column_idx[2]) ? 0 : std::stod(*it);
            points.push_back(new GeoLib::Point(point[0], point[1], point[2]));
        }
        catch (const std::invalid_argument&)
        {
            ERR("Error converting data to coordinates in line {:d}.",
                line_count);
            error_count++;
        }
    }
    return error_count;
}

std::size_t CsvInterface::findColumn(std::string const& line, char delim,
                                     std::string const& column_name)
{
    std::list<std::string> const fields = BaseLib::splitString(line, delim);
    if (fields.empty())
    {
        return std::numeric_limits<std::size_t>::max();
    }

    std::size_t count(0);
    for (const auto& field : fields)
    {
        if (field == column_name)
        {
            break;
        }

        count++;
    }

    if (count == fields.size())
    {
        return std::numeric_limits<std::size_t>::max();
    }

    return count;
}

void CsvInterface::addIndexVectorForWriting(std::size_t s)
{
    std::vector<int> idx_vec(s);
    std::iota(idx_vec.begin(), idx_vec.end(), 0);
    addVectorForWriting("Index", idx_vec);
}

bool CsvInterface::write()
{
    if (_data.empty())
    {
        ERR("CsvInterface::write() - No data to write.");
        return false;
    }

    std::size_t const n_vecs(_data.size());
    std::size_t const vec_size(getVectorSize(0));

    if (_writeCsvHeader)
    {
        out << _vec_names[0];
        for (std::size_t i = 1; i < n_vecs; ++i)
        {
            out << "\t" << _vec_names[i];
        }
        out << "\n";
    }

    for (std::size_t j = 0; j < vec_size; ++j)
    {
        writeValue(0, j);
        for (std::size_t i = 1; i < n_vecs; ++i)
        {
            out << "\t";
            writeValue(i, j);
        }
        out << "\n";
    }
    return true;
}

std::size_t CsvInterface::getVectorSize(std::size_t idx) const
{
    if (_data[idx].type() == typeid(std::vector<std::string>))
    {
        return std::any_cast<std::vector<std::string>>(_data[idx]).size();
    }
    if (_data[idx].type() == typeid(std::vector<double>))
    {
        return std::any_cast<std::vector<double>>(_data[idx]).size();
    }
    if (_data[idx].type() == typeid(std::vector<int>))
    {
        return std::any_cast<std::vector<int>>(_data[idx]).size();
    }
    return 0;
}

void CsvInterface::writeValue(std::size_t vec_idx, std::size_t in_vec_idx)
{
    if (_data[vec_idx].type() == typeid(std::vector<std::string>))
    {
        out << std::any_cast<std::vector<std::string>>(
            _data[vec_idx])[in_vec_idx];
    }
    else if (_data[vec_idx].type() == typeid(std::vector<double>))
    {
        out << std::any_cast<std::vector<double>>(_data[vec_idx])[in_vec_idx];
    }
    else if (_data[vec_idx].type() == typeid(std::vector<int>))
    {
        out << std::any_cast<std::vector<int>>(_data[vec_idx])[in_vec_idx];
    }
}

}  // end namespace FileIO
