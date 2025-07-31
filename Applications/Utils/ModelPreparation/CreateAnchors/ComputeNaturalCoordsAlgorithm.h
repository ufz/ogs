/**
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vtkCell.h>
#include <vtkLine.h>
#include <vtkPointData.h>

#include <algorithm>

#include "ComputeNaturalCoordsSolver.h"
#include "FindCellsForPoints.h"
#include "MathLib/FormattingUtils.h"
#include "MeshLib/IO/VtkIO/VtkMeshConverter.h"
#include "SolverByElementTypeRegistry.h"

namespace ApplicationUtils
{
struct ComputeNaturalCoordsResult
{
    Eigen::MatrixXd natural_coords;
    Eigen::MatrixXd real_coords;
    Eigen::VectorXd initial_anchor_stress;
    Eigen::VectorXd maximum_anchor_stress;
    Eigen::VectorXd residual_anchor_stress;
    Eigen::VectorXd anchor_cross_sectional_area;
    Eigen::VectorXd anchor_stiffness;
    Eigen::VectorX<vtkIdType> bulk_element_ids;
    Eigen::VectorX<vtkIdType> point_cloud_node_ids;
    bool success;
};

class ComputeNaturalCoordsIntermediateResult
{
public:
    void append(Eigen::Vector3d const& natural_coords,
                Eigen::Vector3d const& real_coords,
                vtkIdType bulk_element_id,
                vtkIdType point_cloud_node_id)
    {
        natural_coordss.push_back(natural_coords);
        real_coordss.push_back(real_coords);
        bulk_element_ids.push_back(bulk_element_id);
        point_cloud_node_ids.push_back(point_cloud_node_id);
    }

    void fail() { success = false; }

    ComputeNaturalCoordsResult finished(
        Eigen::VectorXd const& initial_anchor_stress,
        Eigen::VectorXd const& maximum_anchor_stress,
        Eigen::VectorXd const& residual_anchor_stress,
        Eigen::VectorXd const& anchor_cross_sectional_area,
        Eigen::VectorXd const& anchor_stiffness)
    {
        auto bei = Eigen::Map<Eigen::VectorX<vtkIdType>>(
            bulk_element_ids.data(), bulk_element_ids.size());
        auto pcni = Eigen::Map<Eigen::VectorX<vtkIdType>>(
            point_cloud_node_ids.data(), point_cloud_node_ids.size());
        return {conv(natural_coordss),
                conv(real_coordss),
                initial_anchor_stress,
                maximum_anchor_stress,
                residual_anchor_stress,
                anchor_cross_sectional_area,
                anchor_stiffness,
                bei,
                pcni,
                success};
    }

private:
    static Eigen::MatrixXd conv(std::vector<Eigen::Vector3d> const& in)
    {
        Eigen::MatrixXd out(in.size(), 3);
        for (std::size_t i = 0; i < in.size(); ++i)
        {
            out.row(i).noalias() = in[i].transpose();
        }
        return out;
    }

    std::vector<Eigen::Vector3d> natural_coordss;
    std::vector<Eigen::Vector3d> real_coordss;
    std::vector<vtkIdType> bulk_element_ids;
    std::vector<vtkIdType> point_cloud_node_ids;
    bool success = true;
};

ComputeNaturalCoordsResult computeNaturalCoords(
    vtkUnstructuredGrid* const bulk_mesh,
    ComputeNaturalCoordsResult const& json_data, double const tolerance,
    int const max_iter)
{
    // Check inputs ------------------------------------------------------------

    if (json_data.real_coords.cols() != 3)
    {
        OGS_FATAL("Wrong number of input coordinates");
    }

    // General definitions -----------------------------------------------------
    using SolverRegistry =
        SolverByElementTypeRegistry<ComputeNaturalCoordsSolverInterface,
                                    ComputeNaturalCoordsSolverImplementation>;

    double const real_coords_tolerance = tolerance;

    // Initialization ----------------------------------------------------------
    FindCellsForPoints findCellsForPoints;
    findCellsForPoints.initialize(bulk_mesh);

    ComputeNaturalCoordsIntermediateResult result;

    // Do computation for each point ------------------------------------------
    for (Eigen::Index rc_idx = 0; rc_idx < json_data.real_coords.rows();
         ++rc_idx)
    {
        Eigen::RowVector3d const real_coords_single_point =
            json_data.real_coords.row(rc_idx);

        // Find cells for point
        auto const filtered_bulk_mesh =
            findCellsForPoints.find(real_coords_single_point, tolerance);

        // Check multiplicity
        int const actual_multiplicity = filtered_bulk_mesh->GetNumberOfCells();
        if (actual_multiplicity == 0)
        {
            OGS_FATAL(
                "Did not find any cell for the point {} (coordinates: {})",
                rc_idx, real_coords_single_point);
        }
        assert(actual_multiplicity > 0);
        if (actual_multiplicity != 1)
        {
            INFO(
                "Found more then one cell (namely {}) for point #{} "
                "(coordinates: {})",
                actual_multiplicity, rc_idx, real_coords_single_point);
        }

        // Convert to OGS
        std::unique_ptr<MeshLib::Mesh> ogs_mesh(
            MeshLib::VtkMeshConverter::convertUnstructuredGrid(
                filtered_bulk_mesh, "filtered_bulk_mesh"));

        auto const* bulk_element_ids = MeshLib::bulkElementIDs(*ogs_mesh);

        auto const& elements = ogs_mesh->getElements();

        // Compute natural coordinates
        // Found a solution, no need to test other elements ignoring
        // multiplicity.
        constexpr std::size_t elt_idx = 0;
        auto const& element = *elements[elt_idx];
        auto const bulk_element_id = bulk_element_ids->getComponent(elt_idx, 0);

        auto const& solver = SolverRegistry::getFor(element);
        auto const opt_natural_coords = solver.solve(
            element, real_coords_single_point, max_iter, real_coords_tolerance);

        result.append(opt_natural_coords, real_coords_single_point,
                      bulk_element_id, rc_idx);
    }

    return result.finished(
        json_data.initial_anchor_stress, json_data.maximum_anchor_stress,
        json_data.residual_anchor_stress, json_data.anchor_cross_sectional_area,
        json_data.anchor_stiffness);
}

template <typename T>
void addPointData(vtkUnstructuredGrid* grid, std::string const& name,
                  Eigen::MatrixX<T> const& data)
{
    INFO("converting {}", name);
    vtkNew<vtkAOSDataArrayTemplate<T>> array;
    array->SetName(name.c_str());
    array->SetNumberOfComponents(data.cols());
    array->SetNumberOfTuples(grid->GetNumberOfPoints());

    if (grid->GetNumberOfPoints() != data.rows())
    {
        OGS_FATAL(
            "Got {} rows in the table but expected {} rows, same as number of "
            "points in the grid.",
            data.rows(), grid->GetNumberOfPoints());
    }
    for (Eigen::Index i = 0; i < data.rows(); ++i)
    {
        // copy to contiguous storage
        Eigen::RowVectorX<T> const row = data.row(i);

        array->SetTypedTuple(i, row.data());
        DBUG("> {}", row);
    }

    grid->GetPointData()->AddArray(array);
}

template <typename T>
void addCellData(vtkUnstructuredGrid* grid, std::string const& name,
                 Eigen::MatrixX<T> const& data)
{
    INFO("converting {}", name);
    vtkNew<vtkAOSDataArrayTemplate<T>> array;
    array->SetName(name.c_str());
    array->SetNumberOfComponents(data.cols());
    array->SetNumberOfTuples(grid->GetNumberOfCells());

    if (grid->GetNumberOfCells() != data.rows())
    {
        OGS_FATAL(
            "Got {} rows in the table but expected {} rows, same as number of "
            "cells in the grid.",
            data.rows(), grid->GetNumberOfCells());
    }
    for (Eigen::Index i = 0; i < data.rows(); ++i)
    {
        // copy to contiguous storage
        Eigen::RowVectorX<T> const row = data.row(i);

        array->SetTypedTuple(i, row.data());
        DBUG("> {}", row);
    }

    grid->GetCellData()->AddArray(array);
}

/**
 * Creates a mesh from the points described by the passed \c result.
 *
 * The points 2n and 2n+1 will be connected by a line element that could
 * represent an anchor in an FEM model for instance.
 */
vtkSmartPointer<vtkUnstructuredGrid> toVTKGrid(
    ComputeNaturalCoordsResult const& result)
{
    INFO("converting points");
    auto const n_pts = result.natural_coords.rows();
    if (n_pts == 0)
    {
        OGS_FATAL(
            "No point was found in the mesh. Please check whether all anchors "
            "are within the model region.");
    }
    else if (n_pts == 1)
    {
        OGS_FATAL(
            "Only one point was found in the mesh. Please check whether all "
            "anchors are within the model region.");
    }
    else if (n_pts % 2 != 0)
    {
        // number should be even if multiplicity is ignored
        OGS_FATAL(
            "Number of points is not even. Anchors are believed to consist of "
            "a start and an end point.");
    }

    vtkNew<vtkUnstructuredGrid> grid;

    // Points ------------------------------------------------------------------
    vtkNew<vtkPoints> points;
    points->SetNumberOfPoints(n_pts);

    auto const& css = result.real_coords;
    for (Eigen::Index i = 0; i < css.rows(); ++i)
    {
        points->SetPoint(i, css(i, 0), css(i, 1), css(i, 2));

        if ((i % 2) == 1)
        {
            vtkNew<vtkLine> line;
            line->GetPointIds()->SetId(0, i - 1);
            line->GetPointIds()->SetId(1, i);
            grid->InsertNextCell(line->GetCellType(), line->GetPointIds());
        }
    }

    grid->SetPoints(points);

    // Point data --------------------------------------------------------------
    addPointData<double>(grid, "natural_coordinates", result.natural_coords);
    addPointData<std::size_t>(grid, "bulk_element_ids",
                              result.bulk_element_ids.cast<std::size_t>());
    addPointData<std::size_t>(grid, "point_cloud_node_ids",
                              result.point_cloud_node_ids.cast<std::size_t>());
    addCellData<double>(grid, "initial_anchor_stress",
                        result.initial_anchor_stress);
    addCellData<double>(grid, "maximum_anchor_stress",
                        result.maximum_anchor_stress);
    addCellData<double>(grid, "residual_anchor_stress",
                        result.residual_anchor_stress);
    addCellData<double>(grid, "anchor_cross_sectional_area",
                        result.anchor_cross_sectional_area);
    addCellData<double>(grid, "anchor_stiffness", result.anchor_stiffness);
    return grid;
}
}  // namespace ApplicationUtils
