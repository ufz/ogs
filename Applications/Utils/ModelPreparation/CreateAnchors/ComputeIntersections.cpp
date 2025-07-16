/**
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ComputeIntersections.h"

#include <spdlog/spdlog.h>
#include <vtkCellData.h>
#include <vtkCellLocator.h>
#include <vtkExtractEdges.h>
#include <vtkIdTypeArray.h>
#include <vtkLine.h>

#include "BaseLib/Logging.h"

int GetGridDimension(vtkUnstructuredGrid* grid)
{
    vtkIdType numCells = grid->GetNumberOfCells();
    int maxDim = 0;

    vtkSmartPointer<vtkGenericCell> cell =
        vtkSmartPointer<vtkGenericCell>::New();

    for (vtkIdType i = 0; i < numCells; ++i)
    {
        grid->GetCell(i, cell);
        int dim = cell->GetCellDimension();
        if (dim > maxDim)
        {
            maxDim = dim;
        }
    }
    return maxDim;
}

bool isPointClose(const Eigen::Vector3d& a, const Eigen::Vector3d& b,
                  double tol)
{
    return (a - b).squaredNorm() <= tol * tol;
}

// Comparator to sort intersections by t value
bool compareByT(const IntersectionResult& a, const IntersectionResult& b)
{
    return a.t < b.t;
}

std::vector<IntersectionResult> findOrderedIntersections(
    vtkUnstructuredGrid* grid,
    Eigen::Vector3d const& p0,
    Eigen::Vector3d const& p1,
    const double free_fraction,
    double const tol)
{
    std::vector<IntersectionResult> intersections;

    intersections.push_back(IntersectionResult(p0, 0.0));
    intersections.push_back(IntersectionResult(p1, 1.0));
    Eigen::Vector3d transition_point = p0 + (p1 - p0) * free_fraction;
    if (!isPointClose(p0, transition_point, tol) &&
        !isPointClose(p1, transition_point, tol))
    {
        INFO("Adding transition point at t = {:.2f} for free fraction {:.2f}",
             free_fraction, free_fraction);
        intersections.push_back(
            IntersectionResult(transition_point, free_fraction));
    }
    // For cell intersection
    vtkSmartPointer<vtkGenericCell> cell =
        vtkSmartPointer<vtkGenericCell>::New();

    vtkNew<vtkIdList> filteredCellIds;
    vtkNew<vtkCellLocator> locator;
    vtkSmartPointer<vtkDataSet> cells_to_be_intersected;
    if (GetGridDimension(grid) < 3)
    {
        vtkSmartPointer<vtkExtractEdges> extractEdges =
            vtkSmartPointer<vtkExtractEdges>::New();
        extractEdges->SetInputData(grid);
        extractEdges->Update();

        cells_to_be_intersected = extractEdges->GetOutput();
    }
    else
    {
        cells_to_be_intersected = grid;
    }
    locator->SetDataSet(cells_to_be_intersected);
    locator->BuildLocator();
    locator->FindCellsAlongLine(p0.data(), p1.data(), tol, filteredCellIds);
    for (vtkIdType i = 0; i < filteredCellIds->GetNumberOfIds(); ++i)
    {
        vtkIdType cellId = filteredCellIds->GetId(i);
        cells_to_be_intersected->GetCell(cellId, cell);

        double t;
        double x[3];
        double pcoords[3];
        int subId;

        int intersected = cell->IntersectWithLine(p0.data(), p1.data(), tol, t,
                                                  x, pcoords, subId);

        Eigen::Vector3d x_vec(x[0], x[1], x[2]);

        if (!intersected)
        {
            continue;
        }
        bool isDuplicate = false;
        for (const auto& existing : intersections)
        {
            if (isPointClose(existing.point, x_vec, tol))
            {
                isDuplicate = true;
                break;
            }
        }
        if (!isDuplicate && (t < tol || t > free_fraction))
        {
            IntersectionResult result;
            // result.cellId = cellId;
            result.t = t;
            result.point = x_vec;
            intersections.push_back(result);
        }
    }
    std::sort(intersections.begin(), intersections.end(), compareByT);
    return intersections;
}

std::vector<std::vector<IntersectionResult>> getOrderedAnchorCoords(
    vtkUnstructuredGrid* grid,
    Eigen::MatrixX3d const& realcoords,
    Eigen::VectorXd const& free_fraction,
    double const tol)
{
    std::vector<std::vector<IntersectionResult>> ordered_intersections(
        realcoords.rows() / 2);
    assert(realcoords.rows() % 2 == 0);
    assert(free_fraction.size() == realcoords.rows() / 2);
    for (int i = 0; i < realcoords.rows(); i += 2)
    {
        ordered_intersections[i / 2] = findOrderedIntersections(
            grid, realcoords.row(i), realcoords.row(i + 1),
            free_fraction(i / 2), tol);
    }
    return ordered_intersections;
}

AU::ComputeNaturalCoordsResult setPhysicalPropertiesForIntersectionPoints(
    std::vector<std::vector<IntersectionResult>> const& anchor_coords,
    AU::ComputeNaturalCoordsResult const& original_anchor_data)
{
    std::vector<IntersectionResult> all_intersections;
    std::vector<int> anchorids;
    for (auto&& [anchor_id, intersections] :
         ranges::views::enumerate(anchor_coords))
    {
        for (auto&& [index, result] : ranges::views::enumerate(intersections))
        {
            if (index > 0 && index < intersections.size() - 1)
            {
                anchorids.push_back(anchor_id);
                all_intersections.push_back(result);
            }
            anchorids.push_back(anchor_id);
            all_intersections.push_back(result);
        }
    }
    assert(all_intersections.size() == anchorids.size());
    assert(all_intersections.size() % 2 == 0);
    Eigen::MatrixXd realcoords(anchorids.size(), 3);
    auto const number_of_anchors = anchorids.size() / 2;
    Eigen::VectorXd initial_stress(number_of_anchors);
    Eigen::VectorXd maximum_stress(number_of_anchors);
    Eigen::VectorXd residual_stress(number_of_anchors);
    Eigen::VectorXd cross_sectional_area(number_of_anchors);
    Eigen::VectorXd stiffness(number_of_anchors);
    for (std::size_t i = 0; i < all_intersections.size(); ++i)
    {
        realcoords.row(i).noalias() = all_intersections[i].point;
        if (i % 2 == 0)
        {
            initial_stress(i / 2) =
                original_anchor_data.initial_anchor_stress(anchorids[i]);
            maximum_stress(i / 2) =
                original_anchor_data.maximum_anchor_stress(anchorids[i]);
            residual_stress(i / 2) =
                original_anchor_data.residual_anchor_stress(anchorids[i]);
            cross_sectional_area(i / 2) =
                original_anchor_data.anchor_cross_sectional_area(anchorids[i]);
            stiffness(i / 2) =
                original_anchor_data.anchor_stiffness(anchorids[i]);
        }
    }
    AU::ComputeNaturalCoordsResult anchor_data;
    anchor_data.real_coords = realcoords;
    anchor_data.initial_anchor_stress = initial_stress;
    anchor_data.maximum_anchor_stress = maximum_stress;
    anchor_data.residual_anchor_stress = residual_stress;
    anchor_data.anchor_cross_sectional_area = cross_sectional_area;
    anchor_data.anchor_stiffness = stiffness;
    return anchor_data;
}
