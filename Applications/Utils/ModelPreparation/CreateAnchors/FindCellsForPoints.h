/**
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#pragma once

#include <vtkCellData.h>
#include <vtkCellTreeLocator.h>
#include <vtkCharArray.h>
#include <vtkExtractCells.h>
#if VTK_VERSION_STRIPPED < 940
#include <vtkIdFilter.h>
#else
#include <vtkGenerateIds.h>
#endif
#include <vtkIdList.h>
#include <vtkNew.h>
#include <vtkUnstructuredGrid.h>

#include <Eigen/Core>
#include <array>
#include <range/v3/view/iota.hpp>
#include <range/v3/view/transform.hpp>
#include <vector>

#include "BaseLib/Logging.h"

namespace ApplicationUtils
{
/**
 * An algorithm finding mesh elements within a certain tolerance from a given
 * point.
 *
 * The algorithm will find a single element if for points lying within a mesh
 * element (unless the tolerance is very high), and will find multiple elements
 * for points on element faces/edges/corners.
 */
class FindCellsForPoints
{
public:
    void initialize(vtkUnstructuredGrid* const bulk_mesh)
    {
#if VTK_VERSION_STRIPPED < 940
        vtkNew<vtkIdFilter> id_filter;
        id_filter->SetInputData(bulk_mesh);
        id_filter->SetCellIdsArrayName("bulk_element_ids");
        id_filter->SetPointIds(false);
        id_filter->SetCellIds(true);
        id_filter->Update();
        bulk_mesh_ = id_filter->GetOutput();
#else
        vtkNew<vtkGenerateIds> generate_ids;
        generate_ids->SetInputData(bulk_mesh);
        generate_ids->SetCellIdsArrayName("bulk_element_ids");
        generate_ids->PointIdsOff();
        generate_ids->CellIdsOn();
        generate_ids->Update();
        bulk_mesh_ =
            static_cast<vtkUnstructuredGrid*>(generate_ids->GetOutput());
#endif
        locator_ = vtkSmartPointer<vtkCellTreeLocator>::New();
        locator_->SetDataSet(bulk_mesh_);
        locator_->BuildLocator();
    }

    std::vector<vtkIdType> findCellIds(Eigen::Vector3d const& coords,
                                       double const tolerance)
    {
        using namespace ranges;
        auto const cell_ids = findCellIdsImpl(coords, tolerance);
        return views::iota(0, cell_ids->GetNumberOfIds()) |
               views::transform([&](vtkIdType i)
                                { return cell_ids->GetId(i); }) |
               to_vector;
    }

    vtkSmartPointer<vtkUnstructuredGrid> find(Eigen::Vector3d const& coords,
                                              double const tolerance)
    {
        auto const cell_ids = findCellIdsImpl(coords, tolerance);

        vtkNew<vtkExtractCells> extract_cells;
        extract_cells->SetInputData(bulk_mesh_);
        extract_cells->SetCellList(cell_ids);
        extract_cells->Update();
        auto const filtered_bulk_mesh = extract_cells->GetOutput();

        return filtered_bulk_mesh;
    }

private:
    vtkSmartPointer<vtkIdList> findCellIdsImpl(Eigen::Vector3d const& coords,
                                               double const tolerance)
    {
        std::array<double, 6> bbox = {
            coords[0] - tolerance, coords[0] + tolerance,
            coords[1] - tolerance, coords[1] + tolerance,
            coords[2] - tolerance, coords[2] + tolerance};

        vtkNew<vtkIdList> cell_ids;
        locator_->FindCellsWithinBounds(bbox.data(), cell_ids);

        return cell_ids;
    }

    vtkSmartPointer<vtkDataSet> bulk_mesh_;
    vtkSmartPointer<vtkCellTreeLocator> locator_;
};

}  // namespace ApplicationUtils
