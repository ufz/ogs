/**
 * \file
 * \copyright
 * Copyright (c) 2012-2022, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

#include "ReferenceElementUtils.h"

namespace ReferenceElementUtils
{
// Returns the coordinates as a span of dynamic size.
BaseLib::DynamicSpan<const std::array<double, 3>>
getNodeCoordsOfReferenceElement(MeshLib::CellType const cell_type)
{
    using namespace MeshLib;

    switch (cell_type)
    {
        case CellType::POINT1:
            return getNodeCoordsOfReferenceElement<Point>();
        case CellType::LINE2:
            return getNodeCoordsOfReferenceElement<Line>();
        case CellType::LINE3:
            return getNodeCoordsOfReferenceElement<Line3>();
        case CellType::TRI3:
            return getNodeCoordsOfReferenceElement<Tri>();
        case CellType::TRI6:
            return getNodeCoordsOfReferenceElement<Tri6>();
        case CellType::QUAD4:
            return getNodeCoordsOfReferenceElement<Quad>();
        case CellType::QUAD8:
            return getNodeCoordsOfReferenceElement<Quad8>();
        case CellType::QUAD9:
            return getNodeCoordsOfReferenceElement<Quad9>();
        case CellType::TET4:
            return getNodeCoordsOfReferenceElement<Tet>();
        case CellType::TET10:
            return getNodeCoordsOfReferenceElement<Tet10>();
        case CellType::HEX8:
            return getNodeCoordsOfReferenceElement<Hex>();
        case CellType::HEX20:
            return getNodeCoordsOfReferenceElement<Hex20>();
        case CellType::PRISM6:
            return getNodeCoordsOfReferenceElement<Prism>();
        case CellType::PRISM15:
            return getNodeCoordsOfReferenceElement<Prism15>();
        case CellType::PYRAMID5:
            return getNodeCoordsOfReferenceElement<Pyramid>();
        case CellType::PYRAMID13:
            return getNodeCoordsOfReferenceElement<Pyramid13>();
        default:
            OGS_FATAL("Unsupported cell type " + CellType2String(cell_type));
    }
}

std::shared_ptr<MeshLib::Element const> getReferenceElement(
    MeshLib::CellType const cell_type)
{
    auto const create = []<typename ElementType>(Type<ElementType>)
    {
        auto ref_elt = std::make_shared<ReferenceElement<ElementType>>();
        return std::shared_ptr<MeshLib::Element const>{ref_elt,
                                                       &ref_elt->element};
    };

    using namespace MeshLib;

    switch (cell_type)
    {
        case CellType::POINT1:
            return create(Type<Point>{});
        case CellType::LINE2:
            return create(Type<Line>{});
        case CellType::LINE3:
            return create(Type<Line3>{});
        case CellType::TRI3:
            return create(Type<Tri>{});
        case CellType::TRI6:
            return create(Type<Tri6>{});
        case CellType::QUAD4:
            return create(Type<Quad>{});
        case CellType::QUAD8:
            return create(Type<Quad8>{});
        case CellType::QUAD9:
            return create(Type<Quad9>{});
        case CellType::TET4:
            return create(Type<Tet>{});
        case CellType::TET10:
            return create(Type<Tet10>{});
        case CellType::HEX8:
            return create(Type<Hex>{});
        case CellType::HEX20:
            return create(Type<Hex20>{});
        case CellType::PRISM6:
            return create(Type<Prism>{});
        case CellType::PRISM15:
            return create(Type<Prism15>{});
        case CellType::PYRAMID5:
            return create(Type<Pyramid>{});
        case CellType::PYRAMID13:
            return create(Type<Pyramid13>{});
        default:
            OGS_FATAL("Unsupported cell type " + CellType2String(cell_type));
    }
}

std::vector<std::array<double, 3>> getCoordsInReferenceElementForTest(
    MeshLib::Element const& element)
{
    using Result = std::vector<std::array<double, 3>>;

    auto const dim = element.getDimension();

    if (dim == 0)
    {
        return Result{{0, 0, 0}};
    }

    auto const ref_elt_ptr =
        ReferenceElementUtils::getReferenceElement(element.getCellType());

    const std::array coords_1d{1. / 11, 0.1,    1. / 9, 1. / 7, 0.2, 1. / 3,
                               0.5,     2. / 3, 6. / 7, 8. / 9, 0.9, 10. / 11};
    auto const num_coords_1d = coords_1d.size();

    if (dim == 1)
    {
        Result coords;
        coords.reserve(2 * num_coords_1d);
        for (std::size_t i = 0; i < num_coords_1d; ++i)
        {
            for (auto const pm_i : {-1, 1})
            {
                std::array const cs{pm_i * coords_1d[i], 0., 0.};
                if (ref_elt_ptr->isPntInElement(MathLib::Point3d{cs}))
                {
                    coords.push_back(cs);
                }
            }
        }

        EXPECT_LE(num_coords_1d, coords.size())
            << "For all element types at least the positive natural "
               "coordinates should be contained in the result.";
        return coords;
    }

    if (dim == 2)
    {
        Result coords;
        coords.reserve(4 * num_coords_1d * num_coords_1d);
        for (std::size_t i = 0; i < num_coords_1d; ++i)
        {
            for (auto const pm_i : {-1, 1})
            {
                for (std::size_t j = 0; j < num_coords_1d; ++j)
                {
                    for (auto const pm_j : {-1., 1.})
                    {
                        std::array const cs{pm_i * coords_1d[i],
                                            pm_j * coords_1d[j], 0.};
                        if (ref_elt_ptr->isPntInElement(MathLib::Point3d{cs}))
                        {
                            coords.push_back(cs);
                        }
                    }
                }
            }
        }

        auto const rough_correction_for_triangles = 2;
        EXPECT_LE(
            num_coords_1d * num_coords_1d / rough_correction_for_triangles,
            coords.size())
            << "For all element types at least the positive natural "
               "coordinates that sum to <= 1 should be contained in the "
               "result.";
        return coords;
    }

    if (dim == 3)
    {
        Result coords;
        coords.reserve(8 * num_coords_1d * num_coords_1d * num_coords_1d);

        for (std::size_t i = 0; i < num_coords_1d; ++i)
        {
            for (auto const pm_i : {-1, 1})
            {
                for (std::size_t j = 0; j < num_coords_1d; ++j)
                {
                    for (auto const pm_j : {-1., 1.})
                    {
                        for (std::size_t k = 0; k < num_coords_1d; ++k)
                        {
                            for (auto const pm_k : {-1., 1.})
                            {
                                std::array const cs{pm_i * coords_1d[i],
                                                    pm_j * coords_1d[j],
                                                    pm_k * coords_1d[k]};
                                if (ref_elt_ptr->isPntInElement(
                                        MathLib::Point3d{cs}))
                                {
                                    coords.push_back(cs);
                                }
                            }
                        }
                    }
                }
            }
        }

        auto const rough_correction_for_tets = 5;
        EXPECT_LE(num_coords_1d * num_coords_1d * num_coords_1d /
                      rough_correction_for_tets,
                  coords.size())
            << "For all element types at least the positive natural "
               "coordinates that sum to <= 1 should be contained in the "
               "result.";
        return coords;
    }

    OGS_FATAL("Unsupported element dimension " + std::to_string(dim));
}
}  // namespace ReferenceElementUtils
