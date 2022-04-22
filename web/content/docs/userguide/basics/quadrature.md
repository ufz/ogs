+++
date = "2022-03-02T09:30:00+01:00"
title = "Quadrature schemes and extrapolation"
author = "Christoph Lehmann"
weight = 7


[menu]
  [menu.userguide]
    parent = "basics"
+++

OpenGeoSys has implemented certain numerical integration/quadrature rules for
each mesh element type. Currently, all of the implemented quadrature rules are
Gauss-Legendre quadrature or schemes in that spirit. I.e., integration points
are located inside the elements, not at its faces, edges or corners.

The table below gives an overview over the implemented schemes and which
polynomials can be integrated exactly with the respective schemes. The data are
[unit tested](https://gitlab.opengeosys.org/ogs/ogs/-/tree/master/Tests/MathLib/TestGaussLegendreIntegration.cpp).

Integration methods are implemented up to *integration order* 4.
For each integration order *n* there is a maximum polynomial degree *P* that can be
integrated exactly with the respective integration method. For the classical
Gauss-Legendre integration the following holds: *P = 2 · n - 1*. Methods that
fulfill this relation are marked **bold** in column **P** in the table, methods
that are deficient are marked *italic*.

The columns contain the following data:

* #IP: The number of integration points of the integration method.
* P: The maximum polynomial degree that the integration method can integrate
  exactly. The monomials in the polynomial are *x*<sup>i</sup> *y*<sup>j</sup>
  *z*<sup>k</sup> with *i + j + k ≤ P*
* Q: The maximum polynomial degree that the integration method can integrate
  exactly. The monomials in the polynomial are *x*<sup>i</sup> *y*<sup>j</sup>
  *z*<sup>k</sup> with *i ≤ Q*, *j ≤ Q* and *k ≤ Q*. I.e., maximum monomial
  degrees are higher in the Q column than in the P column.


| Integration order → | 1       | 1     | 1     | 2       | 2     | 2     | 3       | 3     | 3     | 4       | 4     | 4     |
| ------------------- | ------: | ----: | ----: | ------: | ----: | ----: | ------: | ----: | ----: | ------: | ----: | ----: |
| **Mesh element ↓**  | **#IP** | **Q** | **P** | **#IP** | **Q** | **P** | **#IP** | **Q** | **P** | **#IP** | **Q** | **P** |
| Point               | 1       | \*    | \*    | 1       | \*    | \*    | 1       | \*    | \*    | 1       | \*    | \*    |
| Line                | 1       | 1     | **1** | 2       | 3     | **3** | 3       | 5     | **5** | 4       | 7     | **7** |
| Quad                | 1       | 1     | **1** | 4       | 3     | **3** | 9       | 5     | **5** | 16      | 7     | **7** |
| Hexahedron          | 1       | 1     | **1** | 8       | 3     | **3** | 27      | 5     | **5** | 64      | 7     | **7** |
| Triangle            | 1       | 0     | **1** | 3       | 1     | **3** | 4       | 1     | *3*   | 7       | 2     | *5*   |
| Tetrahedron         | 1       | 0     | **1** | 5       | 1     | **3** | 14      | 1     | **5** | 20      | 1     | *5*   |
| 3-sided prism       | 1       | 0     | **1** | 6       | 1     | **3** | 21      | 2     | *3*   | 28      | 2     | *5*   |
| Pyramid             | 1       | 1     | **1** | 5       | 1     | **3** | 13      | 3     | *3*   | 13      | 3     | *3*   |

Note, that on pyramids det(*J*) varies over the mesh element, even for linear
elements. Therefore, for pyramids we are actually integrating a polynomial of
higher degree than the degrees P and Q given in the table above.

## Extrapolation of integration point data to mesh nodes

OpenGeoSys can extrapolate integration point data to mesh nodes for easy
postprocessing. **Note, however, that the extrapolation procedures is not exact
and can lead to more or less subtle errors that are hard to find!**

Since OpenGeoSys extrapolates integration point data element-wise, the number of
integration points of each mesh element must be greater or equal to the number of nodes of the element.
Therefore, the ability to extrapolate data is linked to the chosen integration
order, i.e., the number of integration points must be greater or equal to the
number of nodes. The relation is presented in the table below.

The columns contain the following data:

* #IP: The number of integration points of the integration method.
* L:
  Whether extrapolation can be performed on the *linear* elements (e.g. Quad4) of the
  respective type with the given integration order.
* Q:
  Whether extrapolation can be performed on the *quadratic* elements (e.g.
  Hex20) of the respective type with the given integration order.

As you can see, for integration order ≥ 2 extrapolation works for all linear
elements and for integration order ≥ 3 extrapolation works for all element types
implemented in OpenGeoSys.

| Integration order → | 1       | 1     | 1     | 2       | 2     | 2     | 3       | 3     | 3     | 4       | 4     | 4     |
| ------------------- | ------: | :---: | :---: | ------: | :---: | :---: | ------: | :---: | :---: | ------: | :---: | :---: |
| **Mesh element ↓**  | **#IP** | **L** | **Q** | **#IP** | **L** | **Q** | **#IP** | **L** | **Q** | **#IP** | **L** | **Q** |
| Point               | 1       |   ✓   |   ✓   | 1       |   ✓   |   ✓   | 1       |   ✓   |   ✓   | 1       |   ✓   |   ✓   |
| Line                | 1       |   –   |   –   | 2       |   ✓   |   –   | 3       |   ✓   |   ✓   | 4       |   ✓   |   ✓   |
| Quad                | 1       |   –   |   –   | 4       |   ✓   |   –   | 9       |   ✓   |   ✓   | 16      |   ✓   |   ✓   |
| Hexahedron          | 1       |   –   |   –   | 8       |   ✓   |   –   | 27      |   ✓   |   ✓   | 64      |   ✓   |   ✓   |
| Triangle            | 1       |   –   |   –   | 3       |   ✓   |   –   | 4       |   ✓   |   ✓   | 7       |   ✓   |   ✓   |
| Tetrahedron         | 1       |   –   |   –   | 5       |   ✓   |   –   | 14      |   ✓   |   ✓   | 20      |   ✓   |   ✓   |
| 3–sided prism       | 1       |   –   |   –   | 6       |   ✓   |   –   | 21      |   ✓   |   ✓   | 28      |   ✓   |   ✓   |
| Pyramid             | 1       |   –   |   –   | 5       |   ✓   |   –   | 13      |   ✓   |   ✓   | 13      |   ✓   |   ✓   |

