This source term models anchors between two points in the bulk mesh.

Each one-dimensional line element in the `EmbeddedAnchor` source term mesh defines an anchor. It contains `bulk_element_ids`, `natural coordinates` as point data. And it contains`anchor_stiffness` and `anchor_cross_sectional_area`, `initial_anchor_stress`, `maximum_anchor_stress`, and `residual_anchor_stress` as cell data.
Stiffness and stresses are given in pressure units.
