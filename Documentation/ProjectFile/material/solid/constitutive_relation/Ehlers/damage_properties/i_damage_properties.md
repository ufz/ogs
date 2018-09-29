Damage properties for Ehlers with local isotropic damage. This is an optional
subtree; when not set the "usual" Ehlers material model will be used without
damage.

For further information, refer to the implementation manual and to the following
publication: Parisio, F and Laloui, L (2017). Plastic-damage modeling of
saturated quasi-brittle shales. Int J Rock Mech Min Sci, 93:295-306
\cite Parisio2017

\attention This implementation does not contain strain regularization. The model
should be used with care as pathological mesh-dependency may arise. To overcome
this issue, the user can employ the non-local integral formulation.
