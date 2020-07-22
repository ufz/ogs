+++
project = "Mechanics/Linear/PythonHertzContact/hertz_contact.prj"
author = "Christoph Lehmann"
date = "2018-08-06T11:41:00+02:00"
title = "Hertz Contact using Python Boundary Conditions"
weight = 4

[menu]
  [menu.benchmarks]
    parent = "python-bc"

+++

{{< data-link >}}

## Motivation of this test case

The aim of this test is:

* to show that it is possible to apply Dirichlet BCs at a boundary that changes over the course of the simulation
* to give an advanced use-case of the Python BC.
  Here essentially an iterative procedure for a contact problem is implemented
  within the Python BC.

## Problem description

Two elastic spheres of same radius $R$ are brought into contact.
The sphere centers are displaced towards each other by $w_0$, with increasing
values in every load step.
Due to symmetry reasons a flat circular contact area of radius $a$ forms.

{{< img src="../hertz-contact.png" >}}

The contact between the two spheres is modelled as a Dirichlet BC
on a varying boundary. The exact boundary and Dirichlet values for the
$y$ displacements are determined in a Python script.
Compared to the sketch above the prescribed $y$ displacements correspond
to $w_0/2$, because due to symmetry only half of the system (a section of the
lower sphere) is simulated.

## Analytical solution

The derivation of the formulae below can be found, e.g.,
in [this book (in German)](http://www.uni-magdeburg.de/ifme/l-festigkeit/pdf/Bertram-Gluege_Festkoerpermechanik2012.pdf).

The radius of the contact area is given by
$$
\begin{equation}
a = \sqrt{\frac{w_0 R}{2}}
\end{equation}
$$

The average pressure $\bar p$ over a the secant with distance $\xi$ to the
center of the contact area (cf. vertical dashed line in the sketch above) is assumed to be
$$
\begin{equation}
\bar p(\xi) = \kappa \sqrt{a^2 - \xi^2}
\label{eq:bar-p}
\end{equation}
$$
with the prefactor $\kappa$ given by
$$
\begin{equation}
\kappa = \frac{G}{R \cdot (1-\nu)}
\,.
\end{equation}
$$

The total force $F$ exerted on the contact area is given by
$$
\begin{equation}
F = \frac{8 a^3}{3\kappa}
\,.
\end{equation}
$$

## Results

Contact radii:

{{<img src="../contact_radii.png">}}

Average pressure $\bar{p}$:

{{<img src="../stress_at_contact.png">}}

Total force $F$:

{{<img src="../total_force.png">}}

The simulation results for contact radii and total force reproduce the
analytical square root and cubic laws, respectively.
For the average pressure $\bar p$ the analytical form of
Eq.&nbsp;$(\ref{eq:bar-p})$ is reproduced, too.
But for $\bar p$ there are rather strong deviations between the numerical
and analytical results, which might be due to deviations in the
contact radii&nbsp;$a$, due to unsufficient mesh resolution or due to
the chosen linear order of FEM ansatz functions.
