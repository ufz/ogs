+++
date = "2019-08-02T11:33:45+01:00"
title = "Liquid flow with time dependent boundary conditions and source term"
weight = 171
project = "Parabolic/LiquidFlow/TimeDependentHeterogeneousSourceTerm/TimeDependentHeterogeneousSourceTerm.prj"
author = "Thomas Fischer"

[menu]
  [menu.benchmarks]
    parent = "liquid-flow"

+++

{{< data-link >}}

## Motivation

In real world examples the boundary conditions or source terms can vary over time
and can be heterogeneous in space. This behaviour can be modelled using the
TimeDependentHeterogeneousParameter for boundary conditions or source terms.

## Specification in OGS project file

In the parameter specification section of the project file it is possible to add
a parameter type with the type `TimedependentHeterogeneousParameter`.

```xml
<parameter>
    <name>ParameterForSourceTerm</name>
    <type>TimeDependentHeterogeneousParameter</type>
    <time_series>
        <pair>
            <time>0</time>
            <parameter_name>parameter_for_timestep1</parameter_name>
        </pair>
        <pair>
            <time>1</time>
            <parameter_name>parameter_for_timestep2</parameter_name>
        </pair>
        ...
        <pair>
            <time>end_time</time>
            <parameter_name>parameter_for_end_time</parameter_name>
        </pair>
    </time_series>
</parameter>
```

Of course, the referenced parameters for the particular time steps have to be
defined also. Values of the parameter are piecewise linear interpolated.

## Example

This simple example should demonstrate the use of the time dependent
heterogeneous parameter. We start with homogeneous parabolic problem:
$$
\begin{equation}
s\\;\frac{\partial p}{\partial t} + k\; \Delta p = q(t,x) \quad \text{in }\Omega
\end{equation}
$$
w.r.t boundary conditions
$$
\eqalign{
p(t, x) = g_D(t, x) &\quad \text{on }\Gamma_D,\cr
k\\;{\partial p(x) \over \partial n} = g_N(x) &\quad \text{on }\Gamma_N,
}$$

The example the domain $\Omega = [0,1]^2$ is a square. On the left
($x=0$) side and the right ($x=1$) side time dependent Dirichlet-type boundary
conditions are set. Until half of the simulation time high pressure values are
set on the left side and low pressure values on the right side. In the second
half of the simulation there are low pressure values on the left side and high
pressure values on the right side. Additionally, the source term $q$ acts in
the first quarter as a source, in the second quarter as a sink, in the third
quarter as a source, and in the last quarter as a sink again.

## Results

<video width="838" height="762" controls>
<source src="TimeDependentHeterogeneousBoundaryConditionsAndSourceTerm.mp4" type="video/mp4" />
</video>
<p>
<strong>Download Video:</strong>
<a href="TimeDependentHeterogeneousBoundaryConditionsAndSourceTerm.mp4">"MP4"</a>
</p>
