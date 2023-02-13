+++
date = "2018-02-27T11:00:13+01:00"
title = "Time loop"
author = "Feliks Kiszkurno"
weight = 6
+++
<div class="note">

### Work in progress

This page is a work in progress.

It was published in this state to make existing content available to users and highlight missing parts to contributors.

**Contributors:** please see Documentation Contribution Guide to contribute to the documentation.

**Users:** the content of this page has been verified and is correct. Please return later for more content!

</div>

This block defines all aspects related to the progress of time in the simulation and controls when the output files will be written, their naming pattern and content.

## Process

First section of this block, the one controlling time progress has to be set up separately for each process defined in [processes block](/docs/userguide/blocks/processes/#process-variables) separately.
There can be more than one process defined.
They are distinguished with the parameter "ref" in the `<process>` tag:

```xml
<process ref="Process_name">
    ...
</process>
```

### Time stepping

Numerical simulations proceed in time steps.
Their size is important factor impacting if a stable solution will be found.

In OpenGeoSys the following time step definitions are available:

- FixedTimeStepping
- SingleStep
- IterationNumberBasedTimeStepping
- EvolutionaryPIDcontroller

#### Single time step

This is the simplest type of time steppings used for solution of elliptic pdes as in the `SteadyStateDiffusion` process.
It only requires its type to be given:

```xml
<time_stepping>
    <type>SingleStep</type>
</time_stepping>
```

#### Fixed time stepping

This type of time stepping will run the simulation with predetermined number of steps from start to end of the simulation.
The steps are defined with parameter pair of number of repetitions and length. It is allowed to define multiple pairs. This way, the time step can vary throughout the experiment.

In the example below, a duration of $1$ time unit is simulated, with first 5 time steps being $0.1$ time unit long and next 10 being $0.05$ time unit long.

```xml
<time_stepping>
    <type>FixedTimeStepping</type>
    <t_initial>0</t_initial>
    <t_end>1</t_end>
    <timesteps>
        <pair>
            <repeat> 5 </repeat>
            <delta_t> 0.1 </delta_t>
        </pair>
        <pair>
            <repeat> 10 </repeat>
            <delta_t> 0.05 </delta_t>
        </pair>
    </timesteps>
</time_stepping>
```

An arbitrary number of repeat-delta_t pairs can be provided (though at least one has to be defined).
It is not required that sum of duration of all time steps perfectly matches the value of end time.
If this is not the case, a time step at t_end will be added.
The first time step of the simulation is always at t_initial.

The unit for time has to be consistent with the rest of the units used in the experiment.
Following SI system, second is a common choice.
This is how to define experiment simulating a year of time in daily time step using second as time unit:

```xml
<time_stepping>
    <type>FixedTimeStepping</type>
    <t_initial>0</t_initial>
    <!--t_end>365*24*60*60</t_end-->
    <t_end>31536000</t_end>
    <timesteps>
        <pair>
            <repeat>365</repeat>
            <!--delta_t>24*60*60</delta_t-->
            <delta_t>86400</delta_t>
        </pair>
    </timesteps>
</time_stepping>
```

The fixed time stepping is the most commonly used type.
As it is one of the simplest available, it makes a good starting point.

#### Iteration number based time stepping

[See doxygen](https://doxygen.opengeosys.org/df/df4/ogs_file_param__prj__time_loop__processes__process__time_stepping__iterationnumberbasedtimestepping)
The examples discussed above, would be defined in the project file as follows:

```xml
<time_stepping>
    <type>IterationNumberBasedTimeStepping</type>
    <t_initial>0</t_initial>
    <t_end>1e16</t_end>
    <initial_dt>1</initial_dt>
    <minimum_dt>1</minimum_dt>
    <maximum_dt>9</maximum_dt>
    <number_iterations>2 6 8 9</number_iterations>
    <multiplier>1.6 1.0 0.5 0.25</multiplier>
</time_stepping>
```

#### Evolutionary PID controller

[See doxygen](https://doxygen.opengeosys.org/d3/d86/ogs_file_param__prj__time_loop__processes__process__time_stepping__evolutionarypidcontroller)

### Error tolerances

Error tolerances will be applied to the solution vector
There are two two ways of defining error tolerances:

- relative `<reltosl> </reltols>`
- absolute `<abstols> </abstols>`

Both of them can be defined as single value, that will be applied to all process variables, or with multiple ones applied to them individually.

If tolerances per process variable are provided, the order of values defined inside of the tags `<abstols> </abstols>` and `<reltosl> </reltols>` has to match order of process variables defined in [processes](/docs/userguide/blocks/processes/).
If process variable is directional, XYZ pattern is followed.
For example for a 3D THM problem with directional displacement $\mathbf{u}$ following order has to be used: $T$, $p$, $u_x$, $u_y$, $u_z$.
The order of $T$, $p$ and $u$ is prescribed by the THM process.

Depending on process and simulation setup, the number of variables in the solution vector can vary.
In the following example, there are four tolerances given in "abstol": one for T, one for p and two for $\mathbf{u}$ ($u_x$ and $u_y$ as this is a 2D problem):

```xml
<process_variables>
    <displacement>displacement</displacement>
    <pressure>pressure</pressure>
    <temperature>temperature</temperature>
</process_variables>
```

The same order is preserved for [output variables](/docs/userguide/blocks/time_loop/#output).

### Convergence criteria

In this part convergence criterion has to be selected. Following options are available:

- [DeltaX and PerComponentDeltaX](/docs/userguide/blocks/time_loop/#deltax-and-percomponentdeltax)
- [Residual and PerComponentResidual](/docs/userguide/blocks/time_loop/#residual-and-percomponentresidual)

All of the criterions mentioned above compare value quantifying error with user defined tolerances.

There are three way how error tolerances can be set up:

- only relative tolerances $e_{\mathrm{rel}} \le \epsilon_{\mathrm{rel}}$
- only absolute tolerances $e_{\mathrm{abs}} \le \epsilon_{\mathrm{abs}}$
- relative and absolute tolerances $e_{\mathrm{abs}} \le \epsilon_{\mathrm{abs}} \vee e_{\mathrm{rel}} \le \epsilon_{\mathrm{rel}} $

where $e$ indicates error and $\epsilon$ indicates tolerance.

In the last case convergence criterion will be satisfied with **at least one** of the tolerances being met.

Both relative and absolute tolerances can be defined as single value or with separate values for each component, see [error tolerances](/docs/userguide/blocks/time_loop/#error-tolerances) for more details.

#### DeltaX and PerComponentDeltaX

DeltaX criteria are based on the norm value of the update (difference) of the solution between two consecutive iterations and is applied to each process variable.
If process variable is vectorial, criteria will be applied to each component of the vector.

In this criteria, the error is defined as difference between results vector obtained from two consecutive iterations.
If this difference drops below tolerance, the criteria is fulfilled and solver will proceed to the next time step.

The DeltaX criterion usually uses the [euclidean norm](/docs/userguide/blocks/time_loop/#euclidan-norm) as indicated by the `<norm_type>` tag value `NORM2`.
Other possible values are `NORM1` and `INFINITY_N`.  

The norm value of the solution vector is compared with the error tolerances.
They can be provided as absolute, relative or both.

#### Residual and PerComponentResidual

Differently than [DeltaX](/docs/userguide/blocks/time_loop/#deltax-and-percomponentdeltax), the residual based convergence criteria use the residuum values for evaluation, not the increment values.

The residual value is the difference between two sides of modeled equation.
In analytical solution, they would be perfectly equal, but in numerical solution, this will usually not be the case.
The lower the residuals are, the closer numerical solution is to fulfilling the equation that is being solved.

In this criteria, error is defined as imbalance between left and right hand sides of the equation.
Once the residual values drop below user defined tolerances, solution is accepted and solver can proceed to the next time step.

The normal of the residuum vector is compared with the error tolerances.
They can be provided as absolute, relative or both.

### Time discretizationn

TODO: describe different options for time discretization

For the time being only backward Euler time stepping scheme is available.

## Output

In this section all parameters controlling the output of the data from the simulation are defined.

Currently OpenGeoSys supports two output file formats: VTK and XDMF.
The desired format for output files can be declared as follows:

```xml
<type>VTK</type>
```

The number of time steps at which the output files are written doesn't need to match the simulation time steps.
They can be defined in the same way as in [fixed time step time stepping](/docs/userguide/blocks/time_loop/#fixed-time-stepping).

If specific times of the simulation are of a special interest to the user, OpenGeoSys can write the output files at times provided in a list in `<fixed_output_times> </fixed_output_times>` tag:

```xml
<fixed_output_times>
    0
    10
    20
    40
    80
    160
    320  
</fixed_output_times>
```

The naming of output files written at different time steps is handled automatically by OpenGeoSys.
User can use two tags that allow control over how the file names are structured and what they contain: tags `<suffix> </suffix>` and `<prefix> </prefix>`.
The rules regarding their content are the same.
They should contain text string.
Variables can be used to allow those strings to vary between files.
Following variables can be called:

- meshname
- timestep
- time
- iteration

They can be used in the file name with syntax illustrated by following example:

```xml
<suffix>Experiment_name</suffix>
<prefix>_at_time_step_{:timestep}_time_{:time}_iteration_{:iteration}</prefix>
```

Name of each file will be automatically appended by a file extension appropriate to the content of `<type> </type>` tag.
No user action is needed in this regard.

Block `<variables> </variables>` allows user to control what process variables are included in the output files.
If only some of them are of interest, it may be convenient to reduce size of the output files by only writing the necessary ones.
If only temperature and pressure are important, this block can be defined as follows:

```xml
<variables>
    <variable>temperature</variable>
    <variable>pressure</variable>
</variables>
```

The list of available variables differs between processes.
TODO: Create list of variables available in different processes

In order to save space, the compression of output files can be enabled with:

```xml
<compress_output>true</compress_output>
```

Compression will be performed with zlib.

## Norms

In following section $x$ denotes solution vector.

### Manhattan norm

Norm1 is defined as Manhattan norm in [LinAlg library](https://doxygen.opengeosys.org/d6/dcd/namespacemathlib_1_1linalg#ac7415e1254b70c4015ccd3b6f2873338).

Defined with following equation:
$$
x_{norm}=\sum_{i=0}^{n}{|x_i|}
$$

### Euclidan norm

Norm2 in OpenGeoSys implementation in the [LinAlg library](https://doxygen.opengeosys.org/d6/dcd/namespacemathlib_1_1linalg#af298d1ddc92d7ce52046adc669e9904f).

Defined with following equation:
$$
x_{norm} = \sqrt{\sum_{i=0}^{n}{x_i}}
$$

### Max norm

$$
\mathrm{max}_i |x_i|
$$
