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

This block defines all aspects related to the progress of time in the simulation and controls when the output files will be
written, their naming pattern and content.

## Process

First section of this block, the one controlling time progress has to be set up separately for each process defined in
[processes block](/docs/userguide/blocks/processes/#process-variables) separately.
There can be more than one process defined.
They are distinguished with the parameter "ref" in the `<process>` tag:

```xml
<process ref="Process_name">
    ...
</process>
```

### Time stepping

Numerical simulations proceed in time steps, if they are transient.
Their size is an important factor impacting, if a stable solution will be found.

In OpenGeoSys the following time step definitions are available:

- FixedTimeStepping
- SingleStep
- IterationNumberBasedTimeStepping
- EvolutionaryPIDcontroller

#### Single time step

This is the simplest type of time steppings used for solution of elliptic PDEs as in the `SteadyStateDiffusion` process.
It only requires its type to be given and is used for stead-state problems:

```xml
<time_stepping>
    <type>SingleStep</type>
</time_stepping>
```

#### Fixed time stepping

This type of time stepping will run the simulation with a predetermined number of steps from the start to the end of the
simulation.
The steps are defined with parameter pairs giving the number of repetitions and the length of an individual step. It is allowed
to define multiple pairs. This way, the time step can vary throughout the experiment.

In the example below, a duration of $1$ time unit is simulated, with first 5 time steps being $0.1$ time unit long and next 10
being $0.05$ time unit long.

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
If this is not the case, a time step at `t_end` will be added.
The first time step of the simulation is always at `t_initial`.

The unit for time has to be consistent with the rest of the units used in the experiment.
Following the SI system, second is the common choice as time unit.
This is how to define the time-stepping for an experiment simulating a year of time in daily time step using second as time
unit:

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

Fixed time stepping can be also used to create linearly spaced steps between `t_initial` and `t_end`.
The number of steps is defined by `n_steps`.

```xml
<time_stepping>
    <type>FixedTimeStepping</type>
    <t_initial>0</t_initial>
    <t_end>10</t_end>
    <n_steps>10</n_steps>
</time_stepping>
```

will result in 11 time steps at: 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10.

As fixed time stepping is one of the simplest available, it is a good starting point.

#### Iteration number based time stepping

[See Doxygen](https://doxygen.opengeosys.org/df/df4/ogs_file_param__prj__time_loop__processes__process__time_stepping__iterationnumberbasedtimestepping)
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

Iteration number based time stepping corresponds to adaptive time steps. Experiences users may prefer this time-stepping
especially for the solution of non-linear problems.

#### Evolutionary PID controller

[See Doxygen](https://doxygen.opengeosys.org/d3/d86/ogs_file_param__prj__time_loop__processes__process__time_stepping__evolutionarypidcontroller)

### Error tolerances

Error tolerances will be applied to the solution vector.
There are two ways of defining error tolerances:

- relative
- absolute

Both of them can be defined as single value, that will be applied to all process variables, or with multiple ones applied to
them individually ("per component" criteria).

In those two cases appropriate tags have to be used:

- for relative tolerances: `<reltol> </reltol>`
- for absolute tolerances: `<abstol> </abstol>`

and for "per component":

- for relative tolerances: `<reltols> </reltols>`
- for absolute tolerances: `<abstols> </abstols>`

<!-- TODO: Describe the definition of the relative tolerance. -->

If tolerances per process variable are provided, the order of values defined inside of the tags `<abstols> </abstols>` and
`<reltols> </reltols>` has to match order of process variables defined in [processes](/docs/userguide/blocks/processes/).
If process variable is directional, an XYZ order is followed.
For example for a 3D THM problem with directional displacement $\mathbf{u}$ the following order has to be used: $T$, $p$, $u_x$
, $u_y$, $u_z$.
The order of $T$, $p$, and $u$ is prescribed by the THM process.

Depending on process and simulation setup, the number of variables in the solution vector can vary.
For example in a 2D THM simulation, there will four tolerances given in `abstol`: one for $T$, one for $p$, and two for $\mathbf{u}$ ($u_x$
and $u_y$).

The same order is preserved for [output variables](/docs/userguide/blocks/time_loop/#output).

### Convergence criteria

In this part the selection of a convergence criterion is described. The following options are available:

- [DeltaX and PerComponentDeltaX](/docs/userguide/blocks/time_loop/#deltax-and-percomponentdeltax)
- [Residual and PerComponentResidual](/docs/userguide/blocks/time_loop/#residual-and-percomponentresidual)

All of the criteria mentioned above compare a value quantifying the error (a residual or a discrete change in time) with a
user-defined tolerance.

There are three way how error tolerances can be set up:

- only relative tolerances $e_{\mathrm{rel}} \le \epsilon_{\mathrm{rel}}$
- only absolute tolerances $e_{\mathrm{abs}} \le \epsilon_{\mathrm{abs}}$
- relative or absolute tolerances $e_{\mathrm{abs}} \le \epsilon_{\mathrm{abs}} \vee e_{\mathrm{rel}} \le \epsilon_{\mathrm{rel}} $

where $e$ indicates the observed error and $\epsilon$ indicates the defined tolerance.

In the last case the convergence criterion will be satisfied with **at least one** of the tolerances being met.

Both relative and absolute tolerances can be defined as a single value or using separate values for each component, see [error tolerances](/docs/userguide/blocks/time_loop/#error-tolerances) for more details.

#### DeltaX and PerComponentDeltaX

DeltaX criteria are based on a norm applied to the whole solution or increment vector. It reflects the difference of the initial value and its update in time for the next iteration according to a defined norm.

Therefore, the error is defined as difference between the results vector obtained from two consecutive iterations.
If this difference drops below the defined tolerance, the criterion is fulfilled and the solver will proceed to the next time step.

The DeltaX criterion usually uses the [euclidean norm](/docs/userguide/blocks/time_loop/#euclidan-norm) as indicated by the `<norm_type>` tag value `NORM2`.
Other possible values are `NORM1` and `INFINITY_N`.  

The norm of the solution vector is compared with the error tolerances.
They can be provided as absolute, relative, or both.

In opposite to DeltaX, PerComponentDeltaX is applied to an individual variable, or to each component of a vectorial variable. A PerComponentDeltaX criterion will be fulfilled, if the resulting variable under consideration (or all resulting components of a vectorial variable) are close enough to a user-defined criterion of maximum change between two consecutive time steps.

#### Residual and PerComponentResidual

Differently than [DeltaX](/docs/userguide/blocks/time_loop/#deltax-and-percomponentdeltax), the residual based convergence
criteria use the residuum values for evaluation, not the increment values.

The residual value is the difference between two sides of modeled equation.
In analytical solution, they would be perfectly equal, but in numerical solution, this will usually not be the case.
The lower the residuals are, the closer numerical solution is to fulfilling the equation that is being solved.

In this criteria, error is defined as imbalance between left and right hand sides of the equation.
Once the residual values drop below user defined tolerances, solution is accepted and solver can proceed to the next time step.

The normal of the residuum vector is compared with the error tolerances.
They can be provided as absolute, relative or both.

### Norms

In following section $x$ denotes solution vector.

#### Absolute-value norm

`NORM1` is defined as the absolute-value norm in the [LinAlg library](https://doxygen.opengeosys.org/d6/dcd/namespacemathlib_1_1linalg#ac7415e1254b70c4015ccd3b6f2873338).

It is defined by the following equation:
$$
x_{norm}=\sum_{i=0}^{n}{|x_i|}
$$

#### Euclidean norm

`NORM2` in the OpenGeoSys implementation is as well in the [LinAlg library](https://doxygen.opengeosys.org/d6/dcd/namespacemathlib_1_1linalg#af298d1ddc92d7ce52046adc669e9904f).

It is defined by the following equation:
$$
x_{norm} = \sqrt{\sum_{i=0}^{n}{x^2_i}}
$$

#### Infinity norm

Moreover, `INFINITY_N` is the infinity norm (sometimes also called maximum norm), which can be found as well in the [LinAlg Library](https://doxygen.opengeosys.org/d6/dcd/namespacemathlib_1_1linalg#a49abd74780cb8a1a7135a722ca762394) and which reads as follows:

$$
x_{norm} = \mathrm{max} ((|x_i|)_{i = 1,...,n})
$$

### Time discretization

<!-- TODO: describe different options for time discretization -->

For the time being only backward Euler time stepping scheme is available.

## Output

In this section all parameters controlling the output of the data from the simulation are defined.

Currently OpenGeoSys supports two output file formats: VTK and XDMF.
The desired format for output files can be declared as follows:

```xml
<type>VTK</type>
```

The number of time steps at which the output files are written doesn't need to match the simulation time steps.

There are two ways how the times at which the output will be written can be specified: as fixed output times and time steps.

The fixed output times can be defined in the same way as in [fixed time step time stepping](/docs/userguide/blocks/time_loop/#fixed-time-stepping).

This approach is suitable, if specific times of the simulation are of a special interest to the user.
In such a case, OpenGeoSys can write the output files at times provided in a list within the `<fixed_output_times> </fixed_output_times>` tag:

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

If output times are supposed to be cyclic, or follow some patterns, defining them as time steps will be more convenient.
In order to do so, block ```xml <timesteps> </timesteps>```has to be placed in ```<output>```block.
It can contain arbitrary number of blocks ```xml <pair> </pair>```, which define a cyclical pattern each.
The block ```<pair>``` has to contain two tags: ```<repeat>``` and ```<each_steps>```.
```<repeat>``` defines how many outputs with spacing of ```<each_steps>``` time steps should be written.
In a sequence of <pair> blocks the starting time step for a pair is the last time step from the previous one.
Following example illustrates this:

```xml
<timesteps>
    <pair> <!-- First pair -->
        <repeat>1</repeat>
        <each_steps>10</each_steps>
    </pair>
    <pair> <!-- Second pair -->
        <repeat>1</repeat>
        <each_steps>90</each_steps>
    </pair>
    <pair> <!-- Third pair -->
        <repeat>1</repeat>
        <each_steps>900</each_steps>
    </pair>
</timesteps>
```

Using this block in the project file will result in the output being written at the time steps: 10, 100 and 1000.
The first file would be written as a result of the fist ```<pair>```.
It will write output once (as the tag ```<repeat>``` has value of 1) after 10 time steps (10 is the value provided in the tag ```<each_steps>```).
The second pair will write an output file after 90 time steps.
However, the counting starts not at time step 0, but at the last time step resulting from the previous pair.
In this case this is time step 10, hence second pair will write an output at 100 (10 time steps from first pair plus 90 time steps from second pair).
The same applies to the third pair.
It will write an output at 1000 as 1000 is the result of addition of 10, 90, and 900 from first, second and third steps respectively.
Note, that specifically in this case the values in ```<each_steps>``` tag can be summed directly, as each pair is repeated only once.

Now, consider a more complicated pattern.
Let's assume, that the output has to be written at every time step for the range 1-10 time steps, every tenth in the range 10-100 time steps and every hundredth for the range 100-1000.
Block defining this pattern would be written as follows:

```xml
<timesteps>
    <pair> <!-- First pair -->
        <repeat>10</repeat>
        <each_steps>1</each_steps>
    </pair>
    <pair> <!-- Second pair -->
        <repeat>9</repeat>
        <each_steps>10</each_steps>
    </pair>
    <pair> <!-- Third pair -->
        <repeat>9</repeat>
        <each_steps>100</each_steps>
    </pair>
</timesteps>
```

This will result in output being written at time steps: 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900 and 1000.
Note, that since output at 10th time step is the last one resulting from first pair, it will be the starting point for the second pair, therefore the second pair needs to be repeated only 9 times.
Otherwise the last output from the second pair would be written at time step 110.
The same principle applies to the third pair.
Output at time step 100 is written by the second pair, therefore to get to 1000, third pair only needs to be repeated 9 times.

Regardless of the user defined output time steps, OpenGeoSys will write the output files at $t=0$ and $t=t_{end}$.

The naming of output files written at different time steps is handled automatically by OpenGeoSys.
The user can use two tags that allow to control how the file names are structured and what they contain: tags `<suffix> </suffix>` and `<prefix> </prefix>`.
The rules regarding their content are the same.
They should contain text strings.
Variables can be used to allow those strings to vary between files.
The following variables can be called:

- `meshname`
- `timestep`
- `time`
- `iteration`

They can be used in the file name with the syntax illustrated by the following example:

```xml
<suffix>Experiment_name</suffix>
<prefix>_at_time_step_{:timestep}_time_{:time}_iteration_{:iteration}</prefix>
```

The name of each file will be automatically appended by a file extension appropriate to the content of `<type> </type>` tag.
No user action is needed in this regard.

The block `<variables> </variables>` allows users to control what process variables are included in the output files.
If only some of them are of interest, it may be convenient to reduce the size of the output files by only writing the necessary
ones.
If only temperature and pressure are important, this block can be defined as follows:

```xml
<variables>
    <variable>temperature</variable>
    <variable>pressure</variable>
</variables>
```

The list of available variables differs between processes.

<!-- TODO: Create list of variables available in different processes -->

In order to save space, the compression of output files can be enabled with:

```xml
<compress_output>true</compress_output>
```

Compression will be performed with zlib.
