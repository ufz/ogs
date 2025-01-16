The default setting is false. When activated, non-linear iterations will be turned off and the equation system for the soil part is only assembled once.
If the `<use_algebraic_bc>` option is also set to true, the simulation will just go through a single iteration within a time step.
Otherwise iterations are needed to fulfill the boundary conditions of the BHE.
Only when fluid properties do not change much with temperature, this feature can be activated.
The effect of skipping the non-linear iteration is much faster simulation speed.
