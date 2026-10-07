+++
date = "2018-02-27T11:00:13+01:00"
title = "Curves"
author = "Feliks Kiszkurno"
weight = 10
+++
<div class="note">

## Work in progress

This page is a work in progress.

It was published in this state to make existing content available to users and highlight missing parts to contributors.

**Contributors:** please see Documentation Contribution Guide to contribute to the documentation.

**Users:** the content of this page has been verified and is correct. Please return later for more content!

</div>

<!-- TODO: Add general description -->

This block can be used to define one or more ...

They can be accessed from within `<type>Curve</type>` or `<type>CurveScaled</type>` parameters (more details on [curves](/docs/userguide/blocks/parameters/#curve) and [curves scaled](/docs/userguide/blocks/parameters/#curvescaled)).

Following variables are accessible from CurveScaled:

- spatial: x, y, z
- temporal: t

## Example

Curves can be defined using the following block:

```xml
<curves>
    <curve>
        <name>CurveName</name>
        <coords>
            x_1 x_2 ... x_n
        </coords>
        <values>
            y(x_1) y(x_2) ... y(x_n)
        </values>
    </curve>
</curves>
```

or recommended for long curves definitions, they can be read from a binary file (double precision values in little endian format) placed in the project file directory:

```xml
<curves>
    <curve>
        <name>CurveName</name>
        <read_from_file>true</read_from_file>
        <coords>
            coords_filename.bin
        </coords>
        <values>
            values_filename.bin
        </values>
    </curve>
</curves>
```

They can be called in "Properties" and "Parameters" inside "Expression":
`CurveName(evaluation_value)`

where `evaluation_value` is always exactly one value and refers to values provided inside `<coords> </coords>`.
In the `<expression>`s using curves, time $t$ and global coordinate $x$, $y$, and $z$ can be used, and passed as curve's argument for example:

```xml
<parameter>
    <name>density</name>
    <type>Function</type>
    <expression>2000 - CurveName(x^2 + y^2) * TimeDependentCurve(t)</expression>
</parameter>
```

Curves can be useful to include information on variation of time dependent parameters.
Those values can be called from the 'equation' tag in 'parameter' if 'type' 'formula' is selected.

## Interpolation

By default, a curve represents a piece-wise linear function -- the values between coordinates are interpolated linearly.
Outside of the curve's range of definition the curve becomes a constant using the first value on the left and last value on the
right.

With the optional tag `<interpolation>pchip</interpolation>` inside `<curve>`, the values are interpolated by a shape-preserving piecewise cubic Hermite polynomial (PCHIP) instead.
The curve is then continuously differentiable within the range of the supporting points and monotone wherever the data is monotone, i.e., it does not over- or undershoot the supporting values.
At least four supporting points are required.
This gives a continuous derivative, e.g., for the Jacobian of properties defined by curves.

```xml
<curve>
    <name>CurveName</name>
    <interpolation>pchip</interpolation>
    <coords>x_1 x_2 ... x_n</coords>
    <values>y(x_1) y(x_2) ... y(x_n)</values>
</curve>
```

By default the slope at the two end points is the slope of the adjacent interval.
It can be set with the optional tags `<left_end_slope>` and `<right_end_slope>`.
To keep the curve monotone on the end intervals, a given slope must have the sign of the slope of the adjacent interval and at most three times its magnitude (0 for a flat interval); other values are rejected.
A zero slope at the first point is useful for retention curves, which are constant (slope 0) below their first point:
the derivative of the saturation, and thus the capacity of the Richards equation, is then continuous at the saturated edge of the curve instead of jumping there.

```xml
<curve>
    <name>capillary_pressure</name>
    <interpolation>pchip</interpolation>
    <left_end_slope>0</left_end_slope>
    <coords>...</coords>
    <values>...</values>
</curve>
```
