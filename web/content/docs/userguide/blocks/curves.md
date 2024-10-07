+++
date = "2018-02-27T11:00:13+01:00"
title = "Curves"
author = "Feliks Kiszkurno"
weight = 10
+++
<div class="note">

### Work in progress

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

Curve represents a piece-wise linear function -- the values between coordinates are interpolated linearly though.
Outside of the curve's range of definition the curve becomes a constant using the first value on the left and last value on the
right.
