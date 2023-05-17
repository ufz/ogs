+++
date = "2018-02-27T11:00:13+01:00"
title = "Glossary"
author = "Feliks Kiszkurno"
weight = 1
+++

This file will list in alphabetical order some most important terms, that may not be know to wider audience or are used in non
standard way in the OpeGeoSys project.

The terms explained in this section are relevant for the "Beginners Guide" in other sections (for example in process
descriptions), the meaning of some terms used here may vary.

## A

### Aqueous phase

Liquid phase in porous media. Can be used in THM process....

## B

## C

### Constitutive relation

Known also as constitutive equation. It is an equation defining a relation between physical quantities that is specific to
material for which the equation is valid and describes how it responds to external stimulation. For more information click [here](https://en.wikipedia.org/wiki/Constitutive_equation) and for more information on how they are used in OpenGeoSys click [here](/docs/userguide/blocks/processes/#constitutive-relations), for the list of relations available in OpenGeoSys and details on
them click [here](/docs/userguide/blocks/misc/constitutive_relations).

### Constitutive equation

See [constitutive relation](/docs/userguide/troubleshooting/glossary/#constitutive-relation).

## D

## E

## F

## G

## H

## I

### Interface

Boundary between two media (geological layers, construction elements, etc.) defined entirely within the simulation domain.

## J

## K

## L

## M

### Medium

Part of the simulation domain characterized by the same physical properties and phase composition. One medium can consists of multiple phases.  

## N

### NaN

NaN is an acronym for not-a-number. This value is a result of:

- dividing over 0
- taking a square root of negative value

### Native installation

Setup in which program (e.g. OpenGeoSys platform) is installed directly in the system running on a specific machine without
containers or other intermediary.

## O

## P

### Primary variable

### Process variable

### Project file

Input file provided to OGS executable as the required input argument. It defines all aspects of the simulation that is to be
executed.

## Q

## R

### Residuals

Residuals is the difference between left and right hand side of the modeled equation. In theory, if a perfect solution was
available, both sides would be equal. Due to numerical errors, it is usually not the case.

## S

### Solution vector

### Secondary variable

## T

## U

## V

### Vectorial process variable

Is a subtype of [process variable](/docs/userguide/troubleshooting/glossary/#process-variable) which is defined as vector. It
is mostly directional. An example of vector process variable is displacement which can have a different value in each direction:
$$
\mathbf{u} = [u_x u_y u_z]
$$

## X

## Y

## Z
