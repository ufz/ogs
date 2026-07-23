Configuration of the Anderson acceleration for the Picard iteration.

Anderson acceleration mixes the last \f$ m \f$ (damped) Picard steps instead of
taking the most recent one only, which typically reduces the number of
iterations for slowly converging fixed-point problems.

**Example:**

```xml
<anderson>
    <depth>3</depth>
</anderson>
```

Omitting the `anderson` subtree disables the acceleration and gives a plain
(possibly damped) Picard iteration.
