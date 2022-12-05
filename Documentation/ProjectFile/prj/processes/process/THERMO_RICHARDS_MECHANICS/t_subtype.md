Specifies the constitutive setting used by the TRM process.

Currently, two constitutive settings are implemented: The
<tt>Stress_StrainTemperature</tt> one and one named
<tt>StressSaturation_StrainPressureTemperature</tt>. The latter allows for
strong THM coupling effects in the solid material behaviour, which computes
saturation and total stress from strain, liquid pressure and temperature. It has
to be implemented as an MFront material behaviour.

The <tt>Stress_StrainTemperature</tt> behaviour does not allow for such strong
coupling. There the solid material behaviour only computes effective stress from
strain and temperature. It can be implemented as a OGS built-in material
behaviour or via MFront.
