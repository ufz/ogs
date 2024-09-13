<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>geometry</name>
    <points>
         <point id="0"  x="0." y="0."  z="0." />
         <point id="1"  x="48.0e-3" y="44.0e-3"  z="0." />
         <point id="2"  x="48.0e-3" y="60.0e-3"  z="0." name="pnt0" />
         <point id="3"  x="0." y="44.0e-3"  z="0." />
         <point id="4"  x="0." y="0."  z="0.001" />
         <point id="5"  x="48.0e-3" y="44.0e-3"  z="0.001" />
         <point id="6"  x="48.0e-3" y="60.0e-3"  z="0.001" name="pnt2" />
         <point id="7"  x="0." y="44.0e-3"  z="0.001" />
    </points>

    <surfaces>
        <surface id="0" name="left">
            <element p1="0" p2="4" p3="3"/>
            <element p1="4" p2="7" p3="3"/>
        </surface>
        <surface id="1" name="right">
            <element p1="2" p2="6" p3="5"/>
            <element p1="1" p2="2" p3="5"/>
        </surface>
        <surface id="3" name="z0">
            <element p1="0" p2="1" p3="3"/>
            <element p1="1" p2="2" p3="3"/>
        </surface>
        <surface id="4" name="z1">
            <element p1="4" p2="5" p3="7"/>
            <element p1="5" p2="6" p3="7"/>
        </surface>
    </surfaces>
</OpenGeoSysGLI>
