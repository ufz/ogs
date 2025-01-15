<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>geometry</name>
    <points>
         <point id="0"  x="0." y="0."  z="0." />
         <point id="1"  x="0.001" y="0.0"  z="0." />
         <point id="2"  x="0.001" y="0.001"  z="0."/>
         <point id="3"  x="0." y="0.001"  z="0." />
         <point id="4"  x="0." y="0."  z="0.001" />
         <point id="5"  x="0.001" y="0.0"  z="0.001" />
         <point id="6"  x="0.001" y="0.001"  z="0.001"/>
         <point id="7"  x="0." y="0.001"  z="0.001" />
         <point id="8"  x="0.0005" y="0."  z="0.001" />
         <point id="9"  x="0.0005" y="0.0005"  z="0.001" />
         <point id="10" x="0.0" y="0.0005"  z="0.001"/>
    </points>

    <surfaces>
        <surface id="0" name="y0">
            <element p1="0" p2="1" p3="5"/>
            <element p1="0" p2="5" p3="4"/>
        </surface>
        <surface id="1" name="x0">
            <element p1="0" p2="4" p3="3"/>
            <element p1="4" p2="7" p3="3"/>
        </surface>
        <surface id="2" name="z0">
            <element p1="0" p2="3" p3="1"/>
            <element p1="3" p2="2" p3="1"/>
        </surface>
        <surface id="3" name="z1">
            <element p1="4" p2="5" p3="7"/>
            <element p1="5" p2="6" p3="7"/>
        </surface>
        <surface id="4" name="z1_indent">
            <element p1="4" p2="8" p3="9"/>
            <element p1="4" p2="9" p3="10"/>
        </surface>
    </surfaces>
</OpenGeoSysGLI>
