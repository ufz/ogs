<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>theis_geo</name>
    <points>
        <point id="0" x="3.048e-1" y="0" z="0"/>
        <point id="1" x="3.048e-1" y="1" z="0"/>
        <point id="2" x="3.048e+2" y="0" z="0"/>
        <point id="3" x="3.048e+2" y="1" z="0"/>
        <point id="4" x="9.639"   y="1" z="0" name="mpt"/>
    </points>
    <polylines>
        <polyline id="0" name="well">
            <pnt>0</pnt>
            <pnt>1</pnt>
        </polyline>
        <polyline id="1" name="infinit">
            <pnt>2</pnt>
            <pnt>3</pnt>
        </polyline>
    </polylines>
</OpenGeoSysGLI>

