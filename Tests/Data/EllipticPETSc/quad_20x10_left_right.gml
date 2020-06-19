<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>quad_20x10_left_right_geometry</name>
    <points>
        <point id="0" x="0.000000" y="0.000000" z="0.000000"/>
        <point id="1" x="0.000000" y="10.000000" z="0.000000"/>
        <point id="2" x="20.000000" y="0.000000" z="0.000000"/>
        <point id="3" x="20.000000" y="10.000000" z="0.000000"/>
    </points>

    <polylines>
        <polyline id="0" name="left">
            <pnt>0</pnt>
            <pnt>1</pnt>
        </polyline>
        <polyline id="1" name="right">
            <pnt>2</pnt>
            <pnt>3</pnt>
        </polyline>
    </polylines>
</OpenGeoSysGLI>
