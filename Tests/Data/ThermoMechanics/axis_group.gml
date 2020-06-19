<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>square_1x1_geometry</name>
    <points>
        <point id="0" x="10" y="0" z="0"/>
        <point id="1" x="12.5" y="0" z="0"/>
        <point id="2" x="35" y="0" z="0"/>
        <point id="3" x="10" y="0.5" z="0"/>
        <point id="4" x="12.5" y="0.5" z="0"/>
        <point id="5" x="35" y="0.5" z="0"/>
    </points>

    <polylines>
        <polyline id="0" name="in">
            <pnt>0</pnt>
            <pnt>3</pnt>
        </polyline>
        <polyline id="1" name="out">
            <pnt>2</pnt>
            <pnt>5</pnt>
        </polyline>
	<polyline id="2" name="bottom">
            <pnt>0</pnt>
            <pnt>2</pnt>
        </polyline>
        <polyline id="3" name="top">
            <pnt>3</pnt>
            <pnt>5</pnt>
        </polyline>
    </polylines>
</OpenGeoSysGLI>
