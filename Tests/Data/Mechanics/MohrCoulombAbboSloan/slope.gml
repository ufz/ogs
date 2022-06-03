<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>
<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>slope_geometry</name>
    <points>
        <point id="0" x="0" y="0" z="0"/>
        <point id="1" x="18" y="0" z="0"/>
        <point id="2" x="18" y="6.5" z="0"/>
        <point id="3" x="15.5" y="6.5" z="0"/>
        <point id="4" x="12.5" y="6.5" z="0"/>
        <point id="5" x="12" y="6.5" z="0"/>
        <point id="6" x="3" y="2" z="0"/>
        <point id="7" x="0" y="2" z="0"/>
    </points>
    <polylines>
        <polyline id="0" name="bottom">
            <pnt>0</pnt>
            <pnt>1</pnt>
        </polyline>
        <polyline id="1" name="top">
            <pnt>3</pnt>
            <pnt>4</pnt>
        </polyline>
        <polyline id="2" name="right">
            <pnt>1</pnt>
            <pnt>2</pnt>
        </polyline>
        <polyline id="3" name="left">
            <pnt>0</pnt>
            <pnt>7</pnt>
        </polyline>
        <polyline id="4" name="slope">
            <pnt>5</pnt>
            <pnt>6</pnt>
        </polyline>
        <polyline id="5" name="slopefoot">
            <pnt>6</pnt>
            <pnt>7</pnt>
        </polyline>
    </polylines>
</OpenGeoSysGLI>
