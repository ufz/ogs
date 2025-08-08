<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>geometry</name>
    <points>
        <point id="0" x="0."   y="0."   z="0."/>
        <point id="1" x="20."  y="0"    z="0."/>
        <point id="2" x="20."  y="20"   z="0."/>
        <point id="3" x="0."   y="20"   z="0."/>
        <point id="4" x="5.3019734184500145" y="0" z="0"/>
        <point id="5" x="14.628126581549985" y="20" z="0"/>
    </points>
    <polylines>
        <polyline id="0" name="bottom">
            <pnt>0</pnt>
            <pnt>1</pnt>
        </polyline>
        <polyline id="1" name="top">
            <pnt>2</pnt>
            <pnt>3</pnt>
        </polyline>
        <polyline id="2" name="left">
            <pnt>0</pnt>
            <pnt>3</pnt>
        </polyline>
        <polyline id="3" name="right">
            <pnt>1</pnt>
            <pnt>2</pnt>
        </polyline>
        <polyline id="4" name="topleft">
            <pnt>5</pnt>
            <pnt>3</pnt>
        </polyline>
        <polyline id="5" name="topright">
            <pnt>2</pnt>
            <pnt>5</pnt>
        </polyline>
        <polyline id="6" name="bottomleft">
            <pnt>0</pnt>
            <pnt>4</pnt>
        </polyline>
        <polyline id="7" name="bottomright">
            <pnt>4</pnt>
            <pnt>1</pnt>
        </polyline>
</polylines>
</OpenGeoSysGLI>
