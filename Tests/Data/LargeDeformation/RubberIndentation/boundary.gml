<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>geometry</name>
    <points>
         <point id="0"  x="0." y="0."  z="0." />
         <point id="1"  x="0.01" y="0.0"  z="0." />
         <point id="2"  x="0.01" y="0.01"  z="0." />
         <point id="3"  x="0." y="0.01"  z="0." />
         <point id="4"  x="0.005" y="0.01"  z="0." />
    </points>

    <polylines>
        <polyline id="1" name="left">
            <pnt>0</pnt>
            <pnt>3</pnt>
        </polyline>
        <polyline id="2" name="right">
            <pnt>1</pnt>
            <pnt>2</pnt>
        </polyline>
        <polyline id="3" name="bottom">
            <pnt>0</pnt>
            <pnt>1</pnt>
        </polyline>
        <polyline id="4" name="indent_part">
            <pnt>3</pnt>
            <pnt>4</pnt>
        </polyline>
        <polyline id="5" name="top">
            <pnt>3</pnt>
            <pnt>2</pnt>
        </polyline>
    </polylines>
</OpenGeoSysGLI>
