<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>geometry</name>
    <points>
        <point id="0" x="0" y="0" z="0"/>
        <point id="1" x="3000" y="0" z="0"/>
        <point id="2" x="6000" y="0" z="0"/>
        <point id="3" x="0" y="-300" z="0"/>
        <point id="4" x="6000" y="-300" z="0"/>
        <point id="5" x="0" y="-3000" z="0"/>
        <point id="6" x="6000" y="-3000" z="0"/>
        <point id="7" x="6000" y="-1950" z="0"/>
    </points>
    <polylines>
        <polyline id="0" name="left_top_edge">
            <pnt>0</pnt>
            <pnt>3</pnt>
        </polyline>
        <polyline id="0" name="right_top_edge">
            <pnt>2</pnt>
            <pnt>4</pnt>
        </polyline>
        <polyline id="0" name="right_up_edge">
            <pnt>2</pnt>
            <pnt>7</pnt>
        </polyline>
        <polyline id="0" name="top_free">
            <pnt>0</pnt>
            <pnt>1</pnt>
        </polyline>
        <polyline id="0" name="top_ice">
            <pnt>1</pnt>
            <pnt>2</pnt>
        </polyline>
        <polyline id="0" name="top">
            <pnt>0</pnt>
            <pnt>2</pnt>
        </polyline>
        <polyline id="0" name="right">
            <pnt>0</pnt>
            <pnt>5</pnt>
        </polyline>
        <polyline id="0" name="left">
            <pnt>2</pnt>
            <pnt>6</pnt>
        </polyline>
        <polyline id="0" name="bottom">
            <pnt>5</pnt>
            <pnt>6</pnt>
        </polyline>
    </polylines>
</OpenGeoSysGLI>
