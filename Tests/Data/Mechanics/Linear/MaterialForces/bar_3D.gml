<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>bar_geom</name>
    <points>
        <point id="0" x="0" y="0" z="0" name="origin"/>
        <point id="1" x="1.0" y="0" z="0"/>
        <point id="2" x="1.0" y="0.05" z="0"/>
        <point id="3" x="0" y="0.05" z="0"/>

        <point id="4" x="0" y="0" z="0.05"/>
        <point id="5" x="1.0" y="0" z="0.05"/>
        <point id="6" x="1.0" y="0.05" z="0.05"/>
        <point id="7" x="0" y="0.05" z="0.05"/>
    </points>

    <polylines>
        <polyline id="0" name="p_0">
            <pnt>0</pnt>
            <pnt>4</pnt>
        </polyline>
        <polyline id="1" name="p_1">
            <pnt>1</pnt>
            <pnt>5</pnt>
        </polyline>
        <polyline id="2" name="left_vertical">
            <pnt>0</pnt>
            <pnt>3</pnt>
        </polyline>
        <polyline id="3" name="right_vertical">
            <pnt>1</pnt>
            <pnt>2</pnt>
        </polyline>
    </polylines>

    <surfaces>
        <surface id="0" name="left"><!-- x=0 -->
            <element p1="0" p2="3" p3="4"/>
            <element p1="4" p2="3" p3="7"/>
        </surface>
        <surface id="0" name="right"><!-- x=0 -->
            <element p1="1" p2="2" p3="5"/>
            <element p1="5" p2="2" p3="6"/>
        </surface>
    </surfaces>
</OpenGeoSysGLI>
