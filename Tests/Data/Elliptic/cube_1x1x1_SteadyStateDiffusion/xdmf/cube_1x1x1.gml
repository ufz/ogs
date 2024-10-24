<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>cube_1x1x1_geometry</name>
    <points>
        <point id="0" x="0" y="0" z="0" name="origin"/>
        <point id="1" x="0" y="0" z="1" name="left_front_upper"/>
        <point id="2" x="0" y="1" z="1" name="left_back_upper"/>
        <point id="3" x="0" y="1" z="0" name="left_back_lower"/>

        <point id="4" x="1" y="0" z="0" name="right_front_lower"/>
        <point id="5" x="1" y="0" z="1" name="right_front_upper"/>
        <point id="6" x="1" y="1" z="1" name="right_back_upper"/>
        <point id="7" x="1" y="1" z="0" name="right_back_lower"/>

        <point id="8" x="0.1" y="0.1" z="0"/>
        <point id="9" x="0.3" y="0.1" z="0"/>
        <point id="10" x="0.3" y="0.3" z="0"/>
        <point id="11" x="0.1" y="0.3" z="0"/>
    </points>

    <polylines>
        <polyline id="0" name="xy_rectangle_border">
            <pnt>8</pnt>
            <pnt>9</pnt>
            <pnt>10</pnt>
            <pnt>11</pnt>
            <pnt>8</pnt>
        </polyline>
    </polylines>

    <surfaces>
        <surface id="0" name="left"><!-- x=0 -->
            <element p1="0" p2="1" p3="2"/>
            <element p1="0" p2="3" p3="2"/>
        </surface>
        <surface id="1" name="right"><!-- x=1 -->
            <element p1="4" p2="6" p3="5"/>
            <element p1="4" p2="6" p3="7"/>
        </surface>
        <surface id="2" name="top"><!-- z=1 -->
            <element p1="1" p2="2" p3="5"/>
            <element p1="5" p2="2" p3="6"/>
        </surface>
        <surface id="3" name="bottom"><!-- z=0 -->
            <element p1="0" p2="3" p3="4"/>
            <element p1="4" p2="3" p3="7"/>
        </surface>
        <surface id="4" name="front"><!-- y=0 -->
            <element p1="0" p2="1" p3="4"/>
            <element p1="4" p2="1" p3="5"/>
        </surface>
        <surface id="5" name="back"><!-- y=1 -->
            <element p1="2" p2="3" p3="6"/>
            <element p1="6" p2="3" p3="7"/>
        </surface>
    </surfaces>
</OpenGeoSysGLI>
