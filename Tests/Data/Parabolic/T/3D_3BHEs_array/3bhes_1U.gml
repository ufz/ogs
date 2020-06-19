<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://www.opengeosys.org/images/xsd/OpenGeoSysGLI.xsd" xmlns:ogs="http://www.opengeosys.org">
    <name>3bhes_1U</name>
    <points>
        <point id="0" x="-25.0" y="0.0" z="0.0" name="TOP_left_down"/>
        <point id="1" x="-25.0" y="50.0" z="0.0" name="TOP_left_upper"/>
        <point id="2" x="25.0" y="50.0" z="0.0" name="TOP_right_upper"/>
        <point id="3" x="25.0" y="0.0" z="0.0" name="TOP_right_down"/>
        <point id="4" x="-25.0" y="0.0" z="-72.0" name="bottom_left_down"/>
        <point id="5" x="-25.0" y="50.0" z="-72.0" name="bottom_left_upper"/>
        <point id="6" x="25.0" y="50.0" z="-72.0" name="bottom_right_uppert"/>
        <point id="7" x="25.0" y="0.0" z="-72.0" name="bottom_right_down"/>

        <point id="8" x="0.0" y="25.0" z="-52.0" name="BHE1_BOTTOM"/>
        <point id="9" x="0.0" y="25.0" z="-2.0" name="BHE1_TOP"/>
        <point id="10" x="-6.0" y="25.0" z="-52.0" name="BHE2_BOTTOM"/>
        <point id="11" x="-6.0" y="25.0" z="-2.0" name="BHE2_TOP"/>
        <point id="12" x="0.0" y="25.0" z="-52.0" name="BHE3_BOTTOM"/>
        <point id="13" x="0.0" y="25.0" z="-2.0" name="BHE3_TOP"/>
    </points>
    <polylines>
        <polyline id="0" name="BHE_1">
            <pnt>9</pnt>
            <pnt>8</pnt>
        </polyline>
        <polyline id="1" name="BHE_2">
            <pnt>11</pnt>
            <pnt>10</pnt>
        </polyline>
        <polyline id="2" name="BHE_3">
            <pnt>13</pnt>
            <pnt>12</pnt>
        </polyline>
    </polylines>
    <surfaces>
        <surface id="0" name="top">
            <element p1="0" p2="1" p3="2"/>
            <element p1="0" p2="3" p3="2"/>
        </surface>
    </surfaces>
</OpenGeoSysGLI>
