<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://www.opengeosys.org/images/xsd/OpenGeoSysGLI.xsd" xmlns:ogs="http://www.opengeosys.org">
    <name>bar</name>
    <points>
        <point id="0" x="0" y="0" z="0" name="p_0"/>
        <point id="1" x="1.0" y="0" z="0" name="p_1"/>
        <point id="2" x="1.0" y="0.05" z="0" name="p_2"/>
        <point id="3" x="0" y="0.05" z="0" name="p_3"/>
        <point id="4" x="0" y="0" z="0.05" name="p_4"/>
        <point id="5" x="1.0" y="0" z="0.05" name="p_5"/>
        <point id="6" x="1.0" y="0.05" z="0.05" name="p_6"/>
        <point id="7" x="0" y="0.05" z="0.05" name="p_7"/>
    </points>

    <polylines>
    <polyline id="0" name="left_bot">
        <pnt>0</pnt>
        <pnt>4</pnt>
    </polyline>
    <polyline id="1" name="left_front">
        <pnt>0</pnt>
        <pnt>3</pnt>
    </polyline>
    </polylines>


    <surfaces>
        <surface id="0" name="front">
            <element p1="0" p2="1" p3="3"/>
            <element p1="1" p2="3" p3="2"/>
        </surface>
        <surface id="1" name="back">
            <element p1="4" p2="5" p3="7"/>
            <element p1="5" p2="7" p3="6"/>
        </surface>
        <surface id="2" name="left">
            <element p1="0" p2="4" p3="3"/>
            <element p1="4" p2="3" p3="7"/>
        </surface>
        <surface id="3" name="right">
            <element p1="1" p2="5" p3="2"/>
            <element p1="5" p2="2" p3="6"/>
        </surface>
        <surface id="4" name="top">
            <element p1="3" p2="2" p3="7"/>
            <element p1="2" p2="7" p3="6"/>
        </surface>
        <surface id="5" name="bottom">
            <element p1="0" p2="1" p3="4"/>
            <element p1="1" p2="4" p3="5"/>
        </surface>
    </surfaces>

</OpenGeoSysGLI>
