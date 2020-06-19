<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>cube_1x1x1_geometry</name>
    <points>
        <point id="0" x="0" y="0" z="0"/>
        <point id="1" x="0" y="0" z="1"/>
        <point id="2" x="0" y="1" z="1"/>
        <point id="3" x="0" y="1" z="0"/>

        <point id="4" x="1" y="0" z="0"/>
        <point id="5" x="1" y="0" z="1"/>
        <point id="6" x="1" y="1" z="1"/>
        <point id="7" x="1" y="1" z="0"/>

        <point id="8" x="0" y="0" z="0.1"/>
        <point id="9" x="0" y="0" z="0.2"/>
        <point id="10" x="0" y="0" z="0.3"/>
        <point id="11" x="0" y="0" z="0.4"/>
        <point id="12" x="0" y="0" z="0.5"/>
        <point id="13" x="0" y="0" z="0.6"/>
        <point id="14" x="0" y="0" z="0.7"/>
        <point id="15" x="0" y="0" z="0.8"/>
        <point id="16" x="0" y="0" z="0.9"/>

        <point id="17" x="0" y="1" z="0.1"/>
        <point id="18" x="0" y="1" z="0.2"/>
        <point id="19" x="0" y="1" z="0.3"/>
        <point id="20" x="0" y="1" z="0.4"/>
        <point id="21" x="0" y="1" z="0.5"/>
        <point id="22" x="0" y="1" z="0.6"/>
        <point id="23" x="0" y="1" z="0.7"/>
        <point id="24" x="0" y="1" z="0.8"/>
        <point id="25" x="0" y="1" z="0.9"/>
    </points>

    <surfaces>
        <surface id="0" name="left"><!-- x=0 -->
            <element p1="0" p2="1" p3="2"/>
            <element p1="0" p2="3" p3="2"/>
        </surface>
        <surface id="0" name="left0">
            <element p1="0" p2="8" p3="17"/>
            <element p1="0" p2="17" p3="3"/>
        </surface>
        <surface id="0" name="left1">
            <element p1="8" p2="9" p3="18"/>
            <element p1="8" p2="18" p3="17"/>
        </surface>
        <surface id="0" name="left2">
            <element p1="9" p2="10" p3="19"/>
            <element p1="9" p2="19" p3="18"/>
        </surface>
        <surface id="0" name="left3">
            <element p1="10" p2="11" p3="20"/>
            <element p1="10" p2="20" p3="19"/>
        </surface>
        <surface id="0" name="left4">
            <element p1="11" p2="12" p3="21"/>
            <element p1="11" p2="21" p3="20"/>
        </surface>
        <surface id="0" name="left5">
            <element p1="12" p2="13" p3="22"/>
            <element p1="12" p2="22" p3="21"/>
        </surface>
        <surface id="0" name="left6">
            <element p1="13" p2="14" p3="23"/>
            <element p1="13" p2="23" p3="22"/>
        </surface>
        <surface id="0" name="left7">
            <element p1="14" p2="15" p3="24"/>
            <element p1="14" p2="24" p3="23"/>
        </surface>
        <surface id="0" name="left8">
            <element p1="15" p2="16" p3="25"/>
            <element p1="15" p2="25" p3="24"/>
        </surface>
        <surface id="0" name="left9">
            <element p1="16" p2="1" p3="2"/>
            <element p1="16" p2="2" p3="25"/>
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
