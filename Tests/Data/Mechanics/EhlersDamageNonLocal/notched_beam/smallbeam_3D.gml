<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://www.opengeosys.org/images/xsd/OpenGeoSysGLI.xsd" xmlns:ogs="http://www.opengeosys.org">
    <name>smallbeam_3D</name>
    <points>
        <point id="0" x="0" y="0" z="0" name="p1"/>
        <point id="1" x="0" y="0" z="0.02" name="p2"/>
        <point id="2" x="0.1022" y="0" z="0" name="p3"/>
        <point id="3" x="0.1022" y="0" z="0.02" name="p4"/>
        <point id="4" x="0.0511" y="0.0306" z="0" name="p5"/>
        <point id="5" x="0.0511" y="0.0306" z="0.02" name="p6"/>
        <point id="6" x="0.0" y="0.0306" z="0.02" name="p7"/>
        <point id="7" x="0.1022" y="0.0306" z="0.02" name="p8"/>

        <point id="8" x="0.0505" y="0.0306" z="0.0" name="p9"/>
        <point id="9" x="0.0517" y="0.0306" z="0.0" name="p10"/>
        <point id="10" x="0.0505" y="0.0306" z="0.02" name="p11"/>
        <point id="11" x="0.0517" y="0.0306" z="0.02" name="p12"/>

    </points>

    <polylines>
	    <polyline id="0" name="bot_sx">
		    <pnt>0</pnt>
		    <pnt>1</pnt>
	    </polyline>
	    <polyline id="1" name="bot_dx">
		    <pnt>2</pnt>
		    <pnt>3</pnt>
	    </polyline>
	    <polyline id="2" name="top">
		    <pnt>4</pnt>
		    <pnt>5</pnt>
	    </polyline>
    </polylines>

    <surfaces>
        <surface id="0" name="front">
            <element p1="0" p2="2" p3="6"/>
            <element p1="2" p2="6" p3="7"/>
        </surface>
        <surface id="1" name="tamper">
            <element p1="8" p2="9" p3="10"/>
            <element p1="10" p2="9" p3="11"/>
        </surface>
    </surfaces>

</OpenGeoSysGLI>
