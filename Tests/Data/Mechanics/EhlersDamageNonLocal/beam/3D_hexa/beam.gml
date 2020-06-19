<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://www.opengeosys.org/images/xsd/OpenGeoSysGLI.xsd" xmlns:ogs="http://www.opengeosys.org">
    <name>beam</name>
    <points>
        <point id="0" x="0" y="0" z="0" name="p_0"/>
        <point id="1" x="1" y="0" z="0" name="p_1"/>
        <point id="2" x="1" y="0.25" z="0" name="p_2"/>
        <point id="3" x="0" y="0.25" z="0" name="p_3"/>
        <point id="4" x="0.5" y="0.25" z="0" name="p_4"/>

        <point id="5" x="0" y="0" z="0.05" name="p_5"/>
        <point id="6" x="1" y="0" z="0.05" name="p_6"/>
        <point id="7" x="1" y="0.25" z="0.05" name="p_7"/>
        <point id="8" x="0" y="0.25" z="0.05" name="p_8"/>
        <point id="9" x="0.5" y="0.25" z="0.05" name="p_9"/>
    </points>

    <polylines>
    <polyline id="0" name="top">
		<pnt>4</pnt>
		<pnt>9</pnt>
	</polyline>
	<polyline id="1" name="left">
		<pnt>0</pnt>
		<pnt>5</pnt>
	</polyline>
	<polyline id="2" name="right">
		<pnt>1</pnt>
		<pnt>6</pnt>
	</polyline>
    </polylines>

</OpenGeoSysGLI>
