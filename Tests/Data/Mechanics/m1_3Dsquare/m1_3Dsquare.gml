<?xml version="1.0" encoding="utf-8"?>
<OpenGeoSysGLI>
	<name>geometry</name>
	<points>
		<point id="0" x="0" y="-10" z="-2"/>
		<point id="1" x="0" y="10" z="-2"/>
		<point id="2" x="0" y="10" z="0"/>
		<point id="3" x="0" y="-10" z="0"/>
		<point id="4" x="-10" y="0" z="-2"/>
		<point id="5" x="10" y="0" z="-2"/>
		<point id="6" x="10" y="0" z="0"/>
		<point id="7" x="-10" y="0" z="0"/>
        <point id="8" x="0" y="-0.01" z="-0.01"/>
        <point id="9" x="0" y="0.01" z="-0.01"/>
        <point id="10" x="0" y="0.01" z="0.01"/>
        <point id="11" x="0" y="-0.01" z="0.01"/>
        <point id="12" x="-10" y="-10" z="-2"/>
        <point id="13" x="10" y="-10" z="-2"/>
        <point id="14" x="10" y="10" z="-2"/>
        <point id="15" x="-10" y="10" z="-2"/>
		<point id="16" x="0" y="0" z="0" name="POINT_ORIGIN"/>        
	</points>
	<surfaces>
		<surface id="0" name="SURFACE1">
			<element p1="0" p2="2" p3="1"/>
			<element p1="2" p2="0" p3="3"/>
		</surface>
		<surface id="1" name="SURFACE2">
			<element p1="4" p2="6" p3="5"/>
			<element p1="6" p2="4" p3="7"/>
		</surface>
        <surface id="3" name="SURFACE4">
                <element p1="12" p2="14" p3="13"/>
                <element p1="14" p2="12" p3="15"/>
        </surface>
	</surfaces>
</OpenGeoSysGLI>
