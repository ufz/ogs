<?xml version="1.0" encoding="utf-8"?>
<OpenGeoSysGLI>
	<name>geometry</name>
	<points>
		<point id="0" x="100" y="0" z="0" name="POINT0"/>
		<point id="1" x="100" y="10" z="0" name="POINT1"/>
		<point id="2" x="100" y="10" z="10" name="POINT2"/>
		<point id="3" x="100" y="0" z="10" name="POINT3"/>
		<point id="4" x="0" y="0" z="0" name="POINT4"/>
		<point id="5" x="0" y="10" z="0" name="POINT5"/>
		<point id="6" x="0" y="10" z="10" name="POINT6"/>
		<point id="7" x="0" y="0" z="10" name="POINT7"/>
	</points>
	<polylines>
		<polyline id="0" name="POLYLINE1">
			<pnt>0</pnt>
			<pnt>1</pnt>
			<pnt>2</pnt>
			<pnt>3</pnt>
			<pnt>0</pnt>
		</polyline>
		<polyline id="1" name="POLYLINE2">
			<pnt>4</pnt>
			<pnt>5</pnt>
			<pnt>6</pnt>
			<pnt>7</pnt>
			<pnt>4</pnt>
		</polyline>
	</polylines>
	<surfaces>
		<surface id="0" name="SURFACE1">
			<element p1="0" p2="2" p3="1"/>
			<element p1="2" p2="0" p3="3"/>
		</surface>
		<surface id="1" name="SURFACE2">
			<element p1="4" p2="6" p3="5"/>
			<element p1="6" p2="4" p3="7"/>
		</surface>
	</surfaces>
</OpenGeoSysGLI>
