<?xml version="1.0" encoding="utf-8"?>
<OpenGeoSysGLI>
	<name>single_joint</name>
	<points>
		<point id="0" x="0" y="0" z="0" name="POINT0"/>
		<point id="1" x="0.050000000000000003" y="0" z="0"/>
		<point id="2" x="0.050000000000000003" y="0.10000000000000001" z="0"/>
		<point id="3" x="0" y="0.10000000000000001" z="0"/>
		<point id="4" x="0" y="0.050000000000000003" z="0"/>
		<point id="5" x="0.050000000000000003" y="0.050000000000000003" z="0"/>
		<point id="6" x="0.025000000000000001" y="0" z="0"/>
		<point id="7" x="0.025000000000000001" y="0.10000000000000001" z="0"/>
	</points>
	<polylines>
		<polyline id="0" name="PLY_BOTTOM">
			<pnt>0</pnt>
			<pnt>1</pnt>
		</polyline>
		<polyline id="1" name="PLY_TOP">
			<pnt>3</pnt>
			<pnt>2</pnt>
		</polyline>
		<polyline id="1" name="PLY_TOP_LEFT">
			<pnt>3</pnt>
			<pnt>7</pnt>
		</polyline>
		<polyline id="2" name="PLY_LEFT">
			<pnt>0</pnt>
			<pnt>3</pnt>
		</polyline>
		<polyline id="2" name="PLY_RIGHT">
			<pnt>1</pnt>
			<pnt>2</pnt>
		</polyline>
		<polyline id="3" name="PLY_V_PROF">
			<pnt>6</pnt>
			<pnt>7</pnt>
		</polyline>
	</polylines>
</OpenGeoSysGLI>
