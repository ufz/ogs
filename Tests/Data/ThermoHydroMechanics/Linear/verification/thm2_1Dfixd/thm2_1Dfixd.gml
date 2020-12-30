<?xml version="1.0" encoding="utf-8"?>
<OpenGeoSysGLI>
	<name>thm2_1Dfixed_geom</name>
	<points>
		<point id="0" x="0" y="1.0" z="0" name="POINT0"/>
		<point id="1" x="0" y="0" z="0" name="POINT1"/>
		<point id="2" x="0" y="0" z="1.0" name="POINT2"/>
		<point id="3" x="0" y="1.0" z="1.0" name="POINT3"/>
		<point id="4" x="10" y="1.0" z="0" name="POINT4"/>
		<point id="5" x="10" y="0" z="0" name="POINT5"/>
		<point id="6" x="10" y="0" z="1.0" name="POINT6"/>
		<point id="7" x="10" y="1.0" z="1.0" name="POINT7"/>
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
		<polyline id="2" name="POLYLINE3">
			<pnt>1</pnt>
			<pnt>5</pnt>
			<pnt>6</pnt>
			<pnt>2</pnt>
			<pnt>1</pnt>
		</polyline>
		<polyline id="3" name="POLYLINE4">
			<pnt>0</pnt>
			<pnt>4</pnt>
			<pnt>7</pnt>
			<pnt>3</pnt>
			<pnt>0</pnt>
		</polyline>
		<polyline id="4" name="POLYLINE5">
			<pnt>1</pnt>
			<pnt>5</pnt>
			<pnt>4</pnt>
			<pnt>0</pnt>
			<pnt>1</pnt>
		</polyline>
		<polyline id="5" name="POLYLINE6">
			<pnt>2</pnt>
			<pnt>6</pnt>
			<pnt>7</pnt>
			<pnt>3</pnt>
			<pnt>2</pnt>
		</polyline>
		<polyline id="6" name="POLYLINE7">
			<pnt>1</pnt>
			<pnt>0</pnt>
			<pnt>3</pnt>
			<pnt>2</pnt>
			<pnt>1</pnt>
		</polyline>
		<polyline id="7" name="POLYLINE8">
			<pnt>5</pnt>
			<pnt>4</pnt>
			<pnt>7</pnt>
			<pnt>6</pnt>
			<pnt>5</pnt>
		</polyline>
	</polylines>
	<surfaces>
		<surface id="0" name="SURF_X0">
			<element p1="0" p2="2" p3="1"/>
			<element p1="2" p2="0" p3="3"/>
		</surface>
		<surface id="1" name="SURF_X10">
			<element p1="4" p2="6" p3="5"/>
			<element p1="6" p2="4" p3="7"/>
		</surface>
		<surface id="2" name="SURF_Y0">
			<element p1="1" p2="6" p3="5"/>
			<element p1="6" p2="1" p3="2"/>
		</surface>
		<surface id="3" name="SURF_Y1">
			<element p1="0" p2="7" p3="4"/>
			<element p1="7" p2="0" p3="3"/>
		</surface>
		<surface id="4" name="SURF_Z0">
			<element p1="1" p2="4" p3="5"/>
			<element p1="4" p2="1" p3="0"/>
		</surface>
		<surface id="5" name="SURF_Z1">
			<element p1="2" p2="7" p3="6"/>
			<element p1="7" p2="2" p3="3"/>
		</surface>
	</surfaces>
</OpenGeoSysGLI>
