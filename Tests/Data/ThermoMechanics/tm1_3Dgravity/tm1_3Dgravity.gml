<?xml version="1.0" encoding="utf-8"?>
<OpenGeoSysGLI>
	<name>geometry</name>
	<points>
		<point id="0" x="0" y="0" z="30" name="POINT0"/>
		<point id="1" x="60" y="0" z="30" name="POINT1"/>
		<point id="2" x="60" y="30" z="30" name="POINT2"/>
		<point id="3" x="0" y="30" z="30" name="POINT3"/>
		<point id="4" x="0" y="0" z="0" name="POINT4"/>
		<point id="5" x="60" y="0" z="0" name="POINT5"/>
		<point id="6" x="60" y="30" z="0" name="POINT6"/>
		<point id="7" x="0" y="30" z="0" name="POINT7"/>
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
			<pnt>4</pnt>
			<pnt>5</pnt>
			<pnt>6</pnt>
			<pnt>7</pnt>
			<pnt>4</pnt>
		</polyline>
		<polyline id="3" name="POLYLINE4">
			<pnt>4</pnt>
			<pnt>5</pnt>
			<pnt>1</pnt>
			<pnt>0</pnt>
			<pnt>4</pnt>
		</polyline>
		<polyline id="4" name="POLYLINE5">
			<pnt>7</pnt>
			<pnt>6</pnt>
			<pnt>2</pnt>
			<pnt>3</pnt>
			<pnt>7</pnt>
		</polyline>
		<polyline id="5" name="POLYLINE6">
			<pnt>4</pnt>
			<pnt>7</pnt>
			<pnt>3</pnt>
			<pnt>0</pnt>
			<pnt>4</pnt>
		</polyline>
		<polyline id="6" name="POLYLINE7">
			<pnt>5</pnt>
			<pnt>6</pnt>
			<pnt>2</pnt>
			<pnt>1</pnt>
			<pnt>5</pnt>
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
		<surface id="2" name="SURFACE3">
			<element p1="4" p2="6" p3="5"/>
			<element p1="6" p2="4" p3="7"/>
		</surface>
		<surface id="3" name="SURFACE4">
			<element p1="4" p2="1" p3="5"/>
			<element p1="1" p2="4" p3="0"/>
		</surface>
		<surface id="4" name="SURFACE5">
			<element p1="7" p2="2" p3="6"/>
			<element p1="2" p2="7" p3="3"/>
		</surface>
		<surface id="5" name="SURFACE6">
			<element p1="4" p2="3" p3="7"/>
			<element p1="3" p2="4" p3="0"/>
		</surface>
		<surface id="6" name="SURFACE7">
			<element p1="5" p2="2" p3="6"/>
			<element p1="2" p2="5" p3="1"/>
		</surface>
	</surfaces>
</OpenGeoSysGLI>
