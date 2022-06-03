<?xml version="1.0" encoding="utf-8"?>
<OpenGeoSysGLI>
	<name>geometry</name>
	<points>
		<point id="0" x="0" y="-0.5" z="-0.5"/>
		<point id="1" x="0" y="0.5" z="-0.5"/>
		<point id="2" x="0" y="0.5" z="0.5"/>
		<point id="3" x="0" y="-0.5" z="0.5"/>
		<point id="4" x="-10" y="0.5" z="-0.5"/>
		<point id="5" x="10" y="0.5" z="-0.5"/>
		<point id="6" x="10" y="0.5" z="0.5"/>
		<point id="7" x="-10" y="0.5" z="0.5"/>
		<point id="8" x="-10" y="-0.5" z="-0.5"/>
		<point id="9" x="10" y="-0.5" z="-0.5"/>
		<point id="10" x="10" y="-0.5" z="0.5"/>
		<point id="11" x="-10" y="-0.5" z="0.5"/>
		<point id="12" x="0" y="-0.01" z="-0.01"/>
		<point id="13" x="0" y="0.01" z="-0.01"/>
		<point id="14" x="0" y="0.01" z="0.01"/>
		<point id="15" x="0" y="-0.01" z="0.01"/>
        <point id="16" x="0" y="0" z="0" name="POINT_ORIGIN"/>
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
			<pnt>8</pnt>
			<pnt>9</pnt>
			<pnt>10</pnt>
			<pnt>11</pnt>
			<pnt>8</pnt>
		</polyline>
		<polyline id="3" name="POLYLINE4">
			<pnt>12</pnt>
			<pnt>13</pnt>
			<pnt>14</pnt>
			<pnt>15</pnt>
			<pnt>12</pnt>
		</polyline>
		<polyline id="4" name="POLYLINE5">
			<pnt>11</pnt>
			<pnt>10</pnt>
			<pnt>6</pnt>
			<pnt>7</pnt>
			<pnt>11</pnt>
		</polyline>
		<polyline id="5" name="POLYLINE6">
			<pnt>8</pnt>
			<pnt>9</pnt>
			<pnt>5</pnt>
			<pnt>4</pnt>
			<pnt>8</pnt>
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
			<element p1="8" p2="10" p3="9"/>
			<element p1="10" p2="8" p3="11"/>
		</surface>
		<surface id="3" name="SURFACE4">
			<element p1="12" p2="14" p3="13"/>
			<element p1="14" p2="12" p3="15"/>
		</surface>
		<surface id="4" name="SURFACE5">
			<element p1="11" p2="6" p3="10"/>
			<element p1="6" p2="11" p3="7"/>
		</surface>
		<surface id="5" name="SURFACE6">
			<element p1="8" p2="5" p3="9"/>
			<element p1="5" p2="8" p3="4"/>
		</surface>
	</surfaces>
</OpenGeoSysGLI>
