<?xml version="1.0" encoding="ISO-8859-1"?>
<!DOCTYPE OGS-GML-DOM>
<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
 <name>H2_H20_2D</name>
 <points>
  <point x="0.000000" y="0.000000" z="0.000000" id="0"/>
  <point x="200.000000" y="0.000000" z="0.000000" id="1"/>
  <point x="200.000000" y="20.000000" z="0.000000" id="2"/>
  <point x="0.000000" y="20.000000" z="0.000000" id="3"/>
 </points>
 <polylines>
  <polyline id="0" name="PLY_LEFT">
   <pnt>3</pnt>
   <pnt>0</pnt>
  </polyline>
  <polyline id="1" name="PLY_RIGHT">
   <pnt>1</pnt>
   <pnt>2</pnt>
  </polyline>
  <polyline id="2" name="PLY_ALL">
   <pnt>0</pnt>
   <pnt>1</pnt>
   <pnt>2</pnt>
   <pnt>3</pnt>
   <pnt>0</pnt>
  </polyline>
 </polylines>
</OpenGeoSysGLI>
