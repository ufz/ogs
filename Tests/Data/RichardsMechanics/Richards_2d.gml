<?xml version="1.0" encoding="ISO-8859-1"?>
<!DOCTYPE OGS-GML-DOM>
<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
 <name>Richards_2d_geometry</name>
 <points>
  <point x="0.000000" y="2.000000" z="0.000000" id="0"/>
  <point x="0.020000" y="2.000000" z="0.000000" id="1"/>
  <point x="0.040000" y="2.000000" z="0.000000" id="2"/>
  <point x="0.000000" y="0.000000" z="0.000000" id="3"/>
  <point x="0.020000" y="0.000000" z="0.000000" id="4"/>
  <point x="0.040000" y="0.000000" z="0.000000" id="5"/>
  <point x="0.02" y="0.08" z="0" id="6" name="OBSERVATION_POINT"/>
 </points>
 <polylines>
  <polyline id="0" name="TOP">
   <pnt>0</pnt>
   <pnt>1</pnt>
   <pnt>2</pnt>
  </polyline>
  <polyline id="1" name="BOTTOM">
   <pnt>3</pnt>
   <pnt>4</pnt>
   <pnt>5</pnt>
  </polyline>
  <polyline id="2" name="LEFT">
   <pnt>0</pnt>
   <pnt>3</pnt>
  </polyline>
  <polyline id="3" name="RIGHT">
   <pnt>2</pnt>
   <pnt>5</pnt>
  </polyline>
 </polylines>
</OpenGeoSysGLI>
