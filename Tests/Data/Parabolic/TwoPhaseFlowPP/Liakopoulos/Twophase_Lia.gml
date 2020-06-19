<?xml version="1.0" encoding="ISO-8859-1"?>
<!DOCTYPE OGS-GML-DOM>
<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
 <name>Twophase_Lia_Geometry</name>
 <points>
  <point x="0.000000" y="0.000000" z="0.000000" id="0" name="POINT0"/>
  <point x="0.100000" y="0.000000" z="0.000000" id="1"/>
  <point x="0.100000" y="1.000000" z="0.000000" id="2"/>
  <point x="0.000000" y="1.000000" z="0.000000" id="3"/>
 </points>
 <polylines>
  <polyline id="0" name="TOP">
   <pnt>3</pnt>
   <pnt>2</pnt>
  </polyline>
  <polyline id="1" name="BOTTOM">
   <pnt>0</pnt>
   <pnt>1</pnt>
  </polyline>
  <polyline id="2" name="left">
   <pnt>0</pnt>
   <pnt>3</pnt>
  </polyline>
  <polyline id="3" name="right">
   <pnt>2</pnt>
   <pnt>1</pnt>
  </polyline>
 </polylines>
</OpenGeoSysGLI>
