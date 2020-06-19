<?xml version="1.0" encoding="ISO-8859-1"?>
<!DOCTYPE OGS-GML-DOM>
<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
 <name>elder</name>
 <points>
  <point z="75" id="0" x="-150" name="POINT0" y="-1.17"/>
  <point z="75" id="1" x="-150" name="POINT1" y="1.17"/>
  <point z="-75" id="2" x="-150" name="POINT2" y="-1.17"/>
  <point z="-75" id="3" x="-150" name="POINT1.17" y="1.17"/>
  <point z="-75" id="4" x="150" name="POINT4" y="-1.17"/>
  <point z="-75" id="5" x="150" name="POINT5" y="1.17"/>
  <point z="75" id="6" x="150" name="POINT6" y="-1.17"/>
  <point z="75" id="7" x="150" name="POINT7" y="1.17"/>
  <point z="75" id="8" x="0" name="POINT8" y="-1.17"/>
  <point z="75" id="9" x="0" name="POINT9" y="1.17"/>
 </points>
 <polylines>
  <polyline id="0" name="LEFT_TOP_CORNER">
   <pnt>0</pnt>
   <pnt>1</pnt>
  </polyline>
  <polyline id="1" name="TOP_HALF_RIGHT">
   <pnt>8</pnt>
   <pnt>9</pnt>
   <pnt>7</pnt>
   <pnt>6</pnt>
   <pnt>8</pnt>
  </polyline>
  <polyline id="2" name="BOTTOM">
   <pnt>2</pnt>
   <pnt>3</pnt>
   <pnt>5</pnt>
   <pnt>4</pnt>
   <pnt>2</pnt>
  </polyline>
 </polylines>
 <surfaces>
  <surface id="0" name="TOP_HALF_RIGHT_S">
   <element p2="7" p1="8" p3="9"/>
   <element p2="8" p1="7" p3="6"/>
  </surface>
  <surface id="1" name="BOTTOM_S">
   <element p2="5" p1="2" p3="3"/>
   <element p2="2" p1="5" p3="4"/>
  </surface>
 </surfaces>
</OpenGeoSysGLI>
