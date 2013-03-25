<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<!DOCTYPE OGS-GML-DOM>
<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://141.65.34.25/OpenGeoSysCND.xsd" xmlns:ogs="http://www.opengeosys.com">
 <name>TestData</name>
 <points>
  <point x="1.000000" y="1.000000" z="0.000000" id="0"/>
  <point x="1.000000" y="2.000000" z="0.000000" id="1"/>
  <point x="1.000000" y="3.000000" z="0.000000" id="2"/>
  <point x="2.000000" y="1.000000" z="0.000000" id="3"/>
  <point x="2.000000" y="2.000000" z="0.000000" id="4"/>
  <point x="2.000000" y="3.000000" z="0.000000" id="5"/>
  <point x="3.000000" y="1.000000" z="0.000000" id="6"/>
  <point x="3.000000" y="2.000000" z="0.000000" id="7"/>
  <point x="3.000000" y="3.000000" z="0.000000" id="8"/>
 </points>
 <polylines>
  <polyline id="0" name="left">
   <pnt>0</pnt>
   <pnt>1</pnt>
   <pnt>2</pnt>
  </polyline>
  <polyline id="1" name="center">
   <pnt>3</pnt>
   <pnt>4</pnt>
   <pnt>5</pnt>
  </polyline>
  <polyline id="2">
   <pnt>0</pnt>
   <pnt>3</pnt>
  </polyline>
  <polyline id="3">
   <pnt>3</pnt>
   <pnt>6</pnt>
  </polyline>
  <polyline id="4" name="right">
   <pnt>6</pnt>
   <pnt>7</pnt>
   <pnt>8</pnt>
  </polyline>
 </polylines>
 <surfaces>
  <surface id="0">
   <element p1="0" p2="3" p3="1"/>
   <element p1="1" p2="3" p3="4"/>
   <element p1="1" p2="4" p3="2"/>
   <element p1="2" p2="4" p3="5"/>
  </surface>
  <surface id="1">
   <element p1="3" p2="6" p3="8"/>
   <element p1="3" p2="8" p3="5"/>
  </surface>
 </surfaces>
</OpenGeoSysGLI>
