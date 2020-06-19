<?xml version="1.0" encoding="utf-8"?>
<OpenGeoSysGLI>
  <name>plate_ym45</name>
  <points>
    <point id="0" x="0" y="0" z="0"/>
    <point id="1" x="7.07106838616247" y="0" z="7.07106723756844"/>
    <point id="2" x="7.07106838616247" y="10" z="7.07106723756844"/>
    <point id="3" x="0" y="10" z="0"/>
    <point id="4" x="-1.41421344751369" y="0" z="1.41421367723249"/>
    <point id="5" x="5.65685493864878" y="0" z="8.48528091480093"/>
    <point id="6" x="5.65685493864878" y="10" z="8.48528091480093"/>
    <point id="7" x="-1.41421344751369" y="10" z="1.41421367723249"/>
  </points>

  <surfaces>
    <surface id="1" name="SURFACE1">
      <element p1="2" p2="3" p3="1"/>
      <element p1="3" p2="0" p3="1"/>
    </surface>
    <surface id="2" name="SURFACE2">
      <element p1="5" p2="6" p3="4"/>
      <element p1="6" p2="7" p3="4"/>
    </surface>
    <surface id="3" name="SURFACE3">
      <element p1="4" p2="7" p3="3"/>
      <element p1="0" p2="4" p3="3"/>
    </surface>
    <surface id="4" name="SURFACE4">
      <element p1="5" p2="0" p3="1"/>
      <element p1="4" p2="0" p3="5"/>
    </surface>
    <surface id="5" name="SURFACE5">
      <element p1="5" p2="1" p3="2"/>
      <element p1="6" p2="5" p3="2"/>
    </surface>
    <surface id="6" name="SURFACE6">
      <element p1="3" p2="6" p3="2"/>
      <element p1="7" p2="6" p3="3"/>
    </surface>
  </surfaces>

</OpenGeoSysGLI>
