<?xml version="1.0" encoding="utf-8"?>
<OpenGeoSysGLI>
  <name>m3_3Dshearz_rot</name>
  <points>
    <point id="0" x="0" y="0" z="0"/>
    <point id="1" x="4.99999918782138" y="7.07106838616246" z="4.99999999999993"/>
    <point id="2" x=".999999999999987" y="-1.41421344751369" z="1.00000016243572"/>
    <point id="3" x="8.07106838616245" y="-1.41421344751369" z="-6.07106707513271"/>
    <point id="4" x="13.0710675739838" y="5.65685493864878" z="-1.07106707513278"/>
    <point id="5" x="5.99999918782136" y="5.65685493864878" z="6.00000016243566"/>
    <point id="6" x="7.07106838616246" y="0" z="-7.07106723756844"/>
    <point id="7" x="12.0710675739838" y="7.07106838616246" z="-2.07106723756851"/>
  </points>

  <surfaces>
    <surface id="1" name="SURF_YMIN">
      <element p1="6" p2="3" p3="2"/>
      <element p1="0" p2="6" p3="2"/>
    </surface>
    <surface id="2" name="SURF_YMAX">
      <element p1="7" p2="4" p3="5"/>
      <element p1="1" p2="7" p3="5"/>
    </surface>
    <surface id="3" name="SURF_XMAX">
      <element p1="6" p2="7" p3="4"/>
      <element p1="3" p2="6" p3="4"/>
    </surface>
    <surface id="4" name="SURF_XMIN">
      <element p1="5" p2="2" p3="0"/>
      <element p1="1" p2="5" p3="0"/>
    </surface>
    <surface id="5" name="SURF_ZMIN">
      <element p1="6" p2="7" p3="0"/>
      <element p1="7" p2="1" p3="0"/>
    </surface>
    <surface id="6" name="SURF_ZMAX">
      <element p1="4" p2="5" p3="3"/>
      <element p1="5" p2="2" p3="3"/>
    </surface>
  </surfaces>

</OpenGeoSysGLI>
