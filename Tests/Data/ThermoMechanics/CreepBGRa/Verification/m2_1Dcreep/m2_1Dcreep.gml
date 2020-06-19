<?xml version="1.0" encoding="utf-8"?>
<OpenGeoSysGLI>
  <name>m2_1Dcreep</name>
  <points>
    <point id="0" x="0" y="0" z="0"/>
    <point id="1" x="1" y="0" z="0"/>
    <point id="2" x="1" y="1" z="0"/>
    <point id="3" x="0" y="1" z="0"/>
    <point id="4" x="0" y="0" z="0"/>
    <point id="5" x="1" y="0" z="0"/>
    <point id="6" x="1" y="0" z="1"/>
    <point id="7" x="0" y="0" z="1"/>
    <point id="8" x="0" y="0" z="0"/>
    <point id="9" x="0" y="1" z="0"/>
    <point id="10" x="0" y="1" z="1"/>
    <point id="11" x="0" y="0" z="1"/>
    <point id="12" x="0" y="0" z="1"/>
    <point id="13" x="1" y="0" z="1"/>
    <point id="14" x="1" y="1" z="1"/>
    <point id="15" x="0" y="1" z="1"/>
    <point id="16" x="0.5" y="0.5" z="0"/>
    <point id="17" x="0.5" y="0" z="0.5"/>
    <point id="18" x="0" y="0.5" z="0.5"/>
    <point id="19" x="0.5" y="0.5" z="1"/>
  </points>

  <surfaces>
    <surface id="1" name="SURFACE1">
      <element p1="1" p2="16" p3="0"/>
      <element p1="3" p2="16" p3="2"/>
      <element p1="16" p2="3" p3="0"/>
      <element p1="16" p2="1" p3="2"/>
    </surface>
    <surface id="2" name="SURFACE2">
      <element p1="17" p2="5" p3="6"/>
      <element p1="17" p2="7" p3="4"/>
      <element p1="7" p2="17" p3="6"/>
      <element p1="5" p2="17" p3="4"/>
    </surface>
    <surface id="3" name="SURFACE3">
      <element p1="9" p2="18" p3="8"/>
      <element p1="11" p2="18" p3="10"/>
      <element p1="18" p2="11" p3="8"/>
      <element p1="18" p2="9" p3="10"/>
    </surface>
    <surface id="4" name="SURFACE4">
      <element p1="13" p2="19" p3="12"/>
      <element p1="15" p2="19" p3="14"/>
      <element p1="19" p2="15" p3="12"/>
      <element p1="19" p2="13" p3="14"/>
    </surface>
  </surfaces>

</OpenGeoSysGLI>
