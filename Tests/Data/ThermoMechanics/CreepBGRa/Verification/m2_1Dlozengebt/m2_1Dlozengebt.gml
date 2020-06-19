<?xml version="1.0" encoding="utf-8"?>
<OpenGeoSysGLI>
  <name>m2_1Dlozengebt</name>
  <points>
    <point id="0" x="0" y="-1" z="-0.1"/>
    <point id="1" x="0" y="1" z="-0.1"/>
    <point id="2" x="0" y="1" z="0.1"/>
    <point id="3" x="0" y="-1" z="0.1"/>
    <point id="4" x="-1" y="0" z="-0.1"/>
    <point id="5" x="1" y="0" z="-0.1"/>
    <point id="6" x="1" y="0" z="0.1"/>
    <point id="7" x="-1" y="0" z="0.1"/>
    <point id="8" x="-1" y="-1" z="0"/>
    <point id="9" x="1" y="-1" z="0"/>
    <point id="10" x="1" y="1" z="0"/>
    <point id="11" x="-1" y="1" z="0"/>
  </points>
  <polylines>
    <polyline id="1" name="POINTSxzero">
      <pnt>0</pnt>
      <pnt>1</pnt>
      <pnt>2</pnt>
      <pnt>3</pnt>
      <pnt>0</pnt>
    </polyline>
    <polyline id="2" name="POINTSyzero">
      <pnt>4</pnt>
      <pnt>5</pnt>
      <pnt>6</pnt>
      <pnt>7</pnt>
      <pnt>4</pnt>
    </polyline>
  </polylines>
  <surfaces>
    <!-- <surface id="1" name="SURFACE1">
      <element p1="3" p2="0" p3="1"/>
      <element p1="2" p2="3" p3="1"/>
    </surface>
    <surface id="2" name="SURFACE2">
      <element p1="7" p2="4" p3="5"/>
      <element p1="6" p2="7" p3="5"/>
    </surface> -->
    <surface id="3" name="SURFACE3">
      <element p1="9" p2="10" p3="8"/>
      <element p1="10" p2="11" p3="8"/>
    </surface>
    <surface id="4" name="SURFACE4">
      <element p1="7" p2="3" p3="0"/>
      <element p1="4" p2="7" p3="0"/>
    </surface>
    <surface id="5" name="SURFACE5">
      <element p1="4" p2="1" p3="2"/>
      <element p1="7" p2="4" p3="2"/>
    </surface>
    <surface id="6" name="SURFACE6">
      <element p1="1" p2="5" p3="2"/>
      <element p1="5" p2="6" p3="2"/>
    </surface>
    <surface id="7" name="SURFACE7">
      <element p1="0" p2="5" p3="3"/>
      <element p1="5" p2="6" p3="3"/>
    </surface>
    <surface id="8" name="SURFACE8">
      <element p1="7" p2="2" p3="6"/>
      <element p1="3" p2="7" p3="6"/>
    </surface>
    <surface id="9" name="SURFACE9">
      <element p1="1" p2="5" p3="0"/>
      <element p1="4" p2="1" p3="0"/>
    </surface>
  </surfaces>
</OpenGeoSysGLI>
