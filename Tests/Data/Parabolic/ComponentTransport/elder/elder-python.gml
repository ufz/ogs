<?xml version="1.0" encoding="ISO-8859-1"?>
<!-- vim: ft=xml sw=2
-->
<!DOCTYPE OGS-GML-DOM>
<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
  <name>elder</name>
  <points>
    <point id="0" x="-150" y="-1.17" z="-75" />
    <point id="1" x="150" y="-1.17" z="-75" />
    <point id="2" x="150" y="-1.17" z="75" />
    <point id="3" x="-150" y="-1.17" z="75" />
    <point id="4" x="-150" y="1.17" z="-75" />
    <point id="5" x="150" y="1.17" z="-75" />
    <point id="6" x="-150" y="1.17" z="75" />
    <point id="7" x="150" y="1.17" z="75" />
  </points>
  <surfaces>
    <surface id="0" name="whole_domain_boundary">
      <element p1="0" p2="1" p3="3" />
      <element p1="2" p2="3" p3="1" />
      <element p1="0" p2="4" p3="1" />
      <element p1="5" p2="1" p3="4" />
      <element p1="0" p2="3" p3="4" />
      <element p1="6" p2="4" p3="3" />
      <element p1="1" p2="5" p3="2" />
      <element p1="7" p2="2" p3="5" />
      <element p1="4" p2="6" p3="5" />
      <element p1="7" p2="5" p3="6" />
      <element p1="3" p2="2" p3="6" />
      <element p1="7" p2="6" p3="2" />
    </surface>
  </surfaces>
</OpenGeoSysGLI>
