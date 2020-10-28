<?xml version="1.0" encoding="utf-8"?>
<OpenGeoSysGLI>
  <name>2D_meter</name>
  <points>
    <point id="0" x="0" y="0" z="0"/>
    <point id="1" x=".707106838616246" y="-.707106723756844" z="0"/>
    <point id="2" x="7.77817407618469" y="6.36396166240562" z="0"/>
    <point id="3" x="7.07106723756844" y="7.07106838616246" z="0"/>
  </points>

  <polylines>
    <polyline id="1" name = "bottom">
      <pnt>0</pnt>
      <pnt>1</pnt>
    </polyline>
    <polyline id="2" name = "right">
      <pnt>1</pnt>
      <pnt>2</pnt>
    </polyline>
    <polyline id="3" name = "top">
      <pnt>2</pnt>
      <pnt>3</pnt>
    </polyline>
    <polyline id="4" name = "left">
      <pnt>3</pnt>
      <pnt>0</pnt>
    </polyline>
  </polylines>

</OpenGeoSysGLI>
