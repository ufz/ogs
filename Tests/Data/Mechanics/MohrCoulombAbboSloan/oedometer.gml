<?xml version="1.0" encoding="utf-8"?>
<OpenGeoSysGLI>
  <name>2D_meter</name>
  <points>
    <point id="0" x="0" y=".1" z="0"/>
    <point id="1" x="0" y="0" z="0"/>
    <point id="2" x=".025" y="0" z="0"/>
    <point id="3" x=".025" y=".1" z="0"/>
  </points>

  <polylines>
    <polyline id="1" name = "left">
      <pnt>0</pnt>
      <pnt>1</pnt>
    </polyline>
    <polyline id="2" name = "right">
      <pnt>3</pnt>
      <pnt>2</pnt>
    </polyline>
    <polyline id="3" name = "bottom">
      <pnt>1</pnt>
      <pnt>2</pnt>
    </polyline>
    <polyline id="4" name = "top">
      <pnt>0</pnt>
      <pnt>3</pnt>
    </polyline>
  </polylines>

</OpenGeoSysGLI>
