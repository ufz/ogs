<?xml version="1.0" encoding="utf-8"?>
<OpenGeoSysGLI>
    <name>geometry</name>
    <points>
         <point id="0"  x="0." y="0."  z="0." name = "origin"/>
         <point id="1"  x="1.0" y="0.0"  z="0." />
         <point id="2"  x="1.0" y="1.0"  z="0." />
         <point id="3"  x="0." y="1.0"  z="0." />
    </points>

    <polylines>
        <polyline id="1" name="bottom">
            <pnt>0</pnt>
            <pnt>1</pnt>
        </polyline>
        <polyline id="2" name="top">
            <pnt>2</pnt>
            <pnt>3</pnt>
        </polyline>
    </polylines>
</OpenGeoSysGLI>
