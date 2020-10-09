<?xml version="1.0" encoding="utf-8"?>
<OpenGeoSysGLI>
    <name>five_element_two_fractures</name>
    <points>
        <point id="0" x="0" y="0" z="0"/>
        <point id="1" x="0" y="1" z="0" name="FRACTURE1_LEFT"/>
        <point id="2" x="0" y="2" z="0" name="FRACTURE2_LEFT"/>
        <point id="3" x="0" y="3" z="0"/>
        <point id="4" x="1" y="0" z="0"/>
        <point id="5" x="1" y="1" z="0" name="FRACTURE1_RIGHT"/>
        <point id="6" x="1" y="2" z="0" name="FRACTURE2_RIGHT"/>
        <point id="7" x="1" y="3" z="0"/>
    </points>
    <polylines>
        <polyline id="0" name="PLY_BOTTOM">
            <pnt>0</pnt>
            <pnt>4</pnt>
        </polyline>
        <polyline id="1" name="PLY_TOP">
            <pnt>3</pnt>
            <pnt>7</pnt>
        </polyline>
        <polyline id="2" name="PLY_LEFT">
            <pnt>0</pnt>
            <pnt>3</pnt>
        </polyline>
        <polyline id="3" name="PLY_RIGHT">
            <pnt>4</pnt>
            <pnt>7</pnt>
        </polyline>
    </polylines>
</OpenGeoSysGLI>
