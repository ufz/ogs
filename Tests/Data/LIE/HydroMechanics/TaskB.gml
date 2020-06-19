<?xml version="1.0" encoding="utf-8"?>
<OpenGeoSysGLI>
    <name>TaskB</name>
    <points>
        <point id="0" x="-10" y="-10" z="0"/>
        <point id="1" x="10" y="-10" z="0"/>
        <point id="2" x="10" y="10" z="0"/>
        <point id="3" x="-10" y="10" z="0"/>
        <point id="4" x="-4.6323499999999997" y="-10" z="0" name="FRACTURE_BOTTOM"/>
        <point id="5" x="4.6323499999999997" y="10" z="0" name="FRACTURE_TOP"/>
        <point id="6" x="0" y="0" z="0" name="INJECTION"/>
    </points>
    <polylines>
        <polyline id="0" name="PLY_TOP">
            <pnt>3</pnt>
            <pnt>2</pnt>
        </polyline>
        <polyline id="1" name="PLY_BOTTOM">
            <pnt>0</pnt>
            <pnt>1</pnt>
        </polyline>
        <polyline id="2" name="PLY_LEFT">
            <pnt>0</pnt>
            <pnt>3</pnt>
        </polyline>
        <polyline id="3" name="PLY_RIGHT">
            <pnt>1</pnt>
            <pnt>2</pnt>
        </polyline>
    </polylines>
</OpenGeoSysGLI>
