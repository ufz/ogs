<?xml version="1.0" encoding="utf-8"?>
<OpenGeoSysGLI>
    <name>single_fracture_3D</name>
    <points>
        <point id="0" x="0" y="0" z="0"/>
        <point id="1" x="0" y="1" z="0"/>
        <point id="2" x="25" y="1" z="0"/>
        <point id="3" x="25" y="0" z="0"/>
        <point id="4" x="0" y="0.5" z="0"/>
        <point id="5" x="25" y="0.5" z="0"/>

        <point id="6" x="0" y="0" z="1"/>
        <point id="7" x="0" y="1" z="1"/>
        <point id="8" x="25" y="1" z="1"/>
        <point id="9" x="25" y="0" z="1"/>
        <point id="10" x="0" y="0.5" z="1"/>
        <point id="11" x="25" y="0.5" z="1"/>
    </points>
    <polylines>
        <polyline id="0" name="inlet">
            <pnt>4</pnt>
            <pnt>10</pnt>
        </polyline>
        <polyline id="1" name="outlet">
            <pnt>5</pnt>
            <pnt>11</pnt>
        </polyline>
        <polyline id="0" name="fracture_back">
            <pnt>4</pnt>
            <pnt>5</pnt>
        </polyline>
        <polyline id="0" name="fracture_front">
            <pnt>10</pnt>
            <pnt>11</pnt>
        </polyline>
    </polylines>
    <surfaces>
        <surface id="0" name="left">
            <element p1="0" p2="1" p3="6"/>
            <element p1="6" p2="1" p3="7"/>
        </surface>
        <surface id="0" name="right">
            <element p1="9" p2="3" p3="2"/>
            <element p1="9" p2="2" p3="8"/>
        </surface>
        <surface id="0" name="front">
            <element p1="6" p2="9" p3="8"/>
            <element p1="6" p2="8" p3="7"/>
        </surface>
        <surface id="0" name="back">
            <element p1="0" p2="3" p3="2"/>
            <element p1="0" p2="2" p3="1"/>
        </surface>
        <surface id="0" name="bottom">
            <element p1="6" p2="9" p3="3"/>
            <element p1="6" p2="3" p3="0"/>
        </surface>
        <surface id="0" name="top">
            <element p1="7" p2="8" p3="2"/>
            <element p1="7" p2="2" p3="1"/>
        </surface>
    </surfaces>
</OpenGeoSysGLI>
