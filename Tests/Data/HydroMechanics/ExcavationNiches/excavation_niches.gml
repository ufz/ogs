<?xml version="1.0" encoding="utf-8"?>
<OpenGeoSysGLI>
    <name>example</name>
    <points>
        <point id="0" x="-15" y="0" z="0"/>
        <point id="1" x="-7.3" y="0" z="0"/>
        <point id="2" x="-5" y="0" z="0"/>
        <point id="3" x="5" y="0" z="0"/>
        <point id="4" x="7.3" y="0" z="0"/>
        <point id="5" x="15" y="0" z="0"/>
        <point id="6" x="15" y="15" z="0"/>
        <point id="7" x="-15" y="15" z="0"/>
        <point id="8" x="-7.3" y="11" z="0"/>
        <point id="9" x="-5" y="11" z="0"/>
        <point id="10" x="5" y="11" z="0"/>
        <point id="11" x="7.3" y="11" z="0"/>
    </points>
    <polylines>
        <polyline id="1" name="open_niche">
            <pnt>1</pnt>
            <pnt>8</pnt>
            <pnt>9</pnt>
            <pnt>2</pnt>
            <pnt>1</pnt>
        </polyline>
        <polyline id="2" name="closed_niche">
            <pnt>3</pnt>
            <pnt>10</pnt>
            <pnt>11</pnt>
            <pnt>4</pnt>
            <pnt>3</pnt>
        </polyline>
        <polyline id="3" name="left">
            <pnt>0</pnt>
            <pnt>7</pnt>
        </polyline>
        <polyline id="4" name="back">
            <pnt>7</pnt>
            <pnt>6</pnt>
        </polyline>
        <polyline id="5" name="right">
            <pnt>6</pnt>
            <pnt>5</pnt>
        </polyline>
        <polyline id="6" name="front">
            <pnt>0</pnt>
            <pnt>5</pnt>
        </polyline>
    </polylines>
</OpenGeoSysGLI>
