<?xml version="1.0" encoding="utf-8"?>
<OpenGeoSysGLI>
    <name>excavation_niches2</name>
    <points>
        <point id="0" x="-10" y="0" z="0"/>
        <point id="1" x="10" y="0" z="0"/>
        <point id="2" x="10" y="10" z="0"/>
        <point id="3" x="-10" y="10" z="0"/>
        <point id="4" x="-7.5" y="0" z="0"/>
        <point id="5" x="-7.5" y="5" z="0"/>
        <point id="6" x="-5" y="0" z="0"/>
        <point id="7" x="-5" y="5" z="0"/>
        <point id="8" x="-2.5" y="5" z="0"/>
        <point id="9" x="-2.5" y="0" z="0"/>
        <point id="10" x="2.5" y="0" z="0"/>
        <point id="11" x="2.5" y="5" z="0"/>
        <point id="12" x="7.5" y="5" z="0"/>
        <point id="13" x="7.5" y="0" z="0"/>
    </points>
    <polylines>
        <polyline id="1" name="boundary">
            <pnt>0</pnt>
            <pnt>3</pnt>
            <pnt>2</pnt>
            <pnt>1</pnt>
            <pnt>13</pnt>
            <pnt>12</pnt>
            <pnt>11</pnt>
            <pnt>10</pnt>
            <pnt>9</pnt>
            <pnt>8</pnt>
            <pnt>7</pnt>
            <pnt>5</pnt>
            <pnt>4</pnt>
            <pnt>0</pnt>
        </polyline>
        <polyline id="2" name="niche_1">
            <pnt>4</pnt>
            <pnt>5</pnt>
            <pnt>7</pnt>
            <pnt>6</pnt>
            <pnt>4</pnt>
        </polyline>
        <polyline id="3" name="niche_3">
            <pnt>6</pnt>
            <pnt>7</pnt>
            <pnt>8</pnt>
            <pnt>9</pnt>
            <pnt>6</pnt>
        </polyline>
        <polyline id="4" name="niche_2">
            <pnt>10</pnt>
            <pnt>11</pnt>
            <pnt>12</pnt>
            <pnt>13</pnt>
            <pnt>10</pnt>
        </polyline>
        <polyline id="5" name="left">
            <pnt>0</pnt>
            <pnt>3</pnt>
        </polyline>
        <polyline id="6" name="right">
            <pnt>1</pnt>
            <pnt>2</pnt>
        </polyline>
        <polyline id="7" name="back">
            <pnt>3</pnt>
            <pnt>2</pnt>
        </polyline>
    </polylines>
</OpenGeoSysGLI>
