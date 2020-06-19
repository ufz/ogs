<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>geometry</name>
    <points>
        <point id="0" x="0" y="0" z="0" name="bottom_left"/>
        <point id="1" x="100" y="0" z="0"/>
        <point id="2" x="100" y="100" z="0"/>
        <point id="3" x="0" y="100" z="0"/>
        <point id="4" x="50" y="50" z="0" name="po0"/>
        <point id="5" x="40" y="60" z="0" name="po1"/>
        <point id="6" x="45" y="60" z="0" name="po2"/>
        <point id="7" x="50" y="60" z="0" name="po3"/>
        <point id="8" x="55" y="60" z="0" name="po4"/>
        <point id="9" x="60" y="60" z="0" name="po5"/>
        <point id="10" x="40" y="55" z="0" name="po6"/>
        <point id="11" x="45" y="55" z="0" name="po7"/>
        <point id="12" x="50" y="55" z="0" name="po8"/>
        <point id="13" x="55" y="55" z="0" name="po9"/>
        <point id="14" x="60" y="55" z="0" name="po10"/>
        <point id="15" x="40" y="50" z="0" name="po11"/>
        <point id="16" x="45" y="50" z="0" name="po12"/>
        <point id="17" x="55" y="50" z="0" name="po13"/>
        <point id="18" x="60" y="50" z="0" name="po14"/>
        <point id="19" x="40" y="45" z="0" name="po15"/>
        <point id="20" x="45" y="45" z="0" name="po16"/>
        <point id="21" x="50" y="45" z="0" name="po17"/>
        <point id="22" x="55" y="45" z="0" name="po18"/>
        <point id="23" x="60" y="45" z="0" name="po19"/>
        <point id="24" x="40" y="40" z="0" name="po20"/>
        <point id="25" x="45" y="40" z="0" name="po21"/>
        <point id="26" x="50" y="40" z="0" name="po22"/>
        <point id="27" x="55" y="40" z="0" name="po23"/>		
        <point id="28" x="60" y="40" z="0" name="po24"/>		
    </points>

    <polylines>
        <polyline id="0" name="bottom">
            <pnt>0</pnt>
            <pnt>1</pnt>
        </polyline>
        <polyline id="1" name="top">
            <pnt>3</pnt>
            <pnt>2</pnt>
        </polyline>
        <polyline id="2" name="right">
            <pnt>1</pnt>
            <pnt>2</pnt>
        </polyline>
        <polyline id="3" name="left">
            <pnt>0</pnt>
            <pnt>3</pnt>
        </polyline>
    </polylines>
</OpenGeoSysGLI>