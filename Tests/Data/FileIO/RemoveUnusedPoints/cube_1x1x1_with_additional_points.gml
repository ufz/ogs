<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>cube_1x1x1_geometry</name>
    <points>
        <point id="0" x="0" y="0" z="0" name="left_front_bottom"/>
        <point id="1" x="0" y="0" z="1" name="left_front_top"/>
        <point id="2" x="0" y="1" z="1" name="left_back_top"/>
        <point id="3" x="0" y="1" z="0" name="left_back_bottom"/>
        <point id="4" x="0.1" y="0" z="0" name="invalid0"/>
        <point id="5" x="1" y="0" z="0" name="right_front_bottom"/>
        <point id="6" x="1" y="0" z="1" name="right_front_top"/>
        <point id="7" x="0.2" y="0" z="0" name="invalid1"/>
        <point id="8" x="1" y="1" z="1" name="right_back_top"/>
        <point id="9" x="1" y="1" z="0" name="right_back_bottom"/>
        <point id="10" x="0" y="0" z="0" name="invalid2"/>
        <point id="11" x="0" y="0" z="0" name="invalid3"/>
        <point id="12" x="0.1" y="0" z="0" name="invalid4"/>
    </points>

    <polylines>
        <polyline id="0" name="vertical_polyline_left_front">
            <pnt>0</pnt>
            <pnt>1</pnt>
        </polyline>
        <polyline id="1" name="vertical_polyline_right_front">
            <pnt>5</pnt>
            <pnt>6</pnt>
        </polyline>
        <polyline id="2" name="vertical_polyline_left_back">
            <pnt>3</pnt>
            <pnt>2</pnt>
        </polyline>
        <polyline id="3" name="vertical_polyline_right_back">
            <pnt>9</pnt>
            <pnt>8</pnt>
        </polyline>
        <polyline id="4" name="horizontal_polyline_front_bottom">
            <pnt>0</pnt>
            <pnt>5</pnt>
        </polyline>
        <polyline id="5" name="horizontal_polyline_front_top">
            <pnt>1</pnt>
            <pnt>6</pnt>
        </polyline>
        <polyline id="6" name="horizontal_polyline_back_bottom">
            <pnt>3</pnt>
            <pnt>9</pnt>
        </polyline>
        <polyline id="7" name="horizontal_polyline_back_top">
            <pnt>2</pnt>
            <pnt>8</pnt>
        </polyline>
        <polyline id="8" name="horizontal_polyline_left_bottom">
            <pnt>0</pnt>
            <pnt>3</pnt>
        </polyline>
        <polyline id="9" name="horizontal_polyline_left_top">
            <pnt>1</pnt>
            <pnt>2</pnt>
        </polyline>
        <polyline id="10" name="horizontal_polyline_right_top">
            <pnt>6</pnt>
            <pnt>8</pnt>
        </polyline>
        <polyline id="11" name="horizontal_polyline_right_bottom">
            <pnt>5</pnt>
            <pnt>9</pnt>
        </polyline>
    </polylines>

    <surfaces>
        <surface id="0" name="left"><!-- x=0 -->
            <element p1="0" p2="1" p3="2"/>
            <element p1="0" p2="3" p3="2"/>
        </surface>
        <surface id="1" name="right"><!-- x=1 -->
            <element p1="5" p2="8" p3="6"/>
            <element p1="5" p2="8" p3="9"/>
        </surface>
        <surface id="2" name="top"><!-- z=1 -->
            <element p1="1" p2="2" p3="6"/>
            <element p1="6" p2="2" p3="8"/>
        </surface>
        <surface id="3" name="bottom"><!-- z=0 -->
            <element p1="0" p2="3" p3="5"/>
            <element p1="5" p2="3" p3="9"/>
        </surface>
        <surface id="4" name="front"><!-- y=0 -->
            <element p1="0" p2="1" p3="5"/>
            <element p1="5" p2="1" p3="6"/>
        </surface>
        <surface id="5" name="back"><!-- y=1 -->
            <element p1="2" p2="3" p3="8"/>
            <element p1="8" p2="3" p3="9"/>
        </surface>
    </surfaces>
</OpenGeoSysGLI>
