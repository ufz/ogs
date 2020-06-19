<OpenGeoSysGLI>
    <name>geometry</name>
    <points>
        <point id="0" x="0.0" y="0" z="0.0" />
        <point id="1" x="0.0" y="0" z="1.0" />
        <point id="2" x="1.0" y="0" z="0.0" />
        <point id="3" x="1.0" y="0" z="1.0" />
        <point id="4" x="0.9800665778412416" y="0.19866933079506122" z="0.0" />
        <point id="5" x="0.9800665778412416" y="0.19866933079506122" z="1.0" />
    </points>
    <surfaces>
        <surface id="0" name="bottom">
            <element p1="0" p2="4" p3="2" />
        </surface>
        <surface id="1" name="top">
            <element p1="1" p2="3" p3="5" />
        </surface>
        <surface id="2" name="outer">
            <element p1="2" p2="4" p3="5" />
            <element p1="5" p2="3" p3="2" />
        </surface>
        <surface id="3" name="upper">
            <element p1="0" p2="4" p3="1" />
            <element p1="1" p2="5" p3="4" />
        </surface>
        <surface id="4" name="lower">
            <element p1="0" p2="2" p3="1" />
            <element p1="1" p2="3" p3="2" />
        </surface>
    </surfaces>
</OpenGeoSysGLI>
