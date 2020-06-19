<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>
<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>square_1x1_geometry</name>
    <points>
        <point id="0" x="0" y="0" z="0" name="center"/>
        <point id="1" x="0.00204357" y="0" z="0"/>
        <point id="2" x="0" y="0.4" z="0"/>
        <point id="3" x="0" y="0.6" z="0"/>
        <point id="4" x="10" y="0" z="0.0000000000000000e+00"/>
        <point id="5" x="9.95185" y="0.980171" z="0.0000000000000000e+00"/>
        <point id="6" x="9.80785" y="1.9509" z="0.0000000000000000e+00"/>
        <point id="7" x="9.5694" y="2.90285" z="0.0000000000000000e+00"/>
        <point id="8" x="9.2388" y="3.82583" z="0.0000000000000000e+00"/>
        <point id="9" x="8.81921" y="4.71397" z="0.0000000000000000e+00"/>
        <point id="10" x="8.3147" y="5.5557" z="0.0000000000000000e+00"/>
        <point id="11" x="7.7301" y="6.34393" z="0.0000000000000000e+00"/>
        <point id="12" x="7.07107" y="7.07107" z="0.0000000000000000e+00"/>
        <point id="13" x="6.34393" y="7.7301" z="0.0000000000000000e+00"/>
        <point id="14" x="5.5557" y="8.3147" z="0.0000000000000000e+00"/>
        <point id="15" x="4.71397" y="8.81921" z="0.0000000000000000e+00"/>
        <point id="16" x="3.82683" y="9.2388" z="0.0000000000000000e+00"/>
        <point id="17" x="2.90285" y="9.5694" z="0.0000000000000000e+00"/>
        <point id="18" x="1.9509" y="9.80785" z="0.0000000000000000e+00"/>
        <point id="19" x="0.980171" y="9.95185" z="0.0000000000000000e+00"/>
        <point id="20" x="0" y="10" z="0.0000000000000000e+00"/>
    </points>
    <polylines>
        <polyline id="0" name="left">
            <pnt>0</pnt>
            <pnt>20</pnt>
        </polyline>
        <polyline id="1" name="bottom">
            <pnt>0</pnt>
            <pnt>4</pnt>
        </polyline>
        <polyline id="2" name="flux_bc">
            <pnt>0</pnt>
            <pnt>1</pnt>
        </polyline>
        <polyline id="3" name="out">
            <pnt>4</pnt>
            <pnt>5</pnt>
            <pnt>6</pnt>
            <pnt>7</pnt>
            <pnt>8</pnt>
            <pnt>9</pnt>
            <pnt>10</pnt>
            <pnt>11</pnt>
            <pnt>12</pnt>
            <pnt>13</pnt>
            <pnt>14</pnt>
            <pnt>15</pnt>
            <pnt>16</pnt>
            <pnt>17</pnt>
            <pnt>18</pnt>
            <pnt>19</pnt>
            <pnt>20</pnt>
        </polyline>
    </polylines>
</OpenGeoSysGLI>
