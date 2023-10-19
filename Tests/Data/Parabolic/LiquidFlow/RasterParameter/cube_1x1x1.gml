<?xml version="1.0" encoding="ISO-8859-1"?>
<?xml-stylesheet type="text/xsl" href="OpenGeoSysGLI.xsl"?>

<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
    <name>cube_1x1x1_geometry</name>
    <points>
        <point id="0" x="0" y="0" z="0" name="left_front_bottom"/>
        <point id="1" x="0" y="0" z="1" name="left_front_top"/>
        <point id="2" x="0" y="1" z="1" name="left_back_top"/>
        <point id="3" x="0" y="1" z="0" name="left_back_bottom"/>

        <point id="4" x="1" y="0" z="0" name="right_front_bottom"/>
        <point id="5" x="1" y="0" z="1" name="right_front_top"/>
        <point id="6" x="1" y="1" z="1" name="right_back_top"/>
        <point id="7" x="1" y="1" z="0" name="right_back_bottom"/>

        <point id="8" x="0" y="1" z="0.1"/>
        <point id="9" x="0" y="0" z="0.1"/>
        <point id="10" x="0" y="1" z="0.2"/>
        <point id="11" x="0" y="0" z="0.2"/>
        <point id="12" x="0" y="1" z="0.300001"/>
        <point id="13" x="0" y="0" z="0.300001"/>

        <point id="14" x="0" y="1" z="0.3"/>
        <point id="15" x="0" y="0" z="0.3"/>
        <point id="16" x="0" y="1" z="0.400001"/>
        <point id="17" x="0" y="0" z="0.400001"/>

        <point id="18" x="0" y="1" z="0.399999"/>
        <point id="19" x="0" y="0" z="0.399999"/>
        <point id="20" x="0" y="1" z="0.500001"/>
        <point id="21" x="0" y="0" z="0.500001"/>

        <point id="22" x="0" y="1" z="0.499999"/>
        <point id="23" x="0" y="0" z="0.499999"/>
        <point id="24" x="0" y="1" z="0.600001"/>
        <point id="25" x="0" y="0" z="0.600001"/>

        <point id="26" x="0" y="1" z="0.6"/>
        <point id="27" x="0" y="0" z="0.6"/>
        <point id="28" x="0" y="1" z="0.700001"/>
        <point id="29" x="0" y="0" z="0.700001"/>

        <point id="30" x="0" y="1" z="0.7"/>
        <point id="31" x="0" y="0" z="0.7"/>
        <point id="32" x="0" y="1" z="0.800001"/>
        <point id="33" x="0" y="0" z="0.800001"/>

        <point id="34" x="0" y="1" z="0.8"/>
        <point id="35" x="0" y="0" z="0.8"/>
        <point id="36" x="0" y="1" z="0.900001"/>
        <point id="37" x="0" y="0" z="0.900001"/>

        <point id="38" x="0" y="1" z="0.9"/>
        <point id="39" x="0" y="0" z="0.9"/>
    </points>

    <polylines>
        <polyline id="0" name="vertical_polyline_left_front">
            <pnt>0</pnt>
            <pnt>1</pnt>
        </polyline>
        <polyline id="1" name="vertical_polyline_right_front">
            <pnt>4</pnt>
            <pnt>5</pnt>
        </polyline>
        <polyline id="2" name="vertical_polyline_left_back">
            <pnt>3</pnt>
            <pnt>2</pnt>
        </polyline>
        <polyline id="3" name="vertical_polyline_right_back">
            <pnt>7</pnt>
            <pnt>6</pnt>
        </polyline>
        <polyline id="4" name="horizontal_polyline_front_bottom">
            <pnt>0</pnt>
            <pnt>4</pnt>
        </polyline>
        <polyline id="5" name="horizontal_polyline_front_top">
            <pnt>1</pnt>
            <pnt>5</pnt>
        </polyline>
        <polyline id="6" name="horizontal_polyline_back_bottom">
            <pnt>3</pnt>
            <pnt>7</pnt>
        </polyline>
        <polyline id="7" name="horizontal_polyline_back_top">
            <pnt>2</pnt>
            <pnt>6</pnt>
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
            <pnt>5</pnt>
            <pnt>6</pnt>
        </polyline>
        <polyline id="11" name="horizontal_polyline_right_bottom">
            <pnt>4</pnt>
            <pnt>7</pnt>
        </polyline>
        <polyline id="12" name="left0">
            <pnt>0</pnt>
            <pnt>3</pnt>
            <pnt>8</pnt>
            <pnt>9</pnt>
            <pnt>0</pnt>
        </polyline>
        <polyline id="13" name="left1">
            <pnt>9</pnt>
            <pnt>8</pnt>
            <pnt>10</pnt>
            <pnt>11</pnt>
            <pnt>9</pnt>
        </polyline>
        <polyline id="13" name="left2">
            <pnt>11</pnt>
            <pnt>10</pnt>
            <pnt>12</pnt>
            <pnt>13</pnt>
            <pnt>11</pnt>
        </polyline>
        <polyline id="16" name="left3">
            <pnt>15</pnt>
            <pnt>14</pnt>
            <pnt>16</pnt>
            <pnt>17</pnt>
            <pnt>15</pnt>
        </polyline>
        <polyline id="17" name="left4">
            <pnt>19</pnt>
            <pnt>18</pnt>
            <pnt>20</pnt>
            <pnt>21</pnt>
            <pnt>19</pnt>
        </polyline>
        <polyline id="18" name="left5">
            <pnt>23</pnt>
            <pnt>22</pnt>
            <pnt>24</pnt>
            <pnt>25</pnt>
            <pnt>23</pnt>
        </polyline>
        <polyline id="19" name="left6">
            <pnt>27</pnt>
            <pnt>26</pnt>
            <pnt>28</pnt>
            <pnt>29</pnt>
            <pnt>27</pnt>
        </polyline>
        <polyline id="20" name="left7">
            <pnt>31</pnt>
            <pnt>30</pnt>
            <pnt>32</pnt>
            <pnt>33</pnt>
            <pnt>31</pnt>
        </polyline>
        <polyline id="21" name="left8">
            <pnt>35</pnt>
            <pnt>34</pnt>
            <pnt>36</pnt>
            <pnt>37</pnt>
            <pnt>35</pnt>
        </polyline>
        <polyline id="22" name="left9">
            <pnt>39</pnt>
            <pnt>38</pnt>
            <pnt>2</pnt>
            <pnt>1</pnt>
            <pnt>39</pnt>
        </polyline>
    </polylines>

    <surfaces>
        <surface id="0" name="left"><!-- x=0 -->
            <element p1="0" p2="1" p3="2"/>
            <element p1="0" p2="3" p3="2"/>
        </surface>
        <surface id="1" name="right"><!-- x=1 -->
            <element p1="4" p2="6" p3="5"/>
            <element p1="4" p2="6" p3="7"/>
        </surface>
        <surface id="2" name="top"><!-- z=1 -->
            <element p1="1" p2="2" p3="5"/>
            <element p1="5" p2="2" p3="6"/>
        </surface>
        <surface id="3" name="bottom"><!-- z=0 -->
            <element p1="0" p2="3" p3="4"/>
            <element p1="4" p2="3" p3="7"/>
        </surface>
        <surface id="4" name="front"><!-- y=0 -->
            <element p1="0" p2="1" p3="4"/>
            <element p1="4" p2="1" p3="5"/>
        </surface>
        <surface id="5" name="back"><!-- y=1 -->
            <element p1="2" p2="3" p3="6"/>
            <element p1="6" p2="3" p3="7"/>
        </surface>
    </surfaces>
</OpenGeoSysGLI>
