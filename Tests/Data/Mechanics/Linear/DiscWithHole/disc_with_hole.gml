<?xml version="1.0" encoding="ISO-8859-1"?>
<!DOCTYPE OGS-GML-DOM>
<OpenGeoSysGLI xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ogs="http://www.opengeosys.org">
 <name>disc_with_hole</name>
 <points>
  <point x="0.000000" y="2.000000" z="0.000000" id="0"/>
  <point x="0.000000" y="10.000000" z="0.000000" id="1"/>
  <point x="10.000000" y="10.000000" z="0.000000" id="2"/>
  <point x="10.000000" y="0.000000" z="0.000000" id="3"/>
  <point x="2.000000" y="0.000000" z="0.000000" id="4"/>
  <point x="1.997259" y="0.104672" z="0.000000" id="5"/>
  <point x="1.989044" y="0.209057" z="0.000000" id="6"/>
  <point x="1.975377" y="0.312869" z="0.000000" id="7"/>
  <point x="1.956295" y="0.415823" z="0.000000" id="8"/>
  <point x="1.931852" y="0.517638" z="0.000000" id="9"/>
  <point x="1.902113" y="0.618034" z="0.000000" id="10"/>
  <point x="1.867161" y="0.716736" z="0.000000" id="11"/>
  <point x="1.827091" y="0.813473" z="0.000000" id="12"/>
  <point x="1.782013" y="0.907981" z="0.000000" id="13"/>
  <point x="1.732051" y="1.000000" z="0.000000" id="14"/>
  <point x="1.677341" y="1.089278" z="0.000000" id="15"/>
  <point x="1.618034" y="1.175571" z="0.000000" id="16"/>
  <point x="1.554292" y="1.258641" z="0.000000" id="17"/>
  <point x="1.486290" y="1.338261" z="0.000000" id="18"/>
  <point x="1.414214" y="1.414214" z="0.000000" id="19"/>
  <point x="1.338261" y="1.486290" z="0.000000" id="20"/>
  <point x="1.258641" y="1.554292" z="0.000000" id="21"/>
  <point x="1.175571" y="1.618034" z="0.000000" id="22"/>
  <point x="1.089278" y="1.677341" z="0.000000" id="23"/>
  <point x="1.000000" y="1.732051" z="0.000000" id="24"/>
  <point x="0.907981" y="1.782013" z="0.000000" id="25"/>
  <point x="0.813473" y="1.827091" z="0.000000" id="26"/>
  <point x="0.716736" y="1.867161" z="0.000000" id="27"/>
  <point x="0.618034" y="1.902113" z="0.000000" id="28"/>
  <point x="0.517638" y="1.931852" z="0.000000" id="29"/>
  <point x="0.415823" y="1.956295" z="0.000000" id="30"/>
  <point x="0.312869" y="1.975377" z="0.000000" id="31"/>
  <point x="0.209057" y="1.989044" z="0.000000" id="32"/>
  <point x="0.104672" y="1.997259" z="0.000000" id="33"/>
 </points>
 <polylines>
  <polyline id="0" name="BOTTOM">
   <pnt>3</pnt>
   <pnt>4</pnt>
  </polyline>
  <polyline id="1" name="TOP">
   <pnt>1</pnt>
   <pnt>2</pnt>
  </polyline>
  <polyline id="2" name="RIGHT">
   <pnt>2</pnt>
   <pnt>3</pnt>
  </polyline>
  <polyline id="3" name="LEFT">
   <pnt>0</pnt>
   <pnt>1</pnt>
  </polyline>
  <polyline id="4" name="HOLE">
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
   <pnt>21</pnt>
   <pnt>22</pnt>
   <pnt>23</pnt>
   <pnt>24</pnt>
   <pnt>25</pnt>
   <pnt>26</pnt>
   <pnt>27</pnt>
   <pnt>28</pnt>
   <pnt>29</pnt>
   <pnt>30</pnt>
   <pnt>31</pnt>
   <pnt>32</pnt>
   <pnt>33</pnt>
   <pnt>0</pnt>
  </polyline>
  <polyline id="5" name="PLY_DOMAIN">
   <pnt>0</pnt>
   <pnt>1</pnt>
   <pnt>2</pnt>
   <pnt>3</pnt>
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
   <pnt>21</pnt>
   <pnt>22</pnt>
   <pnt>23</pnt>
   <pnt>24</pnt>
   <pnt>25</pnt>
   <pnt>26</pnt>
   <pnt>27</pnt>
   <pnt>28</pnt>
   <pnt>29</pnt>
   <pnt>30</pnt>
   <pnt>31</pnt>
   <pnt>32</pnt>
   <pnt>33</pnt>
   <pnt>0</pnt>
  </polyline>
 </polylines>
 <surfaces>
  <surface id="0" name="DOMAIN">
   <element p1="33" p2="1" p3="0"/>
   <element p1="33" p2="2" p3="1"/>
   <element p1="2" p2="4" p3="3"/>
   <element p1="2" p2="5" p3="4"/>
   <element p1="32" p2="2" p3="33"/>
   <element p1="2" p2="6" p3="5"/>
   <element p1="31" p2="2" p3="32"/>
   <element p1="2" p2="7" p3="6"/>
   <element p1="30" p2="2" p3="31"/>
   <element p1="2" p2="8" p3="7"/>
   <element p1="29" p2="2" p3="30"/>
   <element p1="2" p2="9" p3="8"/>
   <element p1="28" p2="2" p3="29"/>
   <element p1="2" p2="10" p3="9"/>
   <element p1="27" p2="2" p3="28"/>
   <element p1="2" p2="11" p3="10"/>
   <element p1="26" p2="2" p3="27"/>
   <element p1="2" p2="12" p3="11"/>
   <element p1="25" p2="2" p3="26"/>
   <element p1="2" p2="13" p3="12"/>
   <element p1="24" p2="2" p3="25"/>
   <element p1="2" p2="14" p3="13"/>
   <element p1="23" p2="2" p3="24"/>
   <element p1="2" p2="15" p3="14"/>
   <element p1="22" p2="2" p3="23"/>
   <element p1="2" p2="16" p3="15"/>
   <element p1="21" p2="2" p3="22"/>
   <element p1="2" p2="17" p3="16"/>
   <element p1="20" p2="2" p3="21"/>
   <element p1="2" p2="18" p3="17"/>
   <element p1="19" p2="2" p3="20"/>
   <element p1="2" p2="19" p3="18"/>
  </surface>
 </surfaces>
</OpenGeoSysGLI>
