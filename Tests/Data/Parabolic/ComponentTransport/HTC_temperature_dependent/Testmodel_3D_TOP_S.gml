<?xml version="1.0" encoding="utf-8"?>
<OpenGeoSysGLI>
  <name>Testmodel_3D_TOP</name>
  <points>
    <point id="0" x="0" y="0" z="0"/>
    <point id="1" x="1000" y="0" z="0"/>
    <point id="2" x="1000" y="400" z="0"/>
    <point id="3" x="0" y="400" z="0"/>
    <point id="4" x="0" y="400" z="-50"/>
    <point id="5" x="0" y="0" z="-50"/>
    <point id="6" x="390" y="0" z="-50"/>
    <point id="7" x="400" y="0" z="-50"/>
    <point id="8" x="410" y="0" z="-50"/>
    <point id="9" x="590" y="0" z="-50"/>
    <point id="10" x="600" y="0" z="-50"/>
    <point id="11" x="610" y="0" z="-50"/>
    <point id="12" x="1000" y="0" z="-50"/>
    <point id="13" x="1000" y="400" z="-50"/>
    <point id="14" x="1000" y="400" z="-100"/>
    <point id="15" x="0" y="400" z="-100"/>
    <point id="16" x="0" y="0" z="-100"/>
    <point id="17" x="390" y="0" z="-100"/>
    <point id="18" x="400" y="0" z="-100"/>
    <point id="19" x="410" y="0" z="-100"/>
    <point id="20" x="590" y="0" z="-100"/>
    <point id="21" x="600" y="0" z="-100"/>
    <point id="22" x="610" y="0" z="-100"/>
    <point id="23" x="1000" y="0" z="-100"/>
    <point id="24" x="1000" y="0" z="-150"/>
    <point id="25" x="1000" y="400" z="-150"/>
    <point id="26" x="0" y="400" z="-150"/>
    <point id="27" x="0" y="0" z="-150"/>
  </points>

  <polylines>
    <polyline id="1" name = "top_ton">
      <pnt>0</pnt>
      <pnt>1</pnt>
      <pnt>2</pnt>
      <pnt>3</pnt>
      <pnt>0</pnt>
    </polyline>
    <polyline id="2" name = "front_ton">
      <pnt>0</pnt>
      <pnt>1</pnt>
      <pnt>12</pnt>
      <pnt>11</pnt>
      <pnt>10</pnt>
      <pnt>9</pnt>
      <pnt>8</pnt>
      <pnt>7</pnt>
      <pnt>6</pnt>
      <pnt>5</pnt>
      <pnt>0</pnt>
    </polyline>
    <polyline id="3" name = "right_ton">
      <pnt>1</pnt>
      <pnt>2</pnt>
      <pnt>13</pnt>
      <pnt>12</pnt>
      <pnt>1</pnt>
    </polyline>
    <polyline id="4" name = "back_ton">
      <pnt>2</pnt>
      <pnt>3</pnt>
      <pnt>4</pnt>
      <pnt>13</pnt>
      <pnt>2</pnt>
    </polyline>
    <polyline id="5" name = "left_ton">
      <pnt>3</pnt>
      <pnt>4</pnt>
      <pnt>5</pnt>
      <pnt>0</pnt>
      <pnt>3</pnt>
    </polyline>
    <polyline id="6" name = "top_sand">
      <pnt>5</pnt>
      <pnt>6</pnt>
      <pnt>7</pnt>
      <pnt>8</pnt>
      <pnt>9</pnt>
      <pnt>10</pnt>
      <pnt>11</pnt>
      <pnt>12</pnt>
      <pnt>13</pnt>
      <pnt>4</pnt>
      <pnt>5</pnt>
    </polyline>
    <polyline id="7" name = "front_sand">
      <pnt>5</pnt>
      <pnt>6</pnt>
      <pnt>7</pnt>
      <pnt>8</pnt>
      <pnt>9</pnt>
      <pnt>10</pnt>
      <pnt>11</pnt>
      <pnt>12</pnt>
      <pnt>23</pnt>
      <pnt>22</pnt>
      <pnt>21</pnt>
      <pnt>20</pnt>
      <pnt>19</pnt>
      <pnt>18</pnt>
      <pnt>17</pnt>
      <pnt>16</pnt>
      <pnt>5</pnt>
    </polyline>
    <polyline id="8" name = "right_sand">
      <pnt>12</pnt>
      <pnt>13</pnt>
      <pnt>14</pnt>
      <pnt>23</pnt>
      <pnt>12</pnt>
    </polyline>
    <polyline id="9" name = "back_sand">
      <pnt>13</pnt>
      <pnt>4</pnt>
      <pnt>15</pnt>
      <pnt>14</pnt>
      <pnt>13</pnt>
    </polyline>
    <polyline id="10" name = "left_sand">
      <pnt>4</pnt>
      <pnt>5</pnt>
      <pnt>16</pnt>
      <pnt>15</pnt>
      <pnt>4</pnt>
    </polyline>
    <polyline id="11" name = "bottom_sand">
      <pnt>16</pnt>
      <pnt>17</pnt>
      <pnt>18</pnt>
      <pnt>19</pnt>
      <pnt>20</pnt>
      <pnt>21</pnt>
      <pnt>22</pnt>
      <pnt>23</pnt>
      <pnt>14</pnt>
      <pnt>15</pnt>
      <pnt>16</pnt>
    </polyline>
    <polyline id="12" name = "front_ton2">
      <pnt>16</pnt>
      <pnt>17</pnt>
      <pnt>18</pnt>
      <pnt>19</pnt>
      <pnt>20</pnt>
      <pnt>21</pnt>
      <pnt>22</pnt>
      <pnt>23</pnt>
      <pnt>24</pnt>
      <pnt>27</pnt>
      <pnt>16</pnt>
    </polyline>
    <polyline id="13" name = "right_ton2">
      <pnt>23</pnt>
      <pnt>14</pnt>
      <pnt>25</pnt>
      <pnt>24</pnt>
      <pnt>23</pnt>
    </polyline>
    <polyline id="14" name = "back_ton2">
      <pnt>14</pnt>
      <pnt>15</pnt>
      <pnt>26</pnt>
      <pnt>25</pnt>
      <pnt>14</pnt>
    </polyline>
    <polyline id="15" name = "left_ton2">
      <pnt>15</pnt>
      <pnt>26</pnt>
      <pnt>27</pnt>
      <pnt>16</pnt>
      <pnt>15</pnt>
    </polyline>
    <polyline id="16" name = "bottom_ton2">
      <pnt>27</pnt>
      <pnt>24</pnt>
      <pnt>25</pnt>
      <pnt>26</pnt>
      <pnt>27</pnt>
    </polyline>
    <polyline id="17" name = "injection">
      <pnt>7</pnt>
      <pnt>18</pnt>
    </polyline>
    <polyline id="18" name = "pump">
      <pnt>10</pnt>
      <pnt>21</pnt>
    </polyline>
  </polylines>

	<surfaces>
		<surface id="1" name ="top">
			<element p1="0" p2="1" p3="2"  />
			<element p1="2" p2="3" p3="0"  />
		</surface>
		<surface id="2" name="bottom">
			<element p1="27" p2="24" p3="25" />
			<element p1="25" p2="26" p3="27" />
		</surface>
	</surfaces>

</OpenGeoSysGLI>
