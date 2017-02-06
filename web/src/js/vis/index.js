import vtkFullScreenRenderWindow from 'vtk.js/Sources/Rendering/Misc/FullScreenRenderWindow';
import vtkActor           from 'vtk.js/Sources/Rendering/Core/Actor';
import vtkLookupTable     from 'vtk.js/Sources/Common/Core/LookupTable';
import vtkMapper          from 'vtk.js/Sources/Rendering/Core/Mapper';
import { ColorMode, ScalarMode }  from 'vtk.js/Sources/Rendering/Core/Mapper/Constants';
import { AttributeTypes } from 'vtk.js/Sources/Common/DataModel/DataSetAttributes/Constants';
import { FieldDataTypes } from 'vtk.js/Sources/Common/DataModel/DataSet/Constants';
import vtkHttpDataSetReader       from 'vtk.js/Sources/IO/Core/HttpDataSetReader';
import controlPanel from './controller.html';

export function load(container, options) {
// ----------------------------------------------------------------------------
// Standard rendering code setup
// ----------------------------------------------------------------------------
  const fullScreenRenderer = vtkFullScreenRenderWindow.newInstance();
  const renderer = fullScreenRenderer.getRenderer();
  const renderWindow = fullScreenRenderer.getRenderWindow();

// -------------
// DataSetReader
// -------------
  const reader = vtkHttpDataSetReader.newInstance({enableArray: true, fetchGzip: true});
  const mapper = vtkMapper.newInstance({
    interpolateScalarsBeforeMapping: true,
    colorMode: ColorMode.MAP_SCALARS,
    useLookupTableScalarRange: true,
  });
  mapper.setInputConnection(reader.getOutputPort());
  const arrays = [];
  reader.setUrl(options.fileURL).then((reader, dataset) => {
    console.log('Metadata loaded with the geometry', dataset);
    reader.loadData().then((reader, dataset) => {
      renderer.resetCamera();
      renderWindow.render();
      var indexPointData = 0;
      var indexCellData = 0;
      const arraySelect = document.querySelector('#arraySelect');
      reader.getArrays().forEach((array, i) => {
        // UI
        var option = document.createElement('option');
        option.text = array.location + " – " + array.name;
        option.value = array.name;
        arraySelect.add(option);

        var range = [0.0,0.0];
        var pointData = false;
        if (array.location == 'pointData') {
          range[0] = array.ds[0].pointData.arrays[indexPointData].data.ranges[0].min;
          range[1] = array.ds[0].pointData.arrays[indexPointData].data.ranges[0].max;
          pointData = true;
          indexPointData++;
        } else {
          range[0] = array.ds[0].cellData.arrays[indexCellData].data.ranges[0].min;
          range[1] = array.ds[0].cellData.arrays[indexCellData].data.ranges[0].max;
          indexCellData++;
        }
        var lut = vtkLookupTable.newInstance({
          hueRange: [0.0, 0.33],
          mappingRange: range,
        });
        arrays.push({name: array.name, pointData: pointData, lut: lut})
      });
      if (arrays.length > 0)
        setColorArray(0);
    });
  });

  const actor = vtkActor.newInstance();
  var property = actor.getProperty();
  property.setEdgeVisibility(true);

  actor.setMapper(mapper);
  renderer.addActor(actor);

  renderer.resetCamera();
  renderWindow.render();


// -----------------------------------------------------------
// UI control handling
// -----------------------------------------------------------
  fullScreenRenderer.addController(controlPanel);
  const representationSelector = document.querySelector('.representations');
  const resolutionChange = document.querySelector('.resolution');
  representationSelector.addEventListener('change', (e) => {
    const newRepValue = Number(e.target.value);
    actor.getProperty().setRepresentation(newRepValue);
    renderWindow.render();
  });
  resolutionChange.addEventListener('input', (e) => {
    const value = Number(e.target.value);
    actor.getProperty().setOpacity(value / 100.0)
    renderWindow.render();
  });

  const arraySelect = document.querySelector('#arraySelect');

  arraySelect.addEventListener('change', (e) => {
    setColorArray(e.target.selectedIndex);
  });

// -----------------------------------------------------------
// globals for inspecting
// -----------------------------------------------------------
  global.mapper = mapper;
  global.actor = actor;
  global.renderer = renderer;
  global.renderWindow = renderWindow;

  global.reader = reader;

  // -----------------------------------------------------------
  // local functions
  // -----------------------------------------------------------
  function setColorArray(index) {
    var array = arrays[index];
    if (array.pointData) {
      mapper.setScalarModeToUsePointFieldData();
    }
    else {
      mapper.setScalarModeToUseCellFieldData();
    }
    mapper.setLookupTable(array.lut);
    mapper.setColorByArrayName(array.name);
    document.getElementById('range').textContent = array.lut.getMappingRange().toString();
    renderWindow.render();
  }
}

export default { load };
