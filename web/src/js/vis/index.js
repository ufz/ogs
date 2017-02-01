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
  reader.setUrl(options.fileURL).then((reader, dataset) => {
    console.log('Metadata loaded with the geometry', dataset);
    reader.loadData().then((reader, dataset) => {
      renderer.resetCamera();
      renderWindow.render();
      var indexPointData = 0;
      var indexCellData = 0;
      reader.getArrays().forEach((array, i) => {
        console.log('-', array.name, array.location, ':', array.enable);
        if (array.location == 'pointData') {
          console.log('  - Range: ',
            array.ds[0].pointData.arrays[indexPointData].data.ranges[0].min,
            '-',
            array.ds[0].pointData.arrays[indexPointData].data.ranges[0].max
          );
          indexPointData++;
        } else {
          console.log('  - Range: ',
            array.ds[0].cellData.arrays[indexCellData].data.ranges[0].min,
            '-',
            array.ds[0].cellData.arrays[indexCellData].data.ranges[0].max
          );
          indexCellData++;
        }
      });
    });
  });

  const lookupTable = vtkLookupTable.newInstance({
    hueRange: [0.0, 0.33],
    mappingRange: [-1, 1],
  });

  const mapper = vtkMapper.newInstance({
    interpolateScalarsBeforeMapping: true,
    colorMode: ColorMode.MAP_SCALARS,
    scalarMode: ScalarMode.USE_POINT_FIELD_DATA,
    useLookupTableScalarRange: true,
    lookupTable,
  });
  mapper.setInputConnection(reader.getOutputPort());
  mapper.setColorByArrayName('pressure');

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

// -----------------------------------------------------------
// globals for inspecting
// -----------------------------------------------------------
  global.mapper = mapper;
  global.actor = actor;
  global.renderer = renderer;
  global.renderWindow = renderWindow;

  global.reader = reader;
  global.lut = lookupTable;
}

export default { load };
