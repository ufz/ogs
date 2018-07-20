import vtkURLExtract              from 'vtk.js/Sources/Common/Core/URLExtract';
import vis from './vis';

const userParams = vtkURLExtract.extractURLParameters();

if (userParams.url || userParams.fileURL) {
  const exampleContainer = document.querySelector('.content');
  const rootBody = document.querySelector('body');
  const myContainer = exampleContainer || rootBody;
  vis.load(myContainer, userParams);
}
