#--------------------------------------------------------------

# Global timestep output options
timeStepToStartOutputAt = 0
forceOutputAtFirstCall = False

# Global screenshot output options
imageFileNamePadding = 0
rescale_lookuptable = False

# Whether or not to request specific arrays from the adaptor.
requestSpecificArrays = False

# a root directory under which all Catalyst output goes
rootDirectory = ''

# makes a cinema D index table
make_cinema_table = False

#--------------------------------------------------------------
# Code generated from cpstate.py to create the CoProcessor.
# paraview version 5.8.0
#--------------------------------------------------------------

from paraview.simple import *
from paraview import coprocessing

# ----------------------- CoProcessor definition -----------------------


def CreateCoProcessor():
    def _CreatePipeline(coprocessor, datadescription):
        class Pipeline:
            # state file generated using paraview version 5.8.0

            # ----------------------------------------------------------------
            # setup the data processing pipelines
            # ----------------------------------------------------------------

            # trace generated using paraview version 5.8.0
            #
            # To ensure correct image size when batch processing, please search
            # for and uncomment the line `# renderView*.ViewSize = [*,*]`

            #### disable automatic camera reset on 'Show'
            paraview.simple._DisableFirstRenderCameraReset()

            # create a new 'Superquadric'
            ogs_output = coprocessor.CreateProducer(datadescription, 'input')

            # create a new 'Clip'
            clip1 = Clip(Input=ogs_output)
            clip1.ClipType = 'Plane'
            clip1.HyperTreeGridClipper = 'Plane'
            clip1.Scalars = ['POINTS', 'D1_left_bottom_N1_right']
            clip1.Value = 1.3376572416618493

            # init the 'Plane' selected for 'ClipType'
            clip1.ClipType.Origin = [0.5, 0.5, 0.0]

            # init the 'Plane' selected for 'HyperTreeGridClipper'
            clip1.HyperTreeGridClipper.Origin = [0.5, 0.5, 0.0]

            # create a new 'Contour'
            contour1 = Contour(Input=ogs_output)
            contour1.ContourBy = ['POINTS', 'D1_left_bottom_N1_right']
            contour1.Isosurfaces = [
                1.0, 1.075034942591522, 1.1500698851830442, 1.2251048277745662,
                1.3001397703660882, 1.3751747129576102, 1.4502096555491324,
                1.5252445981406544, 1.6002795407321764, 1.6753144833236986
            ]
            contour1.PointMergeMethod = 'Uniform Binning'

            # ----------------------------------------------------------------
            # finally, restore active source
            SetActiveSource(clip1)
            # ----------------------------------------------------------------

            # Now any catalyst writers
            xMLPUnstructuredGridWriter1 = servermanager.writers.XMLPUnstructuredGridWriter(
                Input=clip1)
            coprocessor.RegisterWriter(
                xMLPUnstructuredGridWriter1,
                filename='square_1e1_neumann_clip1_%t.pvtu',
                freq=1,
                paddingamount=0,
                DataMode='Appended',
                HeaderType='UInt64',
                EncodeAppendedData=False,
                CompressorType='None',
                CompressionLevel='6')

            xMLPPolyDataWriter1 = servermanager.writers.XMLPPolyDataWriter(
                Input=contour1)
            coprocessor.RegisterWriter(
                xMLPPolyDataWriter1,
                filename='square_1e1_neumann_contour1_%t.pvtp',
                freq=1,
                paddingamount=0,
                DataMode='Appended',
                HeaderType='UInt64',
                EncodeAppendedData=False,
                CompressorType='None',
                CompressionLevel='6')

        return Pipeline()

    class CoProcessor(coprocessing.CoProcessor):
        def CreatePipeline(self, datadescription):
            self.Pipeline = _CreatePipeline(self, datadescription)

    coprocessor = CoProcessor()
    # these are the frequencies at which the coprocessor updates.
    freqs = {}
    coprocessor.SetUpdateFrequencies(freqs)
    coprocessor.SetInitialOutputOptions(timeStepToStartOutputAt,
                                        forceOutputAtFirstCall)

    if rootDirectory:
        coprocessor.SetRootDirectory(rootDirectory)

    if make_cinema_table:
        coprocessor.EnableCinemaDTable()

    return coprocessor


#--------------------------------------------------------------
# Global variable that will hold the pipeline for each timestep
# Creating the CoProcessor object, doesn't actually create the ParaView pipeline.
# It will be automatically setup when coprocessor.UpdateProducers() is called the
# first time.
coprocessor = CreateCoProcessor()

#--------------------------------------------------------------
# Enable Live-Visualizaton with ParaView and the update frequency
coprocessor.EnableLiveVisualization(False, 1)

# ---------------------- Data Selection method ----------------------


def RequestDataDescription(datadescription):
    "Callback to populate the request for current timestep"
    global coprocessor

    # setup requests for all inputs based on the requirements of the
    # pipeline.
    coprocessor.LoadRequestedData(datadescription)


# ------------------------ Processing method ------------------------


def DoCoProcessing(datadescription):
    "Callback to do co-processing for current timestep"
    global coprocessor

    # Update the coprocessor by providing it the newly generated simulation data.
    # If the pipeline hasn't been setup yet, this will setup the pipeline.
    coprocessor.UpdateProducers(datadescription)

    # Write output data, if appropriate.
    coprocessor.WriteData(datadescription)

    # Write image capture (Last arg: rescale lookup table), if appropriate.
    coprocessor.WriteImages(datadescription,
                            rescale_lookuptable=rescale_lookuptable,
                            image_quality=0,
                            padding_amount=imageFileNamePadding)

    # Live Visualization, if enabled.
    coprocessor.DoLiveVisualization(datadescription, "localhost", 22222)
