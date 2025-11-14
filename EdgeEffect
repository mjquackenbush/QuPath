import static qupath.lib.scripting.QP.*
import static qupath.lib.gui.scripting.QPEx.*

selectCells();
runPlugin('qupath.lib.algorithms.IntensityFeaturesPlugin', '{"pixelSizeMicrons":0.3244,"region":"CIRCLE","tileSizeMicrons":25.0,"channel1":true,"channel2":true,"channel3":true,"channel4":true,"doMean":true,"doStdDev":false,"doMinMax":false,"doMedian":false,"doHaralick":false,"haralickMin":NaN,"haralickMax":NaN,"haralickDistance":1,"haralickBins":32}')

def cells = getCellObjects()

//replace 'CD20-F480' with the channel of interest for background subtraction
def chName = 'CD20-F480'

//measurement will be store as '*channel name* difference'
cells.each {
    double difference=it.measurements[chName+': Cell: Mean'] - it.measurements['Circle: Diameter 25.0 µm: 0.32 µm per pixel: '+chName+': Mean']

    it.measurements[chName+'Background Difference'] = difference
}
