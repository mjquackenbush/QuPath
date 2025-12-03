import static qupath.lib.scripting.QP.*
import static qupath.lib.gui.scripting.QPEx.*
selectAnnotations();
createAnnotationsFromDensityMap("b cell agg", [0: 9.0], "CD20", "SPLIT")
