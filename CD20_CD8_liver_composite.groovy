// Get all annotations classified as "Tissue"
def tissueAnnotations = getAnnotationObjects().findAll {
    it.getPathClass() == getPathClass("Tissue")
}

// Select only those Tissue annotations
clearSelectedObjects()
selectObjects(tissueAnnotations)

// Run the object classifier only on the selected annotations
runObjectClassifier("CD8a_CD20_composite_liver_object_classifier_11102025")
