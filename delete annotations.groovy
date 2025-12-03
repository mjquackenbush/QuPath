
def deleteAnnotation = getAnnotationObjects().findAll{it.getPathClass() == getPathClass("CD20")}
// Remove it
removeObjects(deleteAnnotation, true)