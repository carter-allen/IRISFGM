#' @rdname ProcessData
#' @export
setGeneric(name = "ProcessData", def = function(object, ...) standardGeneric("ProcessData"))

#' @rdname ReadFrom10X_folder
#' @export
setGeneric(name = "ReadFrom10X_folder", def = function(object, ...) standardGeneric("ReadFrom10X_folder"))

#' @rdname ReadFrom10X_h5
#' @export
setGeneric(name = "ReadFrom10X_h5", def = function(object, ...) standardGeneric("ReadFrom10X_h5"))

#' @param object 
#'
#' @param ... other parameter
#'
#' @export
#' @rdname AddMeta
setGeneric(name = "AddMeta", def = function(object, ...) standardGeneric("AddMeta"))

#' @rdname RunLTMG
#' @export
setGeneric(name = "RunLTMG", def = function(object, ...) standardGeneric("RunLTMG"))




#' @export
#' @rdname GetLTMGmatrix
setGeneric(name = "GetLTMGmatrix", def = function(object, ...) standardGeneric("GetLTMGmatrix"))



#' @export
#' @rdname CalBinarySingleSignal
setGeneric(name = "CalBinarySingleSignal", def = function(object) standardGeneric("CalBinarySingleSignal"))


#' @export
#' @rdname GetBinarySingleSignal
setGeneric(name = "GetBinarySingleSignal", def = function(object, ...) standardGeneric("GetBinarySingleSignal"))

#' @export
#' @rdname CalBinaryMultiSignal
setGeneric(name = "CalBinaryMultiSignal", def = function(object) standardGeneric("CalBinaryMultiSignal"))


#' @export
#' @rdname GetBinaryMultiSignal
setGeneric(name = "GetBinaryMultiSignal", def = function(object, ...) standardGeneric("GetBinaryMultiSignal"))


#' @export
#' @rdname RunDiscretization
setGeneric(name = "RunDiscretization", def = function(object, ...) standardGeneric("RunDiscretization"))


#' @export
#' @rdname RunBicluster
setGeneric(name = "RunBicluster", def = function(object, ...) standardGeneric("RunBicluster"))

# LTMG series

#' @export
#' @rdname RunDimensionReduction
setGeneric(name = "RunDimensionReduction", def = function(object, ...) standardGeneric("RunDimensionReduction"))


#' @export
#' @rdname RunClassification
setGeneric(name = "RunClassification", def = function(object, ...) standardGeneric("RunClassification"))


#' @export
#' @rdname FindMarker
setGeneric(name = "FindMarker", def = function(object, ...) standardGeneric("FindMarker"))

#' @export
#' @rdname FindGlobalMarkers
setGeneric(name = "FindGlobalMarkers", def = function(object, ...) standardGeneric("FindGlobalMarkers"))

#' @export
#' @rdname RunPathway
setGeneric(name = "RunPathway", def = function(object, ...) standardGeneric("RunPathway"))

#' @export
#' @rdname PlotHeatmap
setGeneric(name = "PlotHeatmap", def = function(object, ...) standardGeneric("PlotHeatmap"))

#' @export
#' @rdname PlotDimension
setGeneric(name = "PlotDimension", def = function(object, ...) standardGeneric("PlotDimension"))

#' @export
#' @rdname FindClassBasedOnMC
setGeneric(name = "FindClassBasedOnMC", def = function(object, ...) standardGeneric("FindClassBasedOnMC"))

#' @param object 
#'
#' @param ... other parameter
#'
#' @export
#' @rdname DotPlotPathway
setGeneric(name = "DotPlotPathway", def = function(object, ...) standardGeneric("DotPlotPathway"))

#' @export
#' @rdname PlotNetwork
setGeneric(name = "PlotNetwork", def = function(object, ...) standardGeneric("PlotNetwork"))

#' @export
#' @rdname PlotMeta
setGeneric(name = "PlotMeta", def = function(object, ...) standardGeneric("PlotMeta"))

#' @export
#' @rdname SubsetData
setGeneric(name = "SubsetData", def = function(object, ...) standardGeneric("SubsetData"))

#' @export
#' @rdname PlotModuleNetwork
setGeneric(name = "PlotModuleNetwork", def = function(object, ...) standardGeneric("PlotModuleNetwork"))
