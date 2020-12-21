
#' ProcessData
#'
#' @param ... other arguments passed to methods
#' @param object
#' @rdname ProcessData
#' @export
setGeneric(name = "ProcessData", def = function(object, ...) standardGeneric("ProcessData"))

#' ReadFrom10X_folder
#'
#' @param object
#' @rdname ReadFrom10X_folder
#' @export
setGeneric(name = "ReadFrom10X_folder", def = function(object) standardGeneric("ReadFrom10X_folder"))

#' ReadFrom10X_h5
#'
#' @param object
#' @rdname ReadFrom10X_h5
#' @export
setGeneric(name = "ReadFrom10X_h5", def = function(object) standardGeneric("ReadFrom10X_h5"))

#' AddMeta
#'
#' @param ... other arguments passed to methods
#' @param object
#' @export
#' @rdname AddMeta
setGeneric(name = "AddMeta", def = function(object, ...) standardGeneric("AddMeta"))

#' RunLTMG
#'
#' @param ... other arguments passed to methods
#' @param object
#' @rdname RunLTMG
#' @export
setGeneric(name = "RunLTMG", def = function(object, ...) standardGeneric("RunLTMG"))


#' GetLTMGmatrix
#'
#' @param ... other arguments passed to methods
#' @param object
#' @export
#' @rdname GetLTMGmatrix
setGeneric(name = "GetLTMGmatrix", def = function(object, ...) standardGeneric("GetLTMGmatrix"))


#' CalBinarySingleSignal
#'
#' @param object
#' @export
#' @rdname CalBinarySingleSignal
setGeneric(name = "CalBinarySingleSignal", def = function(object) standardGeneric("CalBinarySingleSignal"))

#' GetBinarySingleSignal
#'
#' @param ... other arguments passed to methods
#' @param object
#' @export
#' @rdname GetBinarySingleSignal
setGeneric(name = "GetBinarySingleSignal", def = function(object, ...) standardGeneric("GetBinarySingleSignal"))

#' CalBinaryMultiSignal
#'
#' @param object
#' @export
#' @rdname CalBinaryMultiSignal
setGeneric(name = "CalBinaryMultiSignal", def = function(object) standardGeneric("CalBinaryMultiSignal"))

#' GetBinaryMultiSignal
#'
#' @param ... other arguments passed to methods
#' @param object
#' @export
#' @rdname GetBinaryMultiSignal
setGeneric(name = "GetBinaryMultiSignal", def = function(object, ...) standardGeneric("GetBinaryMultiSignal"))

#' RunDiscretization
#'
#' @param ... other arguments passed to methods
#' @param object
#' @export
#' @rdname RunDiscretization
setGeneric(name = "RunDiscretization", def = function(object, ...) standardGeneric("RunDiscretization"))

#' RunBicluster
#'
#' @param ... other arguments passed to methods
#' @param object
#' @export
#' @rdname RunBicluster
setGeneric(name = "RunBicluster", def = function(object, ...) standardGeneric("RunBicluster"))

# LTMG series

#' RunDimensionReduction
#'
#' @param ... other arguments passed to methods
#' @param object
#' @export
#' @rdname RunDimensionReduction
setGeneric(name = "RunDimensionReduction", def = function(object, ...) standardGeneric("RunDimensionReduction"))

#' RunClassification
#'
#' @param ... other arguments passed to methods
#' @param object
#' @export
#' @rdname RunClassification
setGeneric(name = "RunClassification", def = function(object, ...) standardGeneric("RunClassification"))

#' FindMarker
#'
#' @param ... other arguments passed to methods
#' @param object
#' @export
#' @rdname FindMarker
setGeneric(name = "FindMarker", def = function(object, ...) standardGeneric("FindMarker"))

#' FindGlobalMarkers
#'
#' @param ... other arguments passed to methods
#' @param object
#' @export
#' @rdname FindGlobalMarkers
setGeneric(name = "FindGlobalMarkers", def = function(object, ...) standardGeneric("FindGlobalMarkers"))

#' RunPathway
#'
#' @param ... other arguments passed to methods
#' @param object
#' @export
#' @rdname RunPathway
setGeneric(name = "RunPathway", def = function(object, ...) standardGeneric("RunPathway"))

#' PlotHeatmap
#'
#' @param ... other arguments passed to methods
#' @param object
#' @export
#' @rdname PlotHeatmap
setGeneric(name = "PlotHeatmap", def = function(object, ...) standardGeneric("PlotHeatmap"))

#' PlotDimension
#'
#' @param ... other arguments passed to methods
#' @param object
#' @export
#' @rdname PlotDimension
setGeneric(name = "PlotDimension", def = function(object, ...) standardGeneric("PlotDimension"))

#' FindClassBasedOnMC
#'
#' @param ... other arguments passed to methods
#' @param object 
#' @export
#' @rdname FindClassBasedOnMC
setGeneric(name = "FindClassBasedOnMC", def = function(object, ...) standardGeneric("FindClassBasedOnMC"))

#' DotPlotPathway
#'
#' @param ... other arguments passed to methods
#' @param object 
#' @export
#' @rdname DotPlotPathway
setGeneric(name = "DotPlotPathway", def = function(object, ...) standardGeneric("DotPlotPathway"))

#' PlotNetwork
#'
#' @param ... other arguments passed to methods
#' @param object
#' @export
#' @rdname PlotNetwork
setGeneric(name = "PlotNetwork", def = function(object, ...) standardGeneric("PlotNetwork"))

#' PlotMeta
#'
#' @param ... other arguments passed to methods
#' @param object
#' @export
#' @rdname PlotMeta
setGeneric(name = "PlotMeta", def = function(object, ...) standardGeneric("PlotMeta"))

#' SubsetData
#'
#' @param ... other arguments passed to methods
#' @param object
#' @export
#' @rdname SubsetData
setGeneric(name = "SubsetData", def = function(object, ...) standardGeneric("SubsetData"))

#' PlotModuleNetwork
#'
#' @param ... other arguments passed to methods
#' @param object
#' @export
#' @rdname PlotModuleNetwork
setGeneric(name = "PlotModuleNetwork", def = function(object, ...) standardGeneric("PlotModuleNetwork"))
