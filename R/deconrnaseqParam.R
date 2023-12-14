#!/usr/bin/env R

# Author: Sean Maden

#' deconrnaseqParam-class
#'
#' Applies the DeconRNASeq::DeconRNASeq() deconvolution algorithm.
#' 
#' @details Main constructor for class \linkS4class{deconrnaseqParam}.
#' @rdname deconrnaseqParam-class
#' @seealso 
#' \linkS4class{referencebasedParam-class}
#' 
#' @references
#' Gong, T., et al. (2011) Optimal Deconvolution of Transcriptional Profiling 
#' Data Using Quadratic Programming with Application to Complex Clinical Blood 
#' Samples, PLoS One, 6, e27156.
#' 
#' @examples
#' exampleDataList <- getDeconvolutionExampleData()
#' newParam <- deconrnaseqParam(exampleDataList[["bulkExpression"]],exampleDataList[["referenceExpression"]],exampleDataList[["cellScaleFactors"]])
#' deconvolution(newParam)
#'
setClass("deconrnaseqParam", contains="referencebasedParam", 
         slots=c(use.scale = "logical", returnInfo = "logical"))

#' Make new object of class deconrnaseqParam
#'
#' Main constructor for class \linkS4class{deconrnaseqParam}.
#'
#' @param bulkExpression Bulk mixed signals matrix of samples, which can be 
#' matched to single-cell samples.
#' @param referenceExpression Signature matrix of cell type-specific signals. 
#' If not provided, can be computed from a
#' provided ExpressionSet containing single-cell data.
#' @param cellScaleFactors Cell size factor transformations of length equal to 
#' the K cell types to deconvolve.
#' @param use.scale Whether to scale or center data (default TRUE).
#' @param returnInfo Whether to return metadata and original method outputs 
#' with predicted proportions.
#' 
#' @references
#' Gong, T., et al. (2011) Optimal Deconvolution of Transcriptional Profiling 
#' Data Using Quadratic Programming with Application to Complex Clinical Blood 
#' Samples, PLoS One, 6, e27156.
#'
#' @export
deconrnaseqParam <- function(bulkExpression, referenceExpression, 
                             cellScaleFactors = NULL, use.scale = FALSE, 
                             returnInfo = FALSE) {
  if(is(use.scale, "NULL")){use.scale <- FALSE}
  new("deconrnaseqParam", bulkExpression = bulkExpression, 
      referenceExpression = referenceExpression, 
      cellScaleFactors = cellScaleFactors, 
      use.scale = use.scale, returnInfo = returnInfo)
}


#' Deconvolution method for class \linkS4class{deconrnaseqParam}
#' 
#' Main deconvolution method for the \linkS4class{deconrnaseqParam} to run the 
#' \code{DeconRNASeq::DeconRNASeq()} implementation of the DeconRNASeq 
#' algorithm.
#' 
#' @param object An object of class \linkS4class{deconrnaseqParam}.
#' 
#' @examples
#' exampleDataList <- getDeconvolutionExampleData()
#' newParam <- deconrnaseqParam(exampleDataList[["bulkExpression"]],exampleDataList[["referenceExpression"]],exampleDataList[["cellScaleFactors"]])
#' deconvolution(newParam)
#' 
#' @returns Either a vector of predicted proportions, or a list containing 
#' predictions, metadata, and original outputs.
#'
#' @references
#' Gong, T., et al. (2011) Optimal Deconvolution of Transcriptional Profiling 
#' Data Using Quadratic Programming with Application to Complex Clinical Blood 
#' Samples, PLoS One, 6, e27156.
#'
#' @export
setMethod("deconvolution", signature(object = "deconrnaseqParam"), function(object){
  require(DeconRNASeq)
  lparam <- callNextMethod()
  # instantiate and format objects
  bulkExpression<- lparam[["bulkExpression"]]
  referenceExpression <- lparam[["referenceExpression"]]
  cellScaleFactors <- lparam[["cellScaleFactors"]]
  proportions <- matrix(cellScaleFactors, nrow = 1)
  referenceExpression <- as.data.frame(referenceExpression)
  bulkExpression <- cbind(bulkExpression, bulkExpression) |> as.data.frame()
  use.scale <- object[["use.scale"]]
  result <- lapply(seq(ncol(bulkExpression)), function(bulkSampleIndex){
    bulkExpressionIter <- cbind(bulkExpression[,bulkSampleIndex,drop=F],
                                bulkExpression[,bulkSampleIndex,drop=F]) |>
      as.data.frame()
    DeconRNASeq::DeconRNASeq(datasets = bulkExpressionIter,
                             signatures = referenceExpression,
                             use.scale = use.scale,
                             proportions = NULL)
  }
  )
  names(result) <- colnames(bulkExpression)
  # get standardized return list
  predictions <- lapply(result, function(iter){iter$out.all[1,]})
  names(predictions) <- names(result)
  returnList <- parseDeconvolutionPredictionsResults(
    predictions, colnames(referenceExpression), colnames(bulkExpression))
  if(object[["returnInfo"]]){
    returnList <- list(
      predictions=predictions,
      result.info=result,
      metadata=parametersList[["metadata"]])}
  return(returnList)
})
