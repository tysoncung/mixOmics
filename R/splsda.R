#############################################################################################################
# Authors:
#   Florian Rohart, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
#   Kim-Anh Le Cao, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD

# created: 2011
# last modified: 05-10-2017
#
# Copyright (C) 2011
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#############################################################################################################


# ========================================================================================================
# splsda: perform a sPLS-DA
# this function is a particular setting of .mintBlock, the formatting of the input is checked in .mintWrapper
# ========================================================================================================

#' Sparse Partial Least Squares Discriminant Analysis (sPLS-DA)
#'
#' Function to perform sparse Partial Least Squares to classify samples
#' (supervised analysis) and select variables.
#'
#' \code{splsda} function fits an sPLS model with \eqn{1, \ldots ,}\code{ncomp}
#' components to the factor or class vector \code{Y}. The appropriate indicator
#' (dummy) matrix is created. Logratio transform and multilevel analysis are
#' performed sequentially as internal pre-processing step, through
#' \code{\link{logratio.transfo}} and \code{\link{withinVariation}}
#' respectively.
#'
#' Logratio can only be applied if the data do not contain any 0 value (for
#' count data, we thus advise the normalise raw data with a 1 offset).
#'
#' More details about the PLS modes in \code{?pls}.
## --------------------------------------------------------------------------------------- arguments
#' @inheritParams plsda
#' @inheritParams spls
# @param keepX numeric vector of length \code{ncomp}, the number of variables
# to keep in \eqn{X}-loadings. By default all variables are kept in the model.
## --------------------------------------------------------------------------------------- value
#' @return \code{splsda} returns an object of class \code{"splsda"}, a list that contains the following components:
#' \item{X}{the centered and standardized original predictor matrix.}
#' \item{Y}{the centered and standardized indicator response vector or matrix.}
#' \item{ind.mat}{the indicator matrix.} \item{ncomp}{the number of components
#' included in the model.} \item{keepX}{number of \eqn{X} variables kept in the
#' model on each component.} \item{variates}{list containing the variates.}
#' \item{loadings}{list containing the estimated loadings for the \code{X} and
#' \code{Y} variates.} \item{names}{list containing the names to be used for
#' individuals and variables.} \item{nzv}{list containing the zero- or
#' near-zero predictors information.} \item{tol}{the tolerance used in the
#' iterative algorithm, used for subsequent S3 methods} \item{iter}{Number of
#' iterations of the algorthm for each component} \item{max.iter}{the maximum
#' number of iterations, used for subsequent S3 methods} \item{scale}{boolean
#' indicating whether the data were scaled in MINT S3 methods}
#' \item{logratio}{whether logratio transformations were used for compositional
#' data} \item{explained_variance}{amount of variance explained per component
#' (note that contrary to PCA, this amount may not decrease as the aim of the
#' method is not to maximise the variance, but the covariance between X and the
#' dummy matrix Y).} \item{mat.c}{matrix of coefficients from the regression of
#' X / residual matrices X on the X-variates, to be used internally by
#' \code{predict}.} \item{defl.matrix}{residual matrices X for each dimension.}
#' @author Florian Rohart, Ignacio González, Kim-Anh Lê Cao.
#' @seealso \code{\link{spls}}, \code{\link{summary}}, \code{\link{plotIndiv}},
#' \code{\link{plotVar}}, \code{\link{cim}}, \code{\link{network}},
#' \code{\link{predict}}, \code{\link{perf}}, \code{\link{mint.block.splsda}},
#' \code{\link{block.splsda}} and http://www.mixOmics.org for more details.
#' @references On sPLS-DA: Lê Cao, K.-A., Boitard, S. and Besse, P. (2011).
#' Sparse PLS Discriminant Analysis: biologically relevant feature selection
#' and graphical displays for multiclass problems. \emph{BMC Bioinformatics}
#' \bold{12}:253. On log ratio transformations: Filzmoser, P., Hron, K.,
#' Reimann, C.: Principal component analysis for compositional data with
#' outliers. Environmetrics 20(6), 621-632 (2009) Lê Cao K.-A., Costello ME,
#' Lakis VA, Bartolo, F,Chua XY, Brazeilles R, Rondeau P. MixMC: Multivariate
#' insights into Microbial Communities. PLoS ONE, 11(8): e0160169 (2016). On
#' multilevel decomposition: Westerhuis, J.A., van Velzen, E.J., Hoefsloot,
#' H.C., Smilde, A.K.: Multivariate paired data analysis: multilevel plsda
#' versus oplsda. Metabolomics 6(1), 119-128 (2010) Liquet, B., Lê Cao K.-A.,
#' Hocini, H., Thiebaut, R.: A novel approach for biomarker selection and the
#' integration of repeated measures experiments from two assays. BMC
#' bioinformatics 13(1), 325 (2012)
#' @keywords regression multivariate
#' @examples
#'
#' ## First example
#' X <- breast.tumors$gene.exp
#' # Y will be transformed as a factor in the function,
#' # but we set it as a factor to set up the colors.
#' Y <- as.factor(breast.tumors$sample$treatment)
#'
#' res <- splsda(X, Y, ncomp = 2, keepX = c(25, 25))
#'
#'
#' # individual names appear
#' plotIndiv(res, ind.names = Y, legend = TRUE, ellipse =TRUE)
#'
#' \dontrun{
#' ## Second example: one-factor analysis with sPLS-DA, selecting a subset of variables
#' # as in the paper Liquet et al.
#' #--------------------------------------------------------------
#' X <- vac18$genes
#' Y <- vac18$stimulation
#' # sample indicates the repeated measurements
#' design <- data.frame(sample = vac18$sample)
#' Y = data.frame(stimul = vac18$stimulation)
#'
#' # multilevel sPLS-DA model
#' res.1level <- splsda(X, Y = Y, ncomp = 3, multilevel = design,
#' keepX = c(30, 137, 123))
#'
#' # set up colors for plotIndiv
#' col.stim <- c("darkblue", "purple", "green4","red3")
#' plotIndiv(res.1level, ind.names = Y, col.per.group = col.stim)
#'
#' ## Third example: two-factor analysis with sPLS-DA, selecting a subset of variables
#' # as in the paper Liquet et al.
#' #--------------------------------------------------------------
#'
#' X <- vac18.simulated$genes
#' design <- data.frame(sample = vac18.simulated$sample)
#' Y = data.frame( stimu = vac18.simulated$stimulation,
#' time = vac18.simulated$time)
#'
#' res.2level <- splsda(X, Y = Y, ncomp = 2, multilevel = design,
#' keepX = c(200, 200))
#'
#' plotIndiv(res.2level, group = Y$stimu, ind.names = vac18.simulated$time,
#' legend = TRUE, style = 'lattice')
#'
#'
#'
#' ## Fourth example: with more than two classes
#' # ------------------------------------------------
#'
#' X <- as.matrix(liver.toxicity$gene)
#' # Y will be transformed as a factor in the function,
#' # but we set it as a factor to set up the colors.
#' Y <- as.factor(liver.toxicity$treatment[, 4])
#'
#' splsda.liver <- splsda(X, Y, ncomp = 2, keepX = c(20, 20))
#'
#' # individual name is set to the treatment
#' plotIndiv(splsda.liver, ind.names = Y, ellipse = TRUE, legend = TRUE)
#'
#'
#' ## Fifth example: 16S data with multilevel decomposion and log ratio transformation
#' # ------------------------------------------------
#'
#' splsda.16S = splsda(
#' X = diverse.16S$data.TSS,  # TSS normalised data
#' Y =  diverse.16S$bodysite,
#' multilevel = diverse.16S$sample, # multilevel decomposition
#' ncomp = 2,
#' keepX =  c(10, 150),
#' logratio= 'CLR')  # CLR log ratio transformation
#'
#'
#' plotIndiv(splsda.16S, ind.names = FALSE, pch = 16, ellipse = TRUE, legend = TRUE)
#' #OTUs selected at the family level
#' diverse.16S$taxonomy[selectVar(splsda.16S, comp = 1)$name,'Family']
#'
#' }
#'

###########################################################
## generic function
#############################################################
#' @usage \S4method{splsda}{ANY}(X, Y, ncomp = 2,
#' mode = c("regression", "canonical", "invariant", "classic"),
#' keepX, scale = TRUE, tol = 1e-06, max.iter = 100, near.zero.var = FALSE,
#' logratio = c("none", "CLR"),  multilevel = NULL, all.outputs = TRUE)
## arguemnts for 'ANY' must be copied from internal to @usage ANY,for generic other arguments
## passed to methods can be added plus '...' so the methods can get '...' - if we only
## include X, RStudio won't suggest the rest automatically for autofill
#' @export
#' @rdname splsda
setGeneric("splsda", def = function(X, Y, formula=NULL, ncomp = 2, mode = c("regression", "canonical", "invariant", "classic"), keepX, scale = TRUE, tol = 1e-06, max.iter = 100, near.zero.var = FALSE, logratio = "none", multilevel = NULL, all.outputs = TRUE,...) standardGeneric("splsda"))

#############################################################
## internal function
#############################################################
.splsda = function(X,
                   Y,
                   ncomp = 2,
                   mode = c("regression", "canonical", "invariant", "classic"),
                   keepX,
                   scale = TRUE,
                   tol = 1e-06,
                   max.iter = 100,
                   near.zero.var = FALSE,
                   logratio = "none",   # one of "none", "CLR"
                   multilevel = NULL,
                   all.outputs = TRUE,
                   ret.call=FALSE){
  mc <- match.call.defaults() 
  mc <- .check_plsda(mc)
  mc$mode <- .matchArg(mode)
  logratio <- mc$logratio <- .matchArg(logratio)
  mc$DA <- TRUE
  mc$ret.call <- NULL ## not need by wrapper
  # # call to '.mintWrapper'
  mc[[1L]] <- quote(.mintWrapper)
  result <- eval(mc)
  # choose the desired output from 'result'
  out = list(
    call = match.call(),
    X = result$A[-result$indY][[1]],
    Y = if (is.null(multilevel))
    {
      Y
    } else {
      result$Y.factor
    },
    ind.mat = result$A[result$indY][[1]],
    ncomp = result$ncomp,
    mode = result$mode,
    keepA=result$keepA,
    keepX = result$keepX,
    keepY = result$keepY,
    variates = result$variates,
    loadings = result$loadings,
    loadings.star = result$loadings.star,
    names = result$names,
    tol = result$tol,
    iter = result$iter,
    max.iter = result$max.iter,
    nzv = result$nzv,
    scale = scale,
    logratio = logratio,
    explained_variance = result$explained_variance,#[-result$indY],
    input.X = result$input.X,
    mat.c = result$mat.c#,
  )

  class(out) = c("mixo_splsda","mixo_spls","DA")
  # output if multilevel analysis
  if (!is.null(multilevel))
  {
    out$multilevel = multilevel
    class(out) = c("mixo_mlsplsda",class(out))
  }

  return(invisible(out))
}

## ----------- Methods ----------- 
#### ANY ####
#' @export
#' @rdname splsda
setMethod('splsda', 'ANY', function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) {
  mget(names(formals()), sys.frame(sys.nframe())) ## just to evaluate
  
  ## legacy code
  if ( class(try(data)) %in% c("data.frame", "matrix") )
    .stop(message = "mixOmics arguments have changed.
              Please carefully read the documentation and try to use named 
              arguments such as splsda(X=mat, ...) as opposed to splsda(mat, ...).",
          .subclass = "defunct")
  
  if ( !is_null(data)) { ## data must be NULL
    .stop(message = "data should be a MultiAssayExperiment class, or NULL",
          .subclass = "inv_signature")
  }
  if (!is_null(formula)) { ## formula must be NULL
    .stop(message = "With numerical X and Y, formula should not be provided. 
              See ?splsda",
          .subclass = "inv_signature")
  }
  
  mc <- match.call()
  mc[-1L] <- lapply(mc[-1L], eval.parent)
  mc$data <- mc$formula <- NULL 
  mc[[1L]] <- quote(.splsda)
  result <- eval(mc)
  .call_return(result, mc$ret.call, mcr = match.call(), fun.name = 'splsda')
})

#### signature(data = 'MultiAssayExperiment', formula = "formula") ####
## expect X and Y to be NULL
#' @export
#' @rdname splsda
setMethod('splsda', signature(data = 'MultiAssayExperiment', formula = 'formula'), 
          function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) {
            mget(names(formals()), sys.frame(sys.nframe())) ## just to evaluate
            ## X and Y NULL or missing
            if ( !(is_null(X) && is_null(Y)) )
              .stop(message = "Where 'data' and 'formula' are provided 'X' and 'Y' should be NULL.", 
                    .subclass = "inv_signature")
            mc <- match.call()
            mc[-1L] <- lapply(mc[-1L], eval.parent)
            .sformula_checker(mc) ## check formula validity
            mc[c('Y', 'X')] <- as.character(formula[2:3])
            mc <- .get_xy(mc = mc)
            mc$data <- mc$formula <- NULL 
            mc[[1L]] <- quote(.splsda)
            result <- eval(mc)
            .call_return(result, mc$ret.call, mcr = match.call(), fun.name = 'splsda')
          })


#### signature(data != 'MultiAssayExperiment', formula = "formula") ####
## expect X and Y to be NULL
#' @export
#' @rdname splsda
setMethod('splsda', signature(formula = 'formula'), 
          function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) {
            mget(names(formals()), sys.frame(sys.nframe())) ## just to evaluate
            mc <- match.call()
            mc[-1L] <- lapply(mc[-1L], eval.parent)
            .sformula_checker(mc) ## check formula validity
            mc$X <- eval.parent(as.list(formula)[[3]], n = 2)
            mc$Y <- eval.parent(as.list(formula)[[2]], n = 2)
            # mc <- .get_xy(mc = mc)
            mc$data <- mc$formula <- NULL 
            mc[[1L]] <- quote(.splsda)
            result <- eval(mc)
            .call_return(result, mc$ret.call, mcr = match.call(), fun.name = 'splsda')
          })


#### signature(data = 'MultiAssayExperiment', formula != "formula") ####
## expect X and Y to be valid characters
#' @export
#' @rdname splsda
setMethod('splsda', signature(data = 'MultiAssayExperiment'), 
          function(data=NULL, X=NULL, Y=NULL, formula=NULL, ...) {
            mget(names(formals()), sys.frame(sys.nframe())) ## just to evaluate
            ## X and Y NULL or missing
            if (!is_null(formula)) { ## formula must be NULL
              .stop(message = "With numerical X and Y, formula should not be provided. See ?splsda", 
                    .subclass = "inv_signature")
            }
            mc <- match.call()
            mc[-1L] <- lapply(mc[-1L], eval.parent)
            mc <- .get_xy(mc = mc)
            mc$data <- mc$formula <- NULL 
            mc[[1L]] <- quote(.splsda)
            result <- eval(mc)
            .call_return(result, mc$ret.call, mcr = match.call(), fun.name = 'splsda')
          })
