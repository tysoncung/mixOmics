#' Plot of Individuals (Experimental Units)
#'
#' This function provides scatter plots for individuals (experimental units)
#' representation in (sparse)(I)PCA, (regularized)CCA, (sparse)PLS(DA) and
#' (sparse)(R)GCCA(DA).
#'
#' \code{plotIndiv} method makes scatter plot for individuals representation
#' depending on the subspace of projection. Each point corresponds to an
#' individual.
#'
#' If \code{ind.names=TRUE} and row names is \code{NULL}, then
#' \code{ind.names=1:n}, where \code{n} is the number of individuals. Also, if
#' \code{pch} is an input, then \code{ind.names} is set to FALSE as we do not
#' show both names and shapes.
#'
#' \code{plotIndiv} can have a two layers legend. This is especially convenient
#' when you have two grouping factors, such as a gender effect and a study
#' effect, and you want to highlight both simulatenously on the graphical
#' output. A first layer is coded by the \code{group} factor, the second by the
#' \code{pch} argument. When \code{pch} is missing, a single layer legend is
#' shown. If the \code{group} factor is missing, the \code{col} argument is
#' used to create the grouping factor \code{group}. When a second grouping
#' factor is needed and added via \code{pch}, \code{pch} needs to be a vector
#' of length the number of samples. In the case where \code{pch} is a vector or
#' length the number of groups, then we consider that the user wants a
#' different \code{pch} for each level of \code{group}. This leads to a single
#' layer legend and we merge \code{col} and \code{pch}. In the similar case
#' where \code{pch} is a single value, then this value is used to represent all
#' samples. See examples below for object of class plsda and splsda.
#'
#' In the specific case of a single `omics supervised model
#' (\code{\link{plsda}}, \code{\link{splsda}}), users can overlay prediction
#' results to sample plots in order to visualise the prediction areas of each
#' class, via the \code{background} input parameter. Note that this
#' functionality is only available for models with less than 2 components as
#' the surfaces obtained for higher order components cannot be projected onto a
#' 2D representation in a meaningful way. For more details, see
#' \code{\link{background.predict}}
#'
#' For customized plots (i.e. adding points, text), use the style = 'graphics'
#' (default is ggplot2).
#'
#' Note: the ellipse options were borrowed from the \pkg{ellipse}.
#'
#' @aliases plotIndiv plotIndiv.mixo_pls plotIndiv.mixo_spls plotIndiv.rcc
#' plotIndiv.pca plotIndiv.sipca plotIndiv.sgcca plotIndiv.rgcca
#' plotIndiv.mint.spls plotIndiv.mint.splsda
#' @param object object of class inherited from any \pkg{mixOmics}: \code{PLS,
#' sPLS, PLS-DA, SPLS-DA, rCC, PCA, sPCA, IPCA, sIPCA, rGCCA, sGCCA, sGCCDA}
#' @param comp integer vector of length two (or three to 3d). The components
#' that will be used on the horizontal and the vertical axis respectively to
#' project the individuals.
#' @param rep.space For objects inherited from \code{"rcc"}, \code{"pls"},
#' \code{"spls"}, character string, (partially) matching one of
#' \code{"X-variate"}, \code{"Y-variate"} ,or \code{"XY-variate"}, determining
#' the subspace to project the individuals. Defaults to \code{"X-variate"}
#' \code{"pca"} object and for \code{"plsda"} objects. For objects of class
#' \code{"pls"} and \code{"rcc"},defaults, the tree subspaces represent the
#' individuals. For objects of class \code{"rgcca"} and \code{"sgcca"},
#' numerical value indicating the block data set form which to represent the
#' individuals.
#' @param blocks integer value of name of a block to be plotted using the GCCA
#' module. See examples.
#' @param study Indicates which study-specific outputs to plot. A character
#' vector containing some levels of \code{object$study}, "all.partial" to plot
#' all studies or "global" is expected. Default to "global".
#' @param ind.names either a character vector of names for the individuals to
#' be plotted, or \code{FALSE} for no names. If \code{TRUE}, the row names of
#' the first (or second) data matrix is used as names (see Details).
#' @param group factor indicating the group membership for each sample, useful
#' for ellipse plots. Coded as default for the supervised methods \code{PLS-DA,
#' SPLS-DA,sGCCDA}, but needs to be input for the unsupervised methods
#' \code{PCA, sPCA, IPCA, sIPCA, PLS, sPLS, rCC, rGCCA, sGCCA}
#' @param col.per.group character (or symbol) color to be used when 'group' is
#' defined. Vector of the same length than the number of groups.
#' @param style argument to be set to either \code{'graphics'},
#' \code{'lattice'}, \code{'ggplot2'} or \code{'3d'} for a style of plotting.
#' Default set to 'ggplot2'. See details. \code{3d} is not available for MINT
#' objects.
#' @param ellipse boolean indicating if ellipse plots should be plotted. In the
#' non supervised objects \code{PCA, sPCA, IPCA, sIPCA, PLS, sPLS, rCC, rGCCA,
#' sGCCA} ellipse plot is only be plotted if the argument \code{group} is
#' provided. In the \code{PLS-DA, SPLS-DA,sGCCDA} supervised object, by default
#' the ellipse will be plotted accoding to the outcome \code{Y}.
#' @param ellipse.level Numerical value indicating the confidence level of
#' ellipse being plotted when \code{ellipse =TRUE} (i.e. the size of the
#' ellipse). The default is set to 0.95, for a 95\% region.
#' @param centroid boolean indicating whether centroid points should be
#' plotted. In the non supervised objects \code{PCA, sPCA, IPCA, sIPCA, PLS,
#' sPLS, rCC, rGCCA, sGCCA} the centroid will only be plotted if the argument
#' \code{group} is provided. The centroid will be calculated based on the group
#' categories. In the supervised objects \code{PLS-DA, SPLS-DA,sGCCDA} the
#' centroid will be calculated according to the outcome \code{Y}.
#' @param star boolean indicating whether a star plot should be plotted, with
#' arrows starting from the centroid (see argument \code{centroid}, and ending
#' for each sample belonging to each group or outcome. In the non supervised
#' objects \code{PCA, sPCA, IPCA, sIPCA, PLS, sPLS, rCC, rGCCA, sGCCA} star
#' plot is only be plotted if the argument \code{group} is provided. In the
#' supervised objects \code{PLS-DA, SPLS-DA,sGCCDA} the star plot is plotted
#' according to the outcome \code{Y}.
#' @param title set of characters indicating the title plot.
#' @param subtitle subtitle for each plot, only used when several \code{block}
#' or \code{study} are plotted.
#' @param legend boolean. Whether the legend should be added. Default is FALSE.
#' @param X.label x axis titles.
#' @param Y.label y axis titles.
#' @param Z.label z axis titles (when style = '3d').
#' @param abline should the vertical and horizontal line through the center be
#' plotted? Default set to \code{FALSE}
#' @param xlim,ylim numeric list of vectors of length 2 and length
#' =length(blocks), giving the x and y coordinates ranges.
#' @param col character (or symbol) color to be used, possibly vector.
#' @param cex numeric character (or symbol) expansion, possibly vector.
#' @param pch plot character. A character string or a vector of single
#' characters or integers. See \code{\link{points}} for all alternatives.
#' @param pch.levels Only used when \code{pch} is different from \code{col} or
#' \code{col.per.group}, ie when \code{pch} creates a second factor. Only used
#' for the legend.
#' @param alpha Semi-transparent colors (0 < \code{'alpha'} < 1)
#' @param axes.box for style '3d', argument to be set to either \code{'axes'},
#' \code{'box'}, \code{'bbox'} or \code{'all'}, defining the shape of the box.
#' @param layout layout parameter passed to mfrow. Only used when \code{study}
#' is not "global"
#' @param size.title size of the title
#' @param size.subtitle size of the subtitle
#' @param size.xlabel size of xlabel
#' @param size.ylabel size of ylabel
#' @param size.axis size of the axis
#' @param size.legend size of the legend
#' @param size.legend.title size of the legend title
#' @param legend.title title of the legend
#' @param legend.title.pch title of the second legend created by \code{pch}, if
#' any.
#' @param legend.position position of the legend, one of "bottom", "left",
#' "top" and "right".
#' @param point.lwd \code{lwd} of the points, used when \code{ind.names =
#' FALSE}
#' @param background color the background by the predicted class, see
#' \code{\link{background.predict}}
#' @param ... other arguments passed to function methods.
#' @return none
#' @author Ignacio González, Benoit Gautier, Francois Bartolo, Florian Rohart, Al J Abadi
#' @seealso \code{\link{text}}, \code{\link{background.predict}},
#' \code{\link{points}} and http://mixOmics.org/graphics for more details.
#' @keywords multivariate hplot dplot
#' @example examples/plotIndiv-example.R
## ----------- Generic ----------- 

#' @importFrom ellipse ellipse
#' @export plotIndiv
plotIndiv <- function(object, ...) UseMethod("plotIndiv")
##TODO check tests and examples
## -----------  (s)PLS(DA) ----------- 
#----------------------------------------------------------------------------------------------------------#
#-- Includes plotIndiv for PLS, sPLS, PLS-DA, SPLS-DA,  --#
#----------------------------------------------------------------------------------------------------------#
#' @title PLS sample plot methods
#' @rdname plotIndiv
#' @export
#' @method plotIndiv mixo_pls
plotIndiv.mixo_pls <-  function(object,
                                comp  = NULL,
                                rep.space  = NULL,
                                ind.names  = TRUE,
                                group, # factor indicating the group membership for each sample, useful for ellipse plots. Coded as default for the -da methods, but needs to be input for the unsupervised methods (PCA, IPCA...)
                                col.per.group,
                                style = "ggplot2", # can choose between graphics, 3d, lattice or ggplot2
                                ellipse  = FALSE, #ellipse
                                ellipse.level  = 0.95,
                                centroid = FALSE,  # centroid
                                star = FALSE, # star
                                title = NULL, #title
                                subtitle,
                                legend = FALSE,
                                X.label  = NULL,
                                Y.label  = NULL,
                                Z.label  = NULL,
                                abline  = FALSE, #abline
                                xlim  = NULL,
                                ylim  = NULL,
                                col,
                                cex,
                                pch,
                                pch.levels,
                                alpha = 0.2, # used in shade3d
                                axes.box  = "box",
                                layout = NULL,
                                size.title = rel(2),
                                size.subtitle = rel(1.5),
                                size.xlabel = rel(1),
                                size.ylabel = rel(1),
                                size.axis = rel(0.8),
                                size.legend = rel(1), #size.legend
                                size.legend.title = rel(1.1), #size.legend.title
                                legend.title = "Legend",
                                legend.title.pch = "Legend",
                                legend.position = "right",
                                point.lwd = 1,
                                background = NULL
)
{
    plot_parameters = list(size.title = size.title, size.subtitle = size.subtitle, size.xlabel = size.xlabel, size.ylabel = size.ylabel,
                           size.axis = size.axis, size.legend = size.legend, size.legend.title = size.legend.title, legend.title = legend.title,
                           legend.title.pch = legend.title.pch, legend.position = legend.position, point.lwd = point.lwd)
    
    if (is(object, c("mint.block.pls", "mint.block.spls", "mint.block.plsda", "mint.block.splsda")))
        stop("No plotIndiv for the following functions at this stage: mint.block.pls, mint.block.spls, mint.block.plsda, mint.block.splsda.")
    
    #-- choose rep.space
    if (is.null(rep.space) && is(object, "DA"))#"splsda", "plsda", "mlsplsda")))
    {
        rep.space = "X-variate"
    } else if (is.null(rep.space)) {
        rep.space = "multi"
    }
    rep.space  = match.arg(rep.space, c("XY-variate", "X-variate", "Y-variate", "multi"))
    #c("XY-variate", "X-variate", "Y-variate", "multi")[pmatch(rep.space, c("XY-variate", "X-variate", "Y-variate", "multi"))]
    
    if (rep.space  == "multi")
    {
        blocks = c("X", "Y")
        object$variates  = object$variates[names(object$variates) %in% blocks]
    }
    
    if (rep.space  == "X-variate")
    {
        object$variates  = object$variates["X"]
        blocks  = "X"
    }
    
    if (rep.space  == "Y-variate")
    {
        object$variates  = object$variates["Y"]
        blocks  = "Y"
    }
    
    if (rep.space  == "XY-variate")
    {
        object$variates$XYvariates  = (object$variates$X + object$variates$Y)/2
        object$variates  = object$variates["XYvariates"]
        blocks  = "XY combined"
    }
    
    if (length(blocks)!= length(unique(blocks)))
        stop("Duplicate in 'blocks' not allowed")
    
    if (!isNULL(subtitle))
    {
        if (length(subtitle)!= length(blocks) | length(subtitle)!= length(unique(subtitle)))
            stop("'subtitle' indicates the subtitle of the plot for each 'blocks'; it needs to be the same length as 'blocks' and duplicate are not allowed.")
    }
    
    if(!is.null(background) &&  !is(background, "background.predict"))
        stop("'background' must have been obtained with the 'background.predict' function")
    
    #-- check inputs
    check  = .plotIndivCheckInput(object = object, comp  = comp , blocks  = blocks, ind.names  = ind.names,
                                  style = style, ellipse  = ellipse, ellipse.level  = ellipse.level, centroid = centroid,
                                  star = star, legend = legend, X.label  = X.label, Y.label  = Y.label, Z.label  = Z.label, abline  = abline,
                                  xlim  = xlim, ylim  = ylim, alpha = alpha, axes.box  = axes.box, plot_parameters = plot_parameters)
    #-- retrieve outputs from the checks
    axes.box = check$axes.box
    comp = check$comp
    xlim = check$xlim
    ylim = check$ylim
    ind.names = check$ind.names
    display.names = check$display.names
    
    #-- get the variates
    variate = .getVariatesAndLabels(object, comp, blocks = blocks, rep.space = rep.space, style = style, X.label = X.label,
                                    Y.label = Y.label, Z.label = Z.label)
    #-- retrieve outputs
    x = variate$x
    y = variate$y
    z = variate$z
    X.label = variate$X.label
    Y.label = variate$Y.label
    Z.label = variate$Z.label
    
    n = nrow(object$X)
    
    # create data frame df that contains (almost) all the ploting information
    out = .inputShapePlotIndiv(object = object, n = n, blocks  = blocks, x = x, y = y, z = z, ind.names  = ind.names, group = group,
                               col.per.group = col.per.group, style = style, study = "global", ellipse  = ellipse, ellipse.level  = ellipse.level,
                               centroid = centroid, star = star, title = title, xlim  = xlim, ylim  = ylim,
                               col = col, cex = cex, pch = pch, pch.levels = pch.levels, display.names = display.names, plot_parameters = plot_parameters)
    #-- retrieve outputs
    df = out$df
    df.ellipse = out$df.ellipse
    col.per.group = out$col.per.group
    title = out$title
    display.names = out$display.names
    xlim = out$xlim
    ylim = out$ylim
    #missing.col = out$missing.col
    ellipse = out$ellipse
    centroid = out$centroid
    star = out$star
    plot_parameters = out$plot_parameters
    
    # change the levels of df$Block to "subtitle"
    if (!isNULL(subtitle) & nlevels(df$Block)>1)#& !is.null(title)) # commented so that subtitle can be change without changing the title
    {
        df$Block = factor(df$Block, labels = subtitle)
        if (ellipse)
            df.ellipse$Block = factor(df.ellipse$Block, labels = subtitle)
    }
    
    # match background color to col.per.group, the color of the groups
    if(!is.null(background))
    {
        ind.match = match(names(background), levels(df$group))
        names(background) = adjustcolor(col.per.group[ind.match],alpha.f=0.1)
    }
    
    #save(list = ls(), file = "temp.Rdata")
    
    #call plot module (ggplot2, lattice, graphics, 3d)
    res = .graphicModule(df = df, centroid = centroid, col.per.group = col.per.group, title = title,
                         X.label = X.label, Y.label = Y.label, Z.label = Z.label, xlim = xlim, ylim = ylim, class.object = class(object),
                         display.names = display.names, legend = legend, abline = abline, star = star,
                         ellipse = ellipse, df.ellipse = df.ellipse, style = style, layout = layout, #missing.col = missing.col,
                         axes.box = axes.box, plot_parameters = plot_parameters, alpha = alpha, background = background)
    
    
    return(invisible(list(df = df, df.ellipse = df.ellipse, graph = res)))
    
}

#' @export
#' @method plotIndiv mixo_spls
#' @rdname plotIndiv
plotIndiv.mixo_spls <- plotIndiv.mixo_pls

#' @export
#' @rdname plotIndiv
plotIndiv.rcc <- plotIndiv.mixo_pls

## ----------- PCA family ----------- 
#----------------------------------------------------------------------------------------------------------#
#-- Includes plotIndiv for PCA, sPCA, IPCA, sIPCA --#
#----------------------------------------------------------------------------------------------------------#

#' @export
#' @method plotIndiv pca
#' @rdname plotIndiv
plotIndiv.pca <-
    function(object,
             comp = NULL,
             ind.names = TRUE,
             group,
             # factor indicating the group membership for each sample, useful for ellipse plots. Coded as default for the -da methods, but needs to be input for the unsupervised methods (PCA, IPCA...)
             col.per.group,
             style = "ggplot2",
             # can choose between graphics, 3d, lattice or ggplot2
             ellipse = FALSE,
             ellipse.level = 0.95,
             centroid = FALSE,
             star = FALSE,
             title = NULL,
             legend = FALSE,
             X.label = NULL,
             Y.label = NULL,
             Z.label = NULL,
             abline = FALSE,
             xlim = NULL,
             ylim = NULL,
             col,
             cex,
             pch,
             pch.levels,
             alpha = 0.2,
             axes.box = "box",
             layout = NULL,
             size.title = rel(2),
             size.subtitle = rel(1.5),
             size.xlabel = rel(1),
             size.ylabel = rel(1),
             size.axis = rel(0.8),
             size.legend = rel(1),
             size.legend.title = rel(1.1),
             legend.title = "Legend",
             legend.title.pch = "Legend",
             legend.position = "right",
             point.lwd = 1
    )
    {
        plot_parameters = list(size.title = size.title, size.subtitle = size.subtitle, size.xlabel = size.xlabel, size.ylabel = size.ylabel,
                               size.axis = size.axis, size.legend = size.legend, size.legend.title = size.legend.title, legend.title = legend.title,
                               legend.title.pch = legend.title.pch, legend.position = legend.position, point.lwd = point.lwd)
        
        blocks = "X"
        rep.space = "X-variate"
        
        check = .plotIndivCheckInput(object = object, comp = comp, blocks = blocks, ind.names = ind.names,
                                     style = style, ellipse = ellipse, ellipse.level = ellipse.level, centroid = centroid,
                                     star = star, legend = legend, X.label = X.label, Y.label = Y.label, Z.label = Z.label, abline = abline,
                                     xlim = xlim, ylim = ylim, alpha = alpha, axes.box = axes.box, plot_parameters = plot_parameters)
        
        # retrieve outputs from the checks
        axes.box = check$axes.box
        comp = check$comp
        xlim = check$xlim
        ylim = check$ylim
        ind.names = check$ind.names
        display.names = check$display.names
        
        
        #-- Get variates
        x = y = z = list()
        x[[1]] = object$x[, comp[1]]
        y[[1]] = object$x[, comp[2]]
        if(style == "3d") z[[1]] = object$x[, comp[3]]
        
        
        #-- Variance explained on X, Y and Z labels
        
        if (style ==  "3d")
        {
            inf = object$explained_variance[c(comp[1], comp[2], comp[3])]
            inf = round(inf, 2)
        } else {
            inf = object$explained_variance[c(comp[1], comp[2])]
            inf = round(inf, 2)}
        
        
        if (is.null(X.label))
        {
            X.label = paste("PC", comp[1], sep = '')
            percentage = paste0(inf[1]*100, "% expl. var")
            X.label = paste(X.label, percentage, sep = ": ")
        }
        if (is.null(Y.label))
        {
            Y.label = paste("PC", comp[2], sep = '')
            percentage = paste0(inf[2]*100, "% expl. var")
            Y.label = paste(Y.label, percentage, sep = ": ")
        }
        if (is.null(Z.label)&&style == "3d")
        {
            Z.label = paste("PC", comp[3], sep = '')
            percentage = paste0(inf[3]*100, "% expl. var")
            Z.label = paste(Z.label, percentage, sep = ": ")
        }
        
        
        n = nrow(object$X)
        
        # create data frame df that contains (almost) all the ploting information
        out = .inputShapePlotIndiv(object = object, n = n, blocks = blocks, x = x, y = y, z = z, ind.names = ind.names, group = group,
                                   col.per.group = col.per.group, style = style, study = "global", ellipse = ellipse, ellipse.level = ellipse.level,
                                   centroid = centroid, star = star, title = title, xlim = xlim, ylim = ylim,
                                   col = col, cex = cex, pch = pch, pch.levels = pch.levels, display.names = display.names, plot_parameters = plot_parameters)
        #-- retrieve outputs
        df = out$df
        df.ellipse = out$df.ellipse
        col.per.group = out$col.per.group
        title = out$title
        display.names = out$display.names
        xlim = out$xlim
        ylim = out$ylim
        #missing.col = out$missing.col
        ellipse = out$ellipse
        centroid = out$centroid
        star = out$star
        plot_parameters = out$plot_parameters
        
        #call plot module (ggplot2, lattice, graphics, 3d)
        res = .graphicModule(df = df, centroid = centroid, col.per.group = col.per.group, title = title,
                             X.label = X.label, Y.label = Y.label, Z.label = Z.label, xlim = xlim, ylim = ylim, class.object = class(object),
                             display.names = display.names, legend = legend, abline = abline,
                             star = star, ellipse = ellipse, df.ellipse = df.ellipse, style = style, layout = layout, #missing.col = missing.col,
                             axes.box = axes.box, plot_parameters = plot_parameters, alpha = alpha)
        
        return(invisible(list(df = df, df.ellipse = df.ellipse, graph = res)))
    }

## ----------- MINT ----------- 
#----------------------------------------------------------------------------------------------------------#
#-- Includes plotIndiv for the MINT module --#
#----------------------------------------------------------------------------------------------------------#

#' @export
#' @method plotIndiv mint.pls
#' @rdname plotIndiv
plotIndiv.mint.pls <- function(object,
                               comp = NULL,
                               study = "global",
                               rep.space = NULL,
                               group,
                               # factor indicating the group membership for each sample, useful for ellipse plots. Coded as default for the -da methods, but needs to be input for the unsupervised methods (PCA, IPCA...)
                               col.per.group,
                               style = "ggplot2",
                               # can choose between graphics, lattice or ggplot2
                               ellipse = FALSE,
                               ellipse.level = 0.95,
                               centroid = FALSE,
                               star = FALSE,
                               title = NULL,
                               subtitle,
                               legend = FALSE,
                               X.label = NULL,
                               Y.label = NULL,
                               abline = FALSE,
                               xlim = NULL,
                               ylim = NULL,
                               col,
                               cex,
                               pch,
                               layout = NULL,
                               size.title = rel(2),
                               size.subtitle = rel(1.5),
                               size.xlabel = rel(1),
                               size.ylabel = rel(1),
                               size.axis = rel(0.8),
                               size.legend = rel(1),
                               size.legend.title = rel(1.1),
                               legend.title = "Legend",
                               legend.position = "right",
                               point.lwd = 1
                               
)
    {
        plot_parameters = list(size.title = size.title, size.subtitle = size.subtitle, size.xlabel = size.xlabel, size.ylabel = size.ylabel, size.axis = size.axis,
                               size.legend = size.legend, size.legend.title = size.legend.title, legend.title = legend.title,
                               legend.position = legend.position, point.lwd = point.lwd)
        
        
        if (any(class(object)%in%c("mint.block.pls", "mint.block.spls", "mint.block.plsda", "mint.block.splsda")))
            stop("No plotIndiv for the following functions at this stage: mint.block.pls, mint.block.spls, mint.block.plsda, mint.block.splsda.")
        
        
        #-- rep.space
        if (is.null(rep.space))#"splsda", "plsda", "mlsplsda")))
            rep.space = "X-variate"
        
        rep.space  =  match.arg(rep.space, c("XY-variate", "X-variate", "Y-variate", "multi"))
        
        ind.names = FALSE
        # -------------------------------------------------------------------------------------- #
        #           need study
        # -------------------------------------------------------------------------------------- #
        
        # check study
        #study needs to be either: from levels(object$study), numbers from 1:nlevels(study) or "global"
        if (any(!study%in%c(levels(object$study), "global" , "all.partial")))
            stop("'study' must be one of 'object$study', 'global' or 'all.partial', see help file.")
        
        if (length(study)!=length(unique(study)))
            stop("Duplicate in 'study' not allowed")
        
        if (any(study != "global"))
        {
            if (ellipse == TRUE)
                stop("'ellipse' must be FALSE when study is different from 'global'")
            
            if (star == TRUE)
                stop("'star' must be FALSE when study is different from 'global'")
        }
        
        
        #LOOP ON STUDY, to get a plot with every single one, could be a mixed of numbers and "global", only if there is both "global" and something else.
        
        object.init = object
        study.init = unique(study)
        
        # replace "all.partial" by all levels of object$study
        ind.all.partial = which(study.init == "all.partial")
        if (length(ind.all.partial) > 0)
        {
            if (ind.all.partial > 1 & ind.all.partial < length(study.init))
            {
                # there are things before and after "all.partial"
                study.init = c(study.init[1:(ind.all.partial-1)], levels(object$study), study.init[(ind.all.partial+1) : length(study.init)])
            } else if (ind.all.partial == 1 & ind.all.partial < length(study.init)) {
                # there are only things after "all.partial"
                study.init = c(levels(object$study), study.init[(ind.all.partial+1) : length(study.init)])
            } else if (ind.all.partial > 1 & ind.all.partial == length(study.init)) {
                # there are things only before "all.partial"
                study.init = c(study.init[1:(ind.all.partial-1)], levels(object$study))
            } else if (ind.all.partial == 1 & ind.all.partial == length(study.init)) {
                # there's only "all.partial"
                study.init = levels(object$study)
            }
        }
        
        study.init = unique(study.init) #once again cause we added studies if "all.partial"
        
        if (!isNULL(subtitle))
        {
            if (length(subtitle)!=length(study.init)| length(subtitle)!=length(unique(subtitle)))
                stop("'subtitle' indicates the subtitle of the plot for each study and it needs to be the same length as 'study' (", length(study.init),") and duplicate are not allowed. 'study' includes: ", paste(study.init, collapse = ", "))
        }
        
        df.final = data.frame()
        
        indice.all = grep("global", study.init) # can go faster before and after "global"
        if (length(indice.all)>0)
        {
            study.list = list()
            i = 1
            if (indice.all>1)
            {
                study.list[[1]] = study.init[1:(indice.all-1)]
                i = i+1
            }
            
            study.list[[i]] = study.init[indice.all]
            
            if (indice.all<length(study.init))
                study.list[[i+1]] = study.init[-(1:indice.all)]
        } else {
            study.list = list(study.init)
        }
        
        
        # the following loop consider subset of studies all together, up until "global", and subset of studies after "global"
        for (length.study in 1 : length(study.list))
        {
            object = object.init #reinitialise $variates
            study = study.list[[length.study]]
            
            #-- define 'blocks'
            if (any(study == "global"))
            {
                # can plot both X and Y when one study or when study="global"
                # same as class.object==pls
                
                if (rep.space == "multi")
                {
                    blocks = c("X", "Y")
                    object$variates = object$variates[names(object$variates) %in% blocks]
                }
                
                if (rep.space == "X-variate")
                {
                    object$variates = object$variates["X"]
                    blocks = "X"
                }
                
                if (rep.space == "Y-variate")
                {
                    object$variates = object$variates["Y"]
                    blocks = "Y"
                }
                
                if (rep.space == "XY-variate")
                {
                    object$variates$XYvariates = (object$variates$X + object$variates$Y)/2
                    object$variates = object$variates["XYvariates"]
                    blocks = "XY combined"
                }
                
            } else if (length(study) == 1) {
                # can plot only X, Y or XY variate when more than one study
                # can plot both X and Y when one study or when study="global"
                
                blocks = c("X", "Y")
                
                if (rep.space == "X-variate")
                    blocks = "X"
                
                if (rep.space == "Y-variate")
                    blocks = "Y"
                
                #extract variates for each "blocks" for "study"
                object$variates = lapply(object$variates.partial, function(x){x[[study]]})[names(object$variates) %in% blocks]
                
                #if XY-variate, combine the previous variates (relative to "blocks" and "study")
                if (rep.space == "XY-variate")
                {
                    object$variates$XYvariates = (object$variates$X + object$variates$Y)/2
                    object$variates = object$variates["XYvariates"]
                    blocks = "XY combined"
                }
                blocks.init = blocks #save for ".getVariatesAndLabels"
                blocks = study
                
                
            } else { #length(study)>1
                
                blocks = c("X", "Y")
                
                if (rep.space == "multi")
                {
                    rep.space = "X-variate"
                    warning("More than one study is plotted, 'rep.space' is set to 'X-variate'. Alternatively, you can input 'Y-variate'")
                }
                
                if (rep.space == "X-variate")
                    blocks = "X"
                
                if (rep.space == "Y-variate")
                    blocks = "Y"
                
                
                #extract variates for each "blocks" for "study"
                object$variates = lapply(object$variates.partial, function(x)
                {
                    out = lapply(study, function(y){x[[y]]})
                    names(out) = study
                    out
                })[names(object$variates) %in% blocks]#[[1]]
                
                #if XY-variate, combine the previous variates (relative to "blocks" and "study")
                if (rep.space == "XY-variate")
                {
                    for (i in 1:length(object$variates$X))
                        object$variates$XYvariates[[i]] = (object$variates$X[[i]]+object$variates$Y[[i]])/2
                    
                    names(object$variates$XYvariates) = names(object$variates$X)
                    object$variates = object$variates[["XYvariates"]]
                } else {
                    object$variates = object$variates[[1]] # get rid of the $X or $Y
                }
                
                # blocks becomes study, so each study is plotted
                blocks = study
                object$names$sample = lapply(object$variates, rownames)
                ellipse = FALSE
                star = FALSE
                centroid = FALSE
            }
            
            #-- check inputs
            # check style as we do not do 3d at the moment:
            if (!style %in% c("ggplot2", "lattice", "graphics"))
                stop("'style' must be one of 'ggplot2', 'lattice' or 'graphics'.", call. = FALSE)
            
            check = .plotIndivCheckInput(object = object, comp = comp, blocks = blocks, ind.names = ind.names,
                                         style = style, ellipse = ellipse, ellipse.level = ellipse.level, centroid = centroid,
                                         star = star, legend = legend, X.label = X.label, Y.label = Y.label, abline = abline,
                                         xlim = xlim, ylim = ylim, plot_parameters = plot_parameters)
            #-- retrieve some outputs from the checks
            comp = check$comp
            xlim = check$xlim
            ylim = check$ylim
            ind.names = check$ind.names
            display.names = FALSE#check$display.names
            
            
            #-- get the variates
            variate = .getVariatesAndLabels(object, comp, blocks.init = blocks.init, blocks = blocks, rep.space = rep.space,
                                            style = style, X.label = X.label, Y.label = Y.label, Z.label = NULL)
            #-- retrieve outputs
            x = variate$x
            y = variate$y
            z = variate$z
            X.label = variate$X.label #only the last one of the loop is used
            Y.label = variate$Y.label #only the last one of the loop is used
            
            n = nrow(object$X)
            
            # create data frame df that contains (almost) all the ploting information
            out = .inputShapePlotIndiv(object = object, n = n, blocks = blocks, x = x, y = y, z = z, ind.names = ind.names, group = group,
                                       col.per.group = col.per.group, style = style, study = study, ellipse = ellipse, ellipse.level = ellipse.level,
                                       centroid = centroid, star = star, title = title, xlim = xlim, ylim = ylim,
                                       col = col, cex = cex, pch = pch, display.names = display.names, plot_parameters = plot_parameters)
            #-- retrieve outputs
            df = out$df
            df.ellipse = out$df.ellipse
            col.per.group = out$col.per.group
            title = out$title
            display.names = out$display.names
            xlim = out$xlim
            ylim = out$ylim
            #missing.col = out$missing.col
            plot_parameters = out$plot_parameters
            
            #save(list=ls(),file="temp.Rdata")
            # concatenate results
            df.final = rbind(df.final, df)
        }
        # add study information on df.final, for pch legend
        study.levels = study.init[which(!study.init == "global")]
        if (any(study.init == "global"))
            study.levels = levels(object$study)
        
        
        # change the levels of df.final$Block to "subtitle"
        if (!isNULL(subtitle))
        {
            df.final$Block = factor(df.final$Block, labels = subtitle)
            
            if(ellipse)
                df.ellipse$Block = factor(df.ellipse$Block, labels = subtitle)
        }
        df = df.final
        
        if (style == "ggplot2")
            style = "ggplot2-MINT"
        
        #call plot module (ggplot2, lattice, graphics, 3d)
        res = .graphicModule(df = df, centroid = centroid, col.per.group = col.per.group, title = title,
                             X.label = X.label, Y.label = Y.label, xlim = xlim, ylim = ylim, class.object = class(object),
                             display.names = display.names, legend = legend, abline = abline,
                             star = star, ellipse = ellipse, df.ellipse = df.ellipse, style = style, layout = layout,
                             #missing.col = missing.col,
                             #for ggplot2-MINT
                             study.levels = study.levels, plot_parameters = plot_parameters
        )
        
        return(invisible(list(df = df, graph = res)))
    }

#' @export
#' @method plotIndiv mint.spls
#' @rdname plotIndiv
plotIndiv.mint.spls  <- plotIndiv.mint.pls

#' @export
#' @method plotIndiv mint.plsda
#' @rdname plotIndiv
plotIndiv.mint.plsda  <- plotIndiv.mint.pls

#' @export
#' @method plotIndiv mint.splsda
#' @rdname plotIndiv
plotIndiv.mint.splsda  <- plotIndiv.mint.pls

## ----------- CCA family ----------- 

#----------------------------------------------------------------------------------------------------------#
#-- Includes plotIndiv rGCCA, sGCCA, sGCCDA --#
#----------------------------------------------------------------------------------------------------------#

#' @export
#' @method plotIndiv mint.sgcca
#' @rdname plotIndiv
plotIndiv.sgcca <-  function(object,
                             comp = NULL,
                             blocks = NULL,
                             # to choose which block data to plot, when using GCCA module
                             ind.names = TRUE,
                             group,
                             # factor indicating the group membership for each sample, useful for ellipse plots. Coded as default for the -da methods, but needs to be input for the unsupervised methods (PCA, IPCA...)
                             col.per.group,
                             style = "ggplot2",
                             # can choose between graphics, 3d, lattice or ggplot2
                             ellipse = FALSE,
                             ellipse.level = 0.95,
                             centroid = FALSE,
                             star = FALSE,
                             title = NULL,
                             subtitle,
                             legend = FALSE,
                             X.label = NULL,
                             Y.label = NULL,
                             Z.label = NULL,
                             abline = FALSE,
                             xlim = NULL,
                             ylim = NULL,
                             col,
                             cex,
                             pch,
                             pch.levels,
                             alpha = 0.2,
                             axes.box = "box",
                             layout = NULL,
                             size.title = rel(2),
                             size.subtitle = rel(1.5),
                             size.xlabel = rel(1),
                             size.ylabel = rel(1),
                             size.axis = rel(0.8),
                             size.legend = rel(1),
                             size.legend.title = rel(1.1),
                             legend.title = "Legend",
                             legend.title.pch = "Legend",
                             legend.position = "right",
                             point.lwd = 1
                             
)
{
    plot_parameters = list(size.title = size.title, size.subtitle = size.subtitle, size.xlabel = size.xlabel, size.ylabel = size.ylabel,
                           size.axis = size.axis, size.legend = size.legend, size.legend.title = size.legend.title, legend.title = legend.title,
                           legend.title.pch = legend.title.pch, legend.position = legend.position, point.lwd = point.lwd, alpha = alpha)
    
    if(any(class(object)%in%c("mint.block.pls", "mint.block.spls", "mint.block.plsda", "mint.block.splsda")))
        stop("No plotIndiv for the following functions at this stage: mint.block.pls, mint.block.spls, mint.block.plsda, mint.block.splsda.")
    
    #-- rep.space
    rep.space = "multi" # rep.space is not used afterwards, put to "multi" to plot all blocks
    
    
    if (is.null(blocks))
    {
        blocks = names(object$X)#names$blocks
    } else if (is.numeric(blocks) & min(blocks) > 0 &  max(blocks) <=  length(object$names$blocks)) {
        blocks = object$names$blocks[blocks]
    } else if (is.character(blocks)) {
        if (!any(blocks %in% object$names$blocks))
            stop("One element of 'blocks' does not match with the names of the blocks")
    } else {
        stop("Incorrect value for 'blocks'", call. = FALSE)
    }
    #object$variates = object$variates[names(object$variates) %in% blocks] # reduce the variate to the 'blocks' we are looking at
    object$variates = object$variates[match(blocks, names(object$variates))] # reduce the variate to the 'blocks' we are looking at
    
    if (any(object$ncomp[blocks] ==  1))
        stop(paste("The number of components for one selected block '", paste(blocks, collapse = " - "), "' is 1. The number of components must be superior or equal to 2."), call. = FALSE)
    ncomp = object$ncomp[blocks]
    
    
    if(length(blocks)!= length(unique(blocks)))
        stop("Duplicate in 'blocks' not allowed")
    
    if (!isNULL(subtitle))
    {
        if(length(subtitle)!= length(blocks) | length(subtitle)!= length(unique(subtitle)))
            stop("'subtitle' indicates the subtitle of the plot for each 'blocks'; it needs to be the same length as 'blocks' and duplicate are not allowed.")
    }
    
    
    
    #-- check inputs
    check = .plotIndivCheckInput(object = object, comp = comp, blocks = blocks, ind.names = ind.names,
                                 style = style, ellipse = ellipse, ellipse.level = ellipse.level, centroid = centroid,
                                 star = star, legend = legend, X.label = X.label, Y.label = Y.label, Z.label = Z.label, abline = abline,
                                 xlim = xlim, ylim = ylim, alpha = alpha, axes.box = axes.box, plot_parameters = plot_parameters)
    #-- retrieve outputs from the checks
    axes.box = check$axes.box
    comp = check$comp
    xlim = check$xlim
    ylim = check$ylim
    ind.names = check$ind.names
    display.names = check$display.names
    
    
    #-- get the variates
    variate = .getVariatesAndLabels(object, comp, blocks = blocks, style = style, X.label = X.label, Y.label = Y.label, Z.label = Z.label, rep.space = rep.space)
    #-- retrieve outputs
    x = variate$x
    y = variate$y
    z = variate$z
    X.label = variate$X.label
    Y.label = variate$Y.label
    Z.label = variate$Z.label
    
    n = nrow(object$X[[1]])
    
    # create data frame df that contains (almost) all the ploting information
    out = .inputShapePlotIndiv(object = object, n = n, blocks = blocks, x = x, y = y, z = z, ind.names = ind.names, group, col.per.group = col.per.group,
                               style = style, study = "global", ellipse = ellipse, ellipse.level = ellipse.level,
                               centroid = centroid, star = star, title = title, xlim = xlim, ylim = ylim,
                               col = col, cex = cex, pch = pch, pch.levels = pch.levels, display.names = display.names, plot_parameters = plot_parameters)
    #-- retrieve outputs
    df = out$df
    df.ellipse = out$df.ellipse
    col.per.group = out$col.per.group
    title = out$title
    display.names = out$display.names
    xlim = out$xlim
    ylim = out$ylim
    #missing.col = out$missing.col
    ellipse = out$ellipse
    centroid = out$centroid
    star = out$star
    plot_parameters = out$plot_parameters
    
    
    # change the levels of df.final$Block to "subtitle"
    if (!isNULL(subtitle))
    {
        df$Block = factor(df$Block, labels = subtitle)
        if(ellipse)
            df.ellipse$Block = factor(df.ellipse$Block, labels = subtitle)
    }
    
    #call plot module (ggplot2, lattice, graphics, 3d)
    res = .graphicModule(df = df, centroid = centroid, col.per.group = col.per.group, title = title, X.label = X.label,
                         Y.label = Y.label, Z.label = Z.label, xlim = xlim, ylim = ylim, class.object = class(object), display.names = display.names, legend = legend,
                         abline = abline, star = star, ellipse = ellipse, df.ellipse = df.ellipse, style = style, layout = layout,
                         #missing.col = missing.col,
                         axes.box = axes.box, plot_parameters = plot_parameters, alpha = alpha)
    
    return(invisible(list(df = df, df.ellipse = df.ellipse, graph = res)))
    
}

#' @export
#' @method plotIndiv mint.rgcca
#' @rdname plotIndiv
plotIndiv.rgcca <- plotIndiv.sgcca