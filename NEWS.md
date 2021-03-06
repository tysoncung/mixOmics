## version 6.14.0

### new features / enhancements / changes

* `circosPlot`: The radial location of feature names can now be cutomised using `var.adj`

### bug fixes

* `plotIndiv`: Legend bug which misspecified the groups resolved

-------------------------------------------------------------------------------
## version 6.12.0

### new features / enhancements

* `plotLoadings`'s infamous *figure margins too large* error now handled and informative condition thrown
* `circosPlot`'s `lines` argument default to `FALSE` now
* `circosPlot`'s inconsistentcy of blocks with identical `X` names fixed
*  consensus and weighted consensus plots now supported for `plotIndiv` with relevant block analyses
* `plotLoadings`'s feature name trimming can be customised
* `block.splsda` bug which could drop some Y factors with `near.zero.variance=TRUE` fixed
* `perf.block.splsda` now supports calculation of combined and per-block AUC
* model improvement significance can be custmoised in all `perf` and `tune` functions
* `perf.block.splsda` is now much faster and supports FORK clusters
* `tune.(s)pls(da)`, `perf.(s)plsda` now support FORK clusters

### bug fixes

* `circosPlot`'s faded `lines` bug when many `NA`s present fixed
* `tune.block.splsda()` bug when using fixed `test.keepX` over two or more blocks fixed
* `circosPlot` and `plotLoadings` bug caused by features with `NA`s fixed
* `plotIndiv(..., ind.names = FALSE)` warning/bug fixed
* `tune.block.splsda` bug on Windows parallelisation fixed
* `perf` and `tune` functions' issue  when choosing the optimum component resolved
* added option to suppress `auroc` from printing all the AUCs
-------------------------------------------------------------------------------
## version 6.10.0

### new features / enhancements

* parallel processing on `tune.block.splsda` improved
* `tune.block.splsda` now supports more distances
* You can now customise `auroc` plots. Refer to documentation for more info

### bug fixes

* single factor multilevel error in `pls` fixed
* fixed over-estimated correlation of `cim` for `mixo_(s)pls` objects with single component 
* margin error in `cim` now handled properly
* fixed `plotLoadings` error for very long variable names
* `predict` function bug for single sample prediction fixed
* `plotLoadings` bug for long variable names fixed
* Fixed `tune.spls` and `pef.plsda` bugs when using `cpus` argument for parallel 
processing
* `perf.plot` bug in extracting names fixed
* Few fixes for `tune.splsda` with AUC

### minor improvements

* missing values in `plotIndiv`'s `group` argument no more throws error
* `mixOmics::predict` function documentation now more accessible
* names of `linnerud` datasets fixed
* `plot.perf` now respects `ylim` arguments for custom y range
* package startup message with direct liks to useful resources
* `mixOmics` function documentation disambiguated with instruction on how to get
package help
* Updated onLoad message with discussion forum info, bug reports, and more
* Dropped legacy `comp.tol` argument from `pca`
* `plot.perf` now respects `ylim` arguments for custom y range
* Added Code of Conduct

-------------------------------------------------------------------------------
## version 6.8.0

* NOW HOSTED ON BIOCONDUCTOR

-------------------------------------------------------------------------------

## version < 6.8.0

* Refer to `./inst/legacy/NEWS-old` on GitHub repo