

<!-- README.md is generated from README.qmd. Please edit that file -->

# vsn

**V**ariance **s**tabilization and calibration for microarray data.

`vsn` implements a method for normalising single- and multiple-color
microarray intensities (and, in principle, data from other technologies
with a similar format). The model incorporates data calibration step
(a.k.a. normalization), a model for the dependence of the variance on
the mean intensity and a variance stabilizing data transformation.
Differences between transformed intensities are analogous to “normalized
log-ratios”. However, in contrast to the latter, their variance is
independent of the mean, and they are usually more sensitive and
specific in detecting differential transcription.

## Installation

`vsn` is part of Bioconductor:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("vsn")
```

## Usage

``` r
library("vsn")
```

The simplest way to normalize a dataset is:

``` r
data("kidney")          # ExpressionSet of unnormalised data
xnorm <- justvsn(kidney)
xnorm
#> ExpressionSet (storageMode: lockedEnvironment)
#> assayData: 8704 features, 2 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData
#>   sampleNames: green red
#>   varLabels: channel
#>   varMetadata: labelDescription
#> featureData: none
#> experimentData: use 'experimentData(object)'
#> Annotation:
```

For more control, fit the model and apply it in two steps:

``` r
fit <- vsn2(kidney)
ynorm <- predict(fit, kidney)
ynorm
#> ExpressionSet (storageMode: lockedEnvironment)
#> assayData: 8704 features, 2 samples 
#>   element names: exprs 
#> protocolData: none
#> phenoData
#>   sampleNames: green red
#>   varLabels: channel
#>   varMetadata: labelDescription
#> featureData: none
#> experimentData: use 'experimentData(object)'
#> Annotation:
```

Both are equivalent. The two-step form is useful when fitting on a
subset (e.g. spike-in or control features) and then applying the fit to
the full data, or when you want to inspect the `fit` object. `justvsn()`
and `vsn2()` are also available for `AffyBatch` and `RGList` objects.

## References

``` r
citation(package = "vsn")
#> To cite the vsn package in publications use:
#> 
#>   Huber W, von Heydebreck A, Sueltmann H, Poustka A, Vingron M (2002).
#>   "Variance Stabilization Applied to Microarray Data Calibration and to
#>   the Quantification of Differential Expression." _Bioinformatics_, *18
#>   Suppl. 1*, S96-S104. doi:10.1093/bioinformatics/18.suppl_1.s96
#>   <https://doi.org/10.1093/bioinformatics/18.suppl_1.s96>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Variance Stabilization Applied to Microarray Data Calibration and to the Quantification of Differential Expression},
#>     author = {Wolfgang Huber and Anja {von Heydebreck} and Holger Sueltmann and Annemarie Poustka and Martin Vingron},
#>     doi = {10.1093/bioinformatics/18.suppl_1.s96},
#>     journal = {Bioinformatics},
#>     year = {2002},
#>     volume = {18 Suppl. 1},
#>     pages = {S96-S104},
#>   }
```

## Feedback

This is an approved de.NBI service. Please help us improve by taking our
[short user
survey](https://www.surveymonkey.de/r/denbi-service?sc=hd-hub&tool=vsn).
