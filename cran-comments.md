## Resubmission

This is a resubmission. In this version I have:

* Replaced possibly misspelled words in DESCRIPTION.

* Reduced the running time of the example check.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

I found that certain platforms with development version of R may report an preparation error, because this package imports the 'ComplexHeatmap' package from Bioconductor. The development version of Bioconductor works with R version 4.2.0.
