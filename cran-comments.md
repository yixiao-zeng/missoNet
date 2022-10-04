## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Notes from other test environments 

I have also run R CMD check on other platforms, a total of 2 NOTEs were detected:

* Possibly misspelled words in DESCRIPTION:\
  covariate (16:53)\
  undirected (21:40)
  
  Author: common statistics terms.

* Examples with CPU (user + system) or elapsed time > 10s ...

  Author: some non-essential examples have been put into \\dontrun\{ \}

I also found that certain platforms with development version of R may report a error, because this package imports the 'ComplexHeatmap' package from Bioconductor. The development version of Bioconductor works with R version 4.2.0.
