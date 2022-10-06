## Resubmission

### Resolved problems:

This is a resubmission. In this version I have:

* Reduced the length of the title to less than 65 characters.

* Added the missing \\value{} part to "plot.cv.missoNet.Rd"

* Replaced \\dontrun{} with \\donttest{} in examples.
 
      Maintainer: the examples are executable but > 10 sec, so they were put in \\donttest{}.

* Ensured that no more than 2 cores are used in examples and vignettes.

* Added an argument for user-specified seed in R/missoNet.R, now the function does not set a seed to a specific number.

### Unresolved problems:

* If there are references describing the methods in your package, please add these in the description field of your DESCRIPTION file.

      Maintainer: the reference manuscript is still under construction and not yet publicly accessible, we will add it in the next release.


## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.
