## Test environments
* local R installation, R 4.2.0
* ubuntu 16.04, Windows, macOS (on GitHub Actions), R (devel, release, oldrel)
* win-builder (devel and release) 
* rhub:check (debian-clang-devel, debian-gcc-devel)

## R CMD check results 

0 errors | 0 warnings | 1 notes

Version: 0.3.8
Check: New submission. Package was archived on CRAN.
Result: Note 
* As the old version was achieved, this is a new submission but it's a also a patch submission to fix up some issues picked up by CRAN checks.

Check: Potential misspelled word in DESCRIPTION
    
* Parametrized is an alternative spelling to parameterized. 

There is also an error in running R-Hub checks which showed Bioconductor does not yet build and check packages for R version 4.3 but that should be irrelevant to our package.

## Downstream dependency. 

There are currently no downstream dependencies for this package.