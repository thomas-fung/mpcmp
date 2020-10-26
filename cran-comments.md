## Test environments
* local R installation, R 4.0.2
* ubuntu 16.04 (on travis-ci), R 4.0.3
* ubuntu 16.04, Windows, macOS (on GitHub Actions), R 4.0.3
* win-builder (devel and release) 

This is a patch submission to fix some issues picked up by CRAN checks.

## R CMD check results 

0 errors | 0 warnings | 0 notes

Errors and notes from previous submission.

Version: 0.3.5 
Check: dependencies in R code 
Result: NOTE 
    Namespace in Imports field not imported from: ‘broom’
     All declared Imports should be used. 
Flavors: r-devel-linux-x86_64-fedora-clang, r-devel-linux-x86_64-fedora-gcc, r-release-macos-x86_64, r-oldrel-macos-x86_64

* Moved `broom` to the `Suggests` field in Description. 

Version: 0.3.5 
Check: whether package can be installed 
Result: ERROR 
    Installation failed. 
Flavor: r-patched-solaris-x86

* Cast integer to double for log function in Rcpp code. This should fix the solaris issue. 

## Downstream dependency. 

There are currently no downstream dependencies for this package.

