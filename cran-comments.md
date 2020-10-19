## Test environments
* local R installation, R 4.0.2
* ubuntu 16.04 (on travis-ci), R 4.0.2
* ubuntu 16.04, Windows, macOS (on GitHub Actions), R 4.0.2
* win-builder (devel and release)

## R CMD check results

There were no ERRORs or WARNINGs. 

There was 1 NOTE:

* checking CRAN incoming feasibility ... NOTE
  Maintainer: ‘Thomas Fung <thomas.fung.dr@gmail.com>’
  
  New maintainer:
    Thomas Fung <thomas.fung.dr@gmail.com>
  Old maintainer(s):
    Thomas Fung <thomas.fung@mq.edu.au>

The maintainer changed his email address. 

There are some issues with the ubuntu 16.04 check via GitHub actions: 
"Failed to get R 4.0.3: Failed to get R 4.0.3" but it passes the checks on travis-ci. 

## Downstream dependency. 

There are currently no downstream dependencies for this package.



