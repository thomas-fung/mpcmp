# Sitophilus data set

Ribeiro et al. (2013) carried out an experiment to assess the
bioactivity of extracts from different parts (seeds, leaves and
branches) of Annona mucosa (Annonaceae) to control Sitophilus zeamais
(Coleoptera: Curculionidae), a major pest of stored maize/corn in
Brazil.

## Usage

``` r
data(sitophilus)
```

## Format

A data frame with 40 observations on 2 variables.

- extract:

  the treatment used

- ninsect:

  number emerged insects (progeny)

## Value

A data frame with 40 rows and 2 variables.

## Details

10g of corn and 20 animals adults were placed in each Petri dish.
Extracts prepared with different parts of mucosa or just water (control)
were completely randomized with 10 replicates.

The numbers of emerged insects (progeny) after 60 days and their
corresponding treatments were recorded in this dataset.

This dataset was obtained from the \`cmpreg“package of Ribeiro Jr,
Zeviani & Demétrio (2019), which is based on the work of Ribeiro Junior
et al. (2019).

This data set is also used to illustrate the syntax for regression on
the dispersion.

## References

Ribeiro Junior, E.E., Zeviani, W.M., Bonat, W.H., Demétrio, C.G., &
Hinde, J. (2019). Reparametrization of COM–Poisson regression models
with applications in the analysis of experimental data. Statistical
Modelling. https://doi.org/10.1177%2F1471082X19838651.

Ribeiro, L.P., Vendramim, J.D., Bicalho, K.U., Andrade, M.S., Fernandes,
J.B., Moral, R.A., & Demétrio C.G.B. (2013). Annona mucosa Jacq.
(Annonaceae): A promising source of bioactive compounds against
Sitophilus zeamais Mots. (Coleoptera: Curculionidae). *Journal of Stored
Products Research*, **55**, 6-14.

## Examples

``` r
## For examples see example(glm.cmp)
```
