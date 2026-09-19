# lutefisk

Companion R package to the second edition of *Functional and Phylogenetic
Ecology in R* (Springer, Use R! series). Successor to `lefse` 0.5.

The functions are the hand-built, base-R implementations developed in the
book. Where the mathematics allows, the loops of the first edition have been
replaced by matrix algebra: the pairwise beta diversity matrices, for
example, are single matrix products and run several hundred times faster
than the nested loops they replace. Every function is tested against an
independent implementation (picante, FD, geometry, or the original lefse
0.5 code), which is the book's own advice carried into the test suite.

```r
# install.packages("remotes")
remotes::install_github("NGSwenson/lutefisk")
library(lutefisk)
data(lefse.sample, lefse.phylo)
Fmpd(cophenetic(lefse.phylo), lefse.sample)
```

The example datasets from the book are included both as R objects
(`data(package = "lutefisk")`) and as the original text files
(`system.file("extdata", package = "lutefisk")`).
