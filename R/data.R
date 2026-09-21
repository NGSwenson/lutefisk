#' Ten-species example system
#'
#' A community data matrix of five communities and ten species, a
#' phylogeny of the ten species, and three simulated continuous traits.
#' This is the toy system used throughout chapters 3 to 8 of the book. The
#' same objects are also available under the chapter-specific names used
#' in the text (\code{beta.example.*}, \code{null.example.*},
#' \code{partition.example.*}, \code{comparative.*}).
#'
#' @format \code{lefse.sample} is a data frame with 5 rows (communities)
#'   and 10 columns (species abundances). \code{lefse.phylo} is a
#'   \code{phylo} object with 10 tips. \code{lefse.traits} is a data frame
#'   with 10 rows (species) and 3 numeric columns (\code{trait1} to
#'   \code{trait3}).
#' @source Simulated for the first edition of the book.
#' @name lefse.sample
#' @aliases lefse.phylo lefse.traits beta.example.sample beta.example.phylo
#'   beta.example.traits null.example.sample null.example.phylo
#'   null.example.traits partition.example.sample partition.example.phylo
#'   partition.example.traits comparative.phylo comparative.traits
#' @keywords datasets
NULL

#' Example phylogenies for tree manipulation
#'
#' \code{my.phylo} is a 26-tip ultrametric tree with tips \code{a} to
#' \code{z} and labeled internal nodes, used in chapter 2 to demonstrate
#' reading, plotting, pruning and manipulating trees. \code{example.poly} is
#' a small tree containing polytomies and named genera.
#'
#' @format \code{phylo} objects.
#' @name my.phylo
#' @aliases example.poly
#' @keywords datasets
NULL

#' Phylogenetic diversity example
#'
#' Five communities of fifteen species (\code{spA} to \code{spO}) and their
#' phylogeny, used in chapter 3.
#'
#' @format \code{PD.example.sample} is a data frame with 5 rows and 15
#'   columns. \code{PD.example.phylo} is a \code{phylo} object with 15 tips.
#' @name PD.example.sample
#' @aliases PD.example.phylo
#' @keywords datasets
NULL

#' Functional diversity example
#'
#' Five communities of twelve species and five traits with different
#' scales, used in chapter 4 to illustrate trait standardization and
#' functional diversity metrics.
#'
#' @format \code{FD.example.sample} is a data frame with 5 rows and 12
#'   columns. \code{FD.traits} is a data frame with 12 rows and 5 numeric
#'   columns.
#' @name FD.example.sample
#' @aliases FD.traits
#' @keywords datasets
NULL

#' Null model example with a forty-tip tree
#'
#' A 40-tip simulated phylogeny and presence-only community data in
#' three-column (community, abundance, species) form, used in chapter 6 to
#' illustrate null models for phylogenetic diversity.
#'
#' @format \code{null.pd.example.phylo} is a \code{phylo} object.
#'   \code{null.pd.example.sample} is a data frame with columns
#'   \code{community}, \code{abundance} and \code{species}.
#' @name null.pd.example.phylo
#' @aliases null.pd.example.sample
#' @keywords datasets
NULL

#' Environmental and spatial data for the example communities
#'
#' Simulated environmental variables and spatial coordinates for the
#' communities of the ten-species example system, used in chapter 8.
#'
#' @format \code{env.data} has three numeric columns (\code{env1} to
#'   \code{env3}); \code{xy.data}, \code{xy.beta.data} and
#'   \code{new.my.xys} each hold coordinates in columns \code{x} and
#'   \code{y}. All have communities as row names.
#' @name env.data
#' @aliases xy.data xy.beta.data new.my.xys
#' @keywords datasets
NULL

#' lutefisk: Functional and Phylogenetic Ecology in R
#'
#' Companion package to the second edition of "Functional and Phylogenetic
#' Ecology in R". The functions are the hand-built implementations
#' developed in the book, written in base R so that the arithmetic behind
#' each metric is visible, and vectorized with matrix algebra wherever the
#' mathematics allows. The package succeeds \pkg{lefse} 0.5.
#'
#' @section Alpha diversity:
#' \code{\link{Fmpd}}, \code{\link{Fmntd}}, \code{\link{faith.pd}},
#' \code{\link{weighted.faith}}, \code{\link{com.trait.moments}},
#' \code{\link{trait.range}}, \code{\link{hull.volume}},
#' \code{\link{rank.abund.pd}}.
#' @section Beta diversity:
#' \code{\link{new.comdist}} and relatives, \code{\link{shuffle.dist}}.
#' @section Lineage representation:
#' \code{\link{lineage.hyper}}.
#' @section Utilities:
#' \code{\link{branch.membership}}, \code{\link{to.pa}},
#' \code{\link{weighted.sd}}.
#' @docType package
#' @name lutefisk-package
#' @aliases lutefisk
#' @keywords internal
"_PACKAGE"

#' A larger simulated gradient system for comparing metrics
#'
#' One hundred species on a pure-birth phylogeny, three traits, and forty
#' communities arrayed along an environmental gradient. Trait one evolved
#' under Brownian motion and sets each species' environmental optimum; trait
#' two is weakly conserved (Pagel's lambda of 0.3); trait three is unrelated
#' to the phylogeny. Community abundances are Poisson draws whose means fall
#' off with the distance between a species' trait one value and the
#' community's environment, with a filter width that varies among
#' communities, and with a secondary filter on trait three. The result is a
#' set of communities that range from strongly phylogenetically clustered to
#' overdispersed, which is what one wants for comparing metrics and null
#' models. The script that built the data is in the package's
#' \code{data-raw} directory.
#'
#' @format \code{gradient.sample} is a 40 by 100 data frame of abundances;
#'   \code{gradient.phylo} is a \code{phylo} object with 100 tips;
#'   \code{gradient.traits} is a 100 by 3 data frame; \code{gradient.env}
#'   is a 40 by 4 data frame with coordinates \code{x} and \code{y} and
#'   environmental variables \code{env1} and \code{env2}.
#' @source Simulated for the second edition of the book.
#' @name gradient.sample
#' @aliases gradient.phylo gradient.traits gradient.env
#' @keywords datasets
NULL
