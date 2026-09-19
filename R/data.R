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
