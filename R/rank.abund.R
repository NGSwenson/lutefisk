#' Phylogenetic diversity accumulated across the rank abundance curve
#'
#' For each community, species are ordered from most to least abundant and
#' Faith's phylogenetic diversity is computed for the two most abundant
#' species, then the three most abundant, and so on. The first edition
#' called \code{picante::pd} once per rank inside a loop; here each branch
#' is assigned the rank at which its first descendant enters the community,
#' and the cumulative sum of branch lengths ordered by that rank gives the
#' whole curve in one step. The root path is included, matching
#' \code{picante::pd} with \code{include.root = TRUE}.
#'
#' @param my.sample A community data matrix with communities as rows and
#'   species as columns.
#' @param my.phylo A \code{phylo} object whose tips include every species in
#'   \code{my.sample}.
#' @return A list with one numeric vector per community giving PD at ranks
#'   2 through the community's species richness, or \code{NA} for
#'   communities with fewer than two species.
#' @examples
#' data(lefse.sample, lefse.phylo)
#' rank.abund.pd(lefse.sample, lefse.phylo)
#' @export
rank.abund.pd <- function(my.sample, my.phylo) {
	x <- as.matrix(my.sample)
	storage.mode(x) <- "numeric"
	m <- branch.membership(my.phylo)[colnames(x), , drop = FALSE]
	bl <- my.phylo$edge.length
	out <- lapply(seq_len(nrow(x)), function(i) {
		com <- x[i, ]
		spp <- names(sort(com[com > 0], decreasing = TRUE))
		if (length(spp) < 2) return(NA_real_)
		## rank at which each species enters
		rk <- seq_along(spp)
		## first rank among the descendants of each branch (Inf if none)
		sub <- m[spp, , drop = FALSE] * rk
		sub[sub == 0] <- Inf
		first <- apply(sub, 2, min)
		## cumulative branch length by entry rank
		by.rank <- tapply(bl[is.finite(first)], first[is.finite(first)], sum)
		pd.at <- numeric(length(spp))
		pd.at[as.integer(names(by.rank))] <- by.rank
		cumsum(pd.at)[-1]
	})
	names(out) <- rownames(x)
	out
}

#' @rdname rank.abund.pd
#' @export
rank.abund.pd.null <- function(my.sample, my.phylo) {
	rank.abund.pd(my.sample, phytools_free_tipShuffle(my.phylo))
}

## Shuffle tip labels on a phylogeny without depending on picante.
phytools_free_tipShuffle <- function(my.phylo) {
	my.phylo$tip.label <- sample(my.phylo$tip.label)
	my.phylo
}

#' Convex hull volume accumulated across the rank abundance curve
#'
#' As \code{\link{rank.abund.pd}} but for the convex hull volume of the
#' species in trait space, starting from the \code{ncol(traits) + 1} most
#' abundant species. The first edition called \code{FD::dbFD} at each rank,
#' which recomputes a Gower distance and principal coordinates every time;
#' here the hull is computed directly in the trait space supplied, so the
#' values are those of \code{\link{hull.volume}}.
#'
#' @inheritParams rank.abund.pd
#' @param my.traits A numeric matrix or data frame of traits with species
#'   as row names.
#' @return A list with one numeric vector per community.
#' @examples
#' data(lefse.sample, lefse.traits)
#' rank.abund.fric(lefse.sample, scale(lefse.traits))
#' @export
rank.abund.fric <- function(my.sample, my.traits) {
	if (!requireNamespace("geometry", quietly = TRUE)) {
		stop("the geometry package is required for rank.abund.fric")
	}
	x <- as.matrix(my.sample)
	storage.mode(x) <- "numeric"
	tr <- as.matrix(my.traits)
	storage.mode(tr) <- "numeric"
	k <- ncol(tr)
	out <- lapply(seq_len(nrow(x)), function(i) {
		com <- x[i, ]
		spp <- names(sort(com[com > 0], decreasing = TRUE))
		if (length(spp) < k + 1) return(NA_real_)
		sapply((k + 1):length(spp), function(r) {
			geometry::convhulln(tr[spp[1:r], , drop = FALSE], options = "FA")$vol
		})
	})
	names(out) <- rownames(x)
	out
}

#' @rdname rank.abund.fric
#' @export
rank.abund.fric.null <- function(my.sample, my.traits) {
	tr <- as.matrix(my.traits)
	rownames(tr) <- sample(rownames(tr))
	rank.abund.fric(my.sample, tr)
}
