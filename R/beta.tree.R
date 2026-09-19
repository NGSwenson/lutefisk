#' Tree-based phylogenetic beta diversity: UniFrac and PhyloSor
#'
#' Unweighted UniFrac (the fraction of the branch length spanned by two
#' communities that is unique to one of them) and PhyloSor (twice the
#' shared branch length divided by the summed branch lengths of the two
#' communities) for every pair of communities, together with the
#' abundance-weighted UniFrac of Lozupone et al. (2007). The three are
#' computed for all pairs at once from the tip-by-branch membership matrix.
#'
#' Let \eqn{B} be the community-by-branch matrix of ones and zeros in which
#' an entry is one when at least one species of the community descends from
#' the branch, and \eqn{l} the vector of branch lengths. Then \eqn{B \,
#' \mathrm{diag}(l) \, B'} holds, for every pair of communities, the branch
#' length they share, and its diagonal holds each community's Faith's
#' index (with the root path). UniFrac and PhyloSor are simple functions of
#' those quantities. For the weighted version \eqn{B} is replaced by the
#' community-by-branch matrix of summed relative abundances under each
#' branch, and the result is the sum over branches of the branch length
#' times the absolute difference in those abundances, normalized by the
#' abundance-weighted root-to-tip distance.
#'
#' @inheritParams faith.pd
#' @return A community-by-community matrix (a \code{dist} object can be
#'   made with \code{as.dist}).
#' @references Lozupone, C. and Knight, R. (2005) UniFrac: a new
#'   phylogenetic method for comparing microbial communities. Applied and
#'   Environmental Microbiology 71, 8228-8235. Bryant, J.A. et al. (2008)
#'   Microbes on mountainsides: contrasting elevational patterns of
#'   bacterial and plant diversity. PNAS 105, 11505-11511.
#' @examples
#' data(beta.example.sample, beta.example.phylo)
#' unifrac.pairs(beta.example.phylo, beta.example.sample)
#' phylosor.pairs(beta.example.phylo, beta.example.sample)
#' unifrac.weighted(beta.example.phylo, beta.example.sample)
#' @export
unifrac.pairs <- function(my.phylo, my.sample) {
	sb <- shared.branch.length(my.phylo, my.sample)
	pd <- diag(sb)
	total <- outer(pd, pd, "+") - sb
	out <- (total - sb) / total
	diag(out) <- 0
	out
}

#' @rdname unifrac.pairs
#' @export
phylosor.pairs <- function(my.phylo, my.sample) {
	sb <- shared.branch.length(my.phylo, my.sample)
	pd <- diag(sb)
	2 * sb / outer(pd, pd, "+")
}

#' @rdname unifrac.pairs
#' @export
unifrac.weighted <- function(my.phylo, my.sample) {
	ra <- to.ra(my.sample)
	if (!all(colnames(ra) %in% my.phylo$tip.label)) {
		stop("all species in my.sample must be tips of my.phylo")
	}
	m <- branch.membership(my.phylo)[colnames(ra), , drop = FALSE]
	bl <- my.phylo$edge.length
	## relative abundance under each branch, communities by branches
	a <- ra %*% m
	n <- nrow(a)
	out <- matrix(0, n, n, dimnames = list(rownames(ra), rownames(ra)))
	## root-to-tip distance of each species
	depth <- ape::node.depth.edgelength(my.phylo)[match(colnames(ra), my.phylo$tip.label)]
	for (i in seq_len(n)) for (j in seq_len(n)) if (j > i) {
		u <- sum(bl * abs(a[i, ] - a[j, ]))
		D <- sum(depth * (ra[i, ] + ra[j, ]))
		out[i, j] <- out[j, i] <- u / D
	}
	out
}

## Community-by-community matrix of shared branch length, whose diagonal is
## Faith's index including the root path.
shared.branch.length <- function(my.phylo, my.sample) {
	pa <- to.pa(my.sample)
	if (!all(colnames(pa) %in% my.phylo$tip.label)) {
		stop("all species in my.sample must be tips of my.phylo")
	}
	m <- branch.membership(my.phylo)[colnames(pa), , drop = FALSE]
	on.path <- ((pa %*% m) > 0) * 1
	(on.path * rep(my.phylo$edge.length, each = nrow(on.path))) %*% t(on.path)
}
