#' Mean pairwise distance of communities
#'
#' The mean of the pairwise phylogenetic or functional distances between the
#' species present in each community, either unweighted (\code{Fmpd}) or
#' weighted by the product of the relative abundances of each pair
#' (\code{Fmpd.a}). Both are computed for every community at once with a
#' single matrix product rather than by looping over communities.
#'
#' For a presence-absence matrix \eqn{P} and distance matrix \eqn{D}, the
#' sum of all pairwise distances in community \eqn{i} is the \eqn{i}th
#' diagonal element of \eqn{P D P'}, and dividing by the number of pairs
#' \eqn{n_i (n_i - 1)} gives the mean. The abundance-weighted version uses
#' relative abundances \eqn{R} in place of \eqn{P} and needs no divisor,
#' because the weights already sum to one. This weighted form includes the
#' zero-distance self comparisons in the denominator, as does the
#' \code{mpd} function in \pkg{picante} with \code{abundance.weighted = TRUE}.
#'
#' @param dist.mat A species-by-species distance matrix (a \code{matrix} or
#'   \code{dist}) with dimnames matching the columns of \code{my.sample}.
#' @param my.sample A community data matrix with communities as rows and
#'   species as columns.
#' @return A named numeric vector with one value per community.
#' @examples
#' data(lefse.sample, lefse.phylo)
#' d <- cophenetic(lefse.phylo)
#' Fmpd(d, lefse.sample)
#' Fmpd.a(d, lefse.sample)
#' @export
Fmpd <- function(dist.mat, my.sample) {
	pa <- to.pa(my.sample)
	d <- match.dist(dist.mat, pa)
	n <- rowSums(pa)
	tot <- rowSums((pa %*% d) * pa)
	out <- tot / (n * (n - 1))
	out[n < 2] <- NA
	out
}

#' @rdname Fmpd
#' @export
Fmpd.a <- function(dist.mat, my.sample) {
	ra <- to.ra(my.sample)
	d <- match.dist(dist.mat, ra)
	rowSums((ra %*% d) * ra)
}

## Nearest-neighbor distances for every species in every community. Returns
## a communities-by-species matrix in which entry [i, s] is the distance from
## species s to its nearest neighbor among the other species present in
## community i, or NA when s is absent from i.
nn.within <- function(dist.mat, my.sample) {
	pa <- to.pa(my.sample)
	d <- match.dist(dist.mat, pa)
	diag(d) <- Inf
	out <- t(apply(pa, 1, function(p) {
		present <- p > 0
		if (sum(present) < 2) return(rep(NA_real_, length(p)))
		nn <- apply(d[present, present, drop = FALSE], 1, min)
		res <- rep(NA_real_, length(p))
		res[present] <- nn
		res
	}))
	dimnames(out) <- dimnames(pa)
	out
}

#' Nearest taxon distances of communities
#'
#' The mean (\code{Fmntd}, \code{Fmntd.a}) or standard deviation
#' (\code{Fsntd}, \code{sntd.a}) of the distance from each species present
#' in a community to its nearest neighbor in that community, unweighted or
#' weighted by relative abundance. The abundance-weighted standard deviation
#' uses \code{\link{weighted.sd}}.
#'
#' @inheritParams Fmpd
#' @return A named numeric vector with one value per community.
#' @examples
#' data(lefse.sample, lefse.phylo)
#' d <- cophenetic(lefse.phylo)
#' Fmntd(d, lefse.sample)
#' Fmntd.a(d, lefse.sample)
#' Fsntd(d, lefse.sample)
#' sntd.a(d, lefse.sample)
#' @export
Fmntd <- function(dist.mat, my.sample) {
	nn <- nn.within(dist.mat, my.sample)
	rowMeans(nn, na.rm = TRUE)
}

#' @rdname Fmntd
#' @export
Fmntd.a <- function(dist.mat, my.sample) {
	nn <- nn.within(dist.mat, my.sample)
	ra <- to.ra(my.sample)
	nn[is.na(nn)] <- 0
	rowSums(nn * ra)
}

#' @rdname Fmntd
#' @export
Fsntd <- function(dist.mat, my.sample) {
	nn <- nn.within(dist.mat, my.sample)
	apply(nn, 1, stats::sd, na.rm = TRUE)
}

#' @rdname Fmntd
#' @export
sntd.a <- function(dist.mat, my.sample) {
	nn <- nn.within(dist.mat, my.sample)
	ra <- to.ra(my.sample)
	sapply(seq_len(nrow(nn)), function(i) {
		ok <- !is.na(nn[i, ])
		weighted.sd(nn[i, ok], ra[i, ok])
	})
}

#' Faith's phylogenetic diversity and its abundance-weighted form
#'
#' \code{faith.pd} is the sum of the branch lengths connecting the species
#' in a community, computed for all communities at once as the product of
#' the presence-absence matrix with the tip-by-branch membership matrix. A
#' branch contributes its length when at least one of its descendants is
#' present. Including the branches leading from the community's most recent
#' common ancestor to the root of the full tree (\code{include.root = TRUE})
#' reproduces \code{picante::pd} with its default settings.
#'
#' \code{weighted.faith} is the abundance-weighted index of the first
#' edition: the number of branches in the community's pruned phylogeny
#' multiplied by the mean branch length weighted by the mean relative
#' abundance of the species each branch subtends. The first edition looped
#' over branches to obtain the mean abundances; here they come from one
#' matrix product per community.
#'
#' @param my.phylo A \code{phylo} object whose tips include every species in
#'   \code{my.sample}.
#' @param my.sample A community data matrix with communities as rows and
#'   species as columns.
#' @param include.root Logical. Include the path from the community's
#'   crown node to the root of \code{my.phylo}?
#' @return A named numeric vector with one value per community.
#' @examples
#' data(PD.example.sample, PD.example.phylo)
#' faith.pd(PD.example.phylo, PD.example.sample)
#' weighted.faith(PD.example.phylo, PD.example.sample)
#' @export
faith.pd <- function(my.phylo, my.sample, include.root = TRUE) {
	pa <- to.pa(my.sample)
	if (!all(colnames(pa) %in% my.phylo$tip.label)) {
		stop("all species in my.sample must be tips of my.phylo")
	}
	m <- branch.membership(my.phylo)[colnames(pa), , drop = FALSE]
	n.desc <- (pa %*% m)
	bl <- my.phylo$edge.length
	if (include.root) {
		## every branch with at least one present descendant
		out <- ((n.desc > 0) * 1) %*% bl
	} else {
		## drop branches that subtend all present species, which lie on the
		## path from the crown node to the root
		n <- rowSums(pa)
		keep <- (n.desc > 0) & (n.desc < n)
		out <- (keep * 1) %*% bl
	}
	out <- as.vector(out)
	names(out) <- rownames(pa)
	out
}

#' @rdname faith.pd
#' @export
weighted.faith <- function(my.phylo, my.sample) {
	x <- as.matrix(my.sample)
	storage.mode(x) <- "numeric"
	out <- apply(x, 1, function(com) {
		spp <- names(com)[com > 0]
		if (length(spp) < 2) return(NA_real_)
		tmp.tree <- ape::keep.tip(my.phylo, spp)
		m <- branch.membership(tmp.tree)
		## mean abundance of the species subtended by each branch
		a.i <- as.vector(com[rownames(m)] %*% m) / colSums(m)
		nrow(tmp.tree$edge) * sum(tmp.tree$edge.length * a.i) / sum(a.i)
	})
	names(out) <- rownames(x)
	out
}

#' Phylogenetic Hill numbers
#'
#' Phylogenetic diversity of order \eqn{q} in the attribute diversity
#' framework of Chao et al. (2010, 2014). For a tree of depth \eqn{T}, with
#' branch lengths \eqn{L_i} and \eqn{a_i} the summed relative abundance of
#' the species descending from branch \eqn{i},
#' \deqn{{}^qPD(T) = \left[ \sum_i L_i (a_i / T)^q \right]^{1/(1-q)}}
#' with the limit \eqn{\exp(-\sum_i L_i (a_i/T) \log(a_i/T))} at
#' \eqn{q = 1}. At \eqn{q = 0} this is Faith's index including the root
#' path; dividing by \eqn{T} gives the effective number of equally
#' distinct lineages, which is in the same units as species richness and
#' is what should be compared across trees of different depth.
#'
#' The computation is one matrix product of the relative abundance matrix
#' with the tip-by-branch membership matrix, so all communities and all
#' orders are handled at once.
#'
#' @inheritParams faith.pd
#' @param q A vector of orders; 0, 1 and 2 are the usual choices.
#' @param mean.pd Logical. Divide by the tree depth to return the
#'   effective number of lineages?
#' @return A matrix with one row per community and one column per order.
#' @references Chao, A., Chiu, C.-H. and Jost, L. (2010) Phylogenetic
#'   diversity measures based on Hill numbers. Philosophical Transactions
#'   of the Royal Society B 365, 3599-3609.
#' @examples
#' data(PD.example.sample, PD.example.phylo)
#' phylo.hill(PD.example.phylo, PD.example.sample, q = c(0, 1, 2))
#' phylo.hill(PD.example.phylo, PD.example.sample, q = 0:2, mean.pd = TRUE)
#' @export
phylo.hill <- function(my.phylo, my.sample, q = c(0, 1, 2), mean.pd = FALSE) {
	ra <- to.ra(my.sample)
	if (!all(colnames(ra) %in% my.phylo$tip.label)) {
		stop("all species in my.sample must be tips of my.phylo")
	}
	m <- branch.membership(my.phylo)[colnames(ra), , drop = FALSE]
	bl <- my.phylo$edge.length
	## depth of the tree: the longest root-to-tip path
	depth <- max(ape::node.depth.edgelength(my.phylo))
	## summed relative abundance under each branch, communities by branches
	a <- ra %*% m
	out <- sapply(q, function(qq) {
		apply(a, 1, function(ai) {
			ok <- ai > 0
			p <- ai[ok] / depth
			if (qq == 1) {
				exp(-sum(bl[ok] * p * log(p)))
			} else {
				sum(bl[ok] * p^qq)^(1 / (1 - qq))
			}
		})
	})
	out <- matrix(out, nrow = nrow(ra), dimnames = list(rownames(ra), paste0("q", q)))
	if (mean.pd) out <- out / depth
	out
}
