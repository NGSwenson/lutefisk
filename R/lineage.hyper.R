#' Hypergeometric test of lineage representation in assemblages
#'
#' For every internal node of a phylogeny for a global species pool, test
#' whether the number of that node's descendant species present in an
#' assemblage is larger or smaller than expected if the assemblage were a
#' random draw of species from the pool. This is the analytical form of the
#' node-by-node test that the \code{nodesig} routine of Phylocom performed
#' by randomization.
#'
#' If the pool contains \eqn{N} species of which \eqn{K} descend from the
#' focal node, and the assemblage contains \eqn{n} species, then the number
#' \eqn{k} of those \eqn{n} that descend from the node follows a
#' hypergeometric distribution, and the probability of \eqn{k} or more
#' (over-representation) or \eqn{k} or fewer (under-representation) comes
#' directly from \code{\link[stats]{phyper}}. The local counts for every
#' node and every assemblage are obtained by one matrix product of the
#' presence-absence matrix with the tip-by-node membership matrix, so the
#' test runs in a fraction of a second on trees with thousands of tips.
#'
#' Nested nodes are not independent: an over-represented clade will tend to
#' pull its parent node along with it, and with hundreds of nodes some
#' small probabilities are expected by chance, so the probabilities are best
#' read as a ranking of lineages rather than as a set of independent tests.
#'
#' @param my.phylo A \code{phylo} object for the global species pool.
#' @param my.sample A community data matrix with assemblages as rows and
#'   species as columns. Only presence is used.
#' @param min.size Ignore nodes subtending fewer than this many species in
#'   the pool.
#' @param alternative \code{"both"}, \code{"over"} or \code{"under"};
#'   determines the ordering of the output and, for \code{"both"}, adds the
#'   smaller of the two one-tailed probabilities as \code{p.min}.
#' @return A list with one data frame per assemblage. Each has one row per
#'   tested node with the node number, the global and local descendant
#'   counts, the expected local count, the one-tailed probabilities, the
#'   direction of the deviation and, when the tree has node labels, the
#'   label.
#' @examples
#' ## five of the ten species of a genus in a plot of 100 species drawn
#' ## from a pool of 1000
#' phyper(4, 10, 990, 100, lower.tail = FALSE)
#'
#' data(null.example.sample, null.example.phylo)
#' res <- lineage.hyper(null.example.phylo, null.example.sample)
#' head(res[[1]])
#' @export
lineage.hyper <- function(my.phylo, my.sample, min.size = 2,
	alternative = c("both", "over", "under")) {
	alternative <- match.arg(alternative)

	## Species present in the phylogeny and the sample, in phylogeny order.
	spp <- my.phylo$tip.label
	if (!all(colnames(my.sample) %in% spp)) {
		stop("all species in my.sample must be tips of my.phylo")
	}
	pa <- (my.sample > 0) * 1
	pa <- pa[, colnames(pa) %in% spp, drop = FALSE]

	N <- length(spp)
	n.tip <- N

	## Tip-by-node membership matrix: entry [i, j] is 1 when tip i
	## descends from internal node j.
	member <- node.membership(my.phylo)

	## Global descendant count for every node (K) and the nodes to test.
	K <- colSums(member)
	keep <- K >= min.size & K < N
	member <- member[, keep, drop = FALSE]
	K <- K[keep]

	## Local descendant counts for every node in every community (k), by a
	## single matrix product: communities by species times species by nodes.
	k <- pa[, spp, drop = FALSE] %*% member
	n <- rowSums(pa)

	## Hypergeometric probabilities. stats::phyper(q, m, n, k) gives P(X <= q) for
	## m white balls (K), n black balls (N - K) and k draws (n species).
	one.com <- function(i) {
		ki <- k[i, ]
		ni <- n[i]
		p.under <- stats::phyper(ki, K, N - K, ni)
		p.over <- stats::phyper(ki - 1, K, N - K, ni, lower.tail = FALSE)
		expected <- ni * K / N
		out <- data.frame(node = as.integer(colnames(member)),
			global.n = as.integer(K), local.n = as.integer(ki),
			expected = expected, p.over = p.over, p.under = p.under,
			direction = ifelse(ki > expected, "over",
				ifelse(ki < expected, "under", "none")),
			stringsAsFactors = FALSE)
		if (!is.null(my.phylo$node.label)) {
			out$label <- my.phylo$node.label[out$node - n.tip]
		}
		if (alternative == "over") out <- out[order(out$p.over), ]
		if (alternative == "under") out <- out[order(out$p.under), ]
		if (alternative == "both") {
			out$p.min <- pmin(out$p.over, out$p.under)
			out <- out[order(out$p.min), ]
		}
		rownames(out) <- NULL
		out
	}
	res <- lapply(seq_len(nrow(pa)), one.com)
	names(res) <- rownames(pa)
	res
}
