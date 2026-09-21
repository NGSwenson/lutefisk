#' Paint over- and under-represented lineages onto the phylogeny
#'
#' Runs \code{\link{lineage.hyper}} for one assemblage and draws the pool
#' phylogeny with each internal node colored by the direction and strength
#' of its deviation from the hypergeometric expectation, so that the
#' lineages an assemblage draws on more or less heavily than a random draw
#' from the pool can be seen at a glance. Node size scales with
#' \eqn{-\log_{10}} of the smaller one-tailed probability, warm colors mark
#' over-representation and cool colors under-representation, and nodes
#' whose probability exceeds \code{alpha} are left unmarked. The tips that
#' occur in the assemblage are marked as well.
#'
#' The function returns the per-node table invisibly, with a column of
#' colors and one of signed scores (\eqn{-\log_{10} p}, negative for
#' under-representation), so that the same information can be drawn with
#' \code{ggtree} by joining the table to the tree on \code{node}.
#'
#' @inheritParams lineage.hyper
#' @param community The row name or number of the assemblage to paint.
#' @param alpha Nodes with both one-tailed probabilities above this value
#'   are not marked.
#' @param plot Draw the tree, or only return the table.
#' @param tip.cex,node.cex Sizes of the tip markers and the largest node
#'   marker.
#' @param over.col,under.col Colors for over- and under-represented nodes.
#' @param ... Further arguments to \code{\link[ape]{plot.phylo}}.
#' @return Invisibly, a data frame with one row per internal node of
#'   \code{my.phylo}, in node order, containing the output of
#'   \code{lineage.hyper} plus \code{score} and \code{color}.
#' @examples
#' data(null.pd.example.phylo)
#' set.seed(2)
#' clade <- geiger::tips(null.pd.example.phylo, 55)
#' com <- c(sample(clade, min(6, length(clade))),
#'          sample(setdiff(null.pd.example.phylo$tip.label, clade), 4))
#' cdm <- matrix(as.integer(null.pd.example.phylo$tip.label %in% com), 1,
#'               dimnames = list("plot1", null.pd.example.phylo$tip.label))
#' res <- lineage.paint(null.pd.example.phylo, cdm, cex = 0.6)
#' head(res[order(res$p.min), ])
#' @export
lineage.paint <- function(my.phylo, my.sample, community = 1, min.size = 2,
	alpha = 0.05, plot = TRUE, tip.cex = 0.8, node.cex = 2.5,
	over.col = "firebrick", under.col = "steelblue", ...) {
	pa <- to.pa(my.sample)
	pa <- pa[community, , drop = FALSE]
	res <- lineage.hyper(my.phylo, pa, min.size = min.size, alternative = "both")[[1]]
	n.tip <- length(my.phylo$tip.label)
	## one row per internal node, in node order, NA where untested
	full <- data.frame(node = seq_len(my.phylo$Nnode) + n.tip)
	full <- merge(full, res, by = "node", all.x = TRUE, sort = TRUE)
	full$score <- ifelse(full$direction == "under", 1, -1) * log10(full$p.min)
	full$score[is.na(full$score)] <- 0
	sig <- !is.na(full$p.min) & full$p.min <= alpha
	full$color <- NA_character_
	full$color[sig & full$direction == "over"] <- over.col
	full$color[sig & full$direction == "under"] <- under.col
	if (plot) {
		ape::plot.phylo(my.phylo, ...)
		present <- my.phylo$tip.label %in% colnames(pa)[pa[1, ] > 0]
		ape::tiplabels(pch = 15, col = ifelse(present, "black", "grey80"),
			cex = tip.cex, adj = 0.5)
		if (any(sig)) {
			sz <- node.cex * pmin(abs(full$score[sig]), 4) / 4
			ape::nodelabels(node = full$node[sig], pch = 21, bg = full$color[sig],
				col = "black", cex = pmax(sz, 0.8))
		}
	}
	invisible(full)
}
