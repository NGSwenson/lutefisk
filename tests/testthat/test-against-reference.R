## Every hand-built function is checked against an independent
## implementation: the picante or FD package where one exists, and the
## original loop-based code of lefse 0.5 (kept in tests/testthat/lefse05.R)
## otherwise. This carries the book's argument, that you should not trust a
## number you have not computed a second way, into the test suite.

source(test_path("lefse05.R"))
library(ape)

data(lefse.sample, lefse.phylo, lefse.traits)
data(PD.example.sample, PD.example.phylo)
d <- cophenetic(lefse.phylo)
cdm <- as.matrix(lefse.sample)

test_that("Fmpd and Fmpd.a match picante::mpd", {
	skip_if_not_installed("picante")
	expect_equal(unname(Fmpd(d, cdm)), picante::mpd(cdm, d))
	expect_equal(unname(Fmpd.a(d, cdm)), picante::mpd(cdm, d, abundance.weighted = TRUE))
})

test_that("Fmntd and Fmntd.a match picante::mntd", {
	skip_if_not_installed("picante")
	expect_equal(unname(Fmntd(d, cdm)), picante::mntd(cdm, d))
	expect_equal(unname(Fmntd.a(d, cdm)), picante::mntd(cdm, d, abundance.weighted = TRUE))
})

test_that("nearest taxon functions match lefse 0.5", {
	expect_equal(Fmntd(d, cdm), lefse05$Fmntd(d, cdm))
	expect_equal(Fmntd.a(d, cdm), lefse05$Fmntd.a(d, cdm))
	expect_equal(Fsntd(d, cdm), lefse05$Fsntd(d, cdm))
	expect_equal(unname(sntd.a(d, cdm)), unname(lefse05$sntd.a(d, cdm)))
})

test_that("faith.pd matches picante::pd with and without the root", {
	skip_if_not_installed("picante")
	s <- as.matrix(PD.example.sample)
	expect_equal(unname(faith.pd(PD.example.phylo, s)),
		picante::pd(s, PD.example.phylo, include.root = TRUE)$PD)
	expect_equal(unname(faith.pd(PD.example.phylo, s, include.root = FALSE)),
		picante::pd(s, PD.example.phylo, include.root = FALSE)$PD)
})

test_that("weighted.faith matches lefse 0.5", {
	skip_if_not_installed("geiger")
	s <- as.matrix(PD.example.sample)
	expect_equal(weighted.faith(PD.example.phylo, s), lefse05$weighted.faith(PD.example.phylo, s))
})

test_that("trait moments match lefse 0.5 and SDMTools-style weighted sd", {
	expect_equal(com.trait.moments(cdm, lefse.traits), lefse05$com.trait.moments(cdm, lefse.traits))
	expect_equal(com.trait.weighted(cdm, lefse.traits), lefse05$com.trait.weighted(cdm, lefse.traits))
	expect_equal(trait.range(cdm, lefse.traits), t(lefse05$trait.range(cdm, lefse.traits)))
})

test_that("weighted.sd reproduces SDMTools::wt.sd", {
	x <- c(2.1, 3.4, 5.5, 1.2, 4.4); w <- c(1, 3, 2, 5, 1)
	xbar <- sum(w * x) / sum(w)
	ref <- sqrt(sum(w * (x - xbar)^2) * (sum(w) / (sum(w)^2 - sum(w^2))))
	expect_equal(weighted.sd(x, w), ref)
})

test_that("beta diversity matrices match picante::comdist and comdistnt", {
	skip_if_not_installed("picante")
	cd <- as.matrix(picante::comdist(cdm, d))
	cda <- as.matrix(picante::comdist(cdm, d, abundance.weighted = TRUE))
	cn <- as.matrix(picante::comdistnt(cdm, d))
	cna <- as.matrix(picante::comdistnt(cdm, d, abundance.weighted = TRUE))
	off <- upper.tri(cd)
	expect_equal(new.comdist(cdm, d)[off], cd[off])
	expect_equal(new.comdist.prime(cdm, d)[off], cda[off])
	expect_equal(new.comdistnn(cdm, d)[off], cn[off])
	expect_equal(new.comdistnn.prime(cdm, d)[off], cna[off])
})

test_that("beta diversity matrices match lefse 0.5", {
	expect_equal(new.comdist(cdm, d), lefse05$new.comdist(cdm, d), ignore_attr = TRUE)
	expect_equal(new.comdist.prime(cdm, d), lefse05$new.comdist.prime(cdm, d), ignore_attr = TRUE)
	expect_equal(new.comdistnn(cdm, d), lefse05$new.comdistnn(cdm, d), ignore_attr = TRUE)
	expect_equal(new.comdistnn.prime(cdm, d), lefse05$new.comdistnn.prime(cdm, d), ignore_attr = TRUE)
})

test_that("rank.abund.pd matches the loop over picante::pd", {
	skip_if_not_installed("picante")
	r <- rank.abund.pd(cdm, lefse.phylo)
	ref <- lefse05$rank.abund.pd(cdm, lefse.phylo)
	for (i in seq_along(r)) expect_equal(r[[i]], as.vector(ref[[i]]))
})

test_that("hull.volume matches geometry directly", {
	skip_if_not_installed("geometry")
	tr <- as.matrix(scale(lefse.traits))[, 1:2]
	spp <- names(cdm[4, ])[cdm[4, ] > 0]
	expect_equal(unname(hull.volume(cdm, tr)[4]),
		geometry::convhulln(tr[spp, ], options = "FA")$vol)
	expect_true(is.na(hull.volume(cdm, as.matrix(scale(lefse.traits)))[1]))
})

test_that("lineage.hyper agrees with phyper on a hand example and with randomization", {
	expect_equal(stats::phyper(4, 10, 990, 100, lower.tail = FALSE), 0.001521422, tolerance = 1e-6)
	set.seed(3)
	tr <- ape::rcoal(200)
	tr$tip.label <- paste0("s", 1:200)
	m <- node.membership(tr)
	node <- colnames(m)[which(colSums(m) == 8)[1]]
	clade <- rownames(m)[m[, node] == 1]
	spp <- c(sample(clade, 4), sample(setdiff(tr$tip.label, clade), 36))
	s <- matrix(0, 1, 200, dimnames = list("p", tr$tip.label)); s[1, spp] <- 1
	res <- lineage.hyper(tr, s)[[1]]
	p.an <- res$p.over[res$node == as.integer(node)]
	sim <- replicate(20000, sum(sample(tr$tip.label, 40) %in% clade) >= 4)
	expect_equal(p.an, mean(sim), tolerance = 0.15)
	expect_equal(res$local.n[res$node == as.integer(node)], 4L)
})
