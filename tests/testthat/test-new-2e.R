test_that("exact null moments match randomization", {
  set.seed(5)
  tr <- ape::rcoal(40)
  d <- stats::cophenetic(tr)
  k <- 8
  pa <- t(replicate(5000, { x <- integer(40); x[sample(40, k)] <- 1L; x }))
  colnames(pa) <- tr$tip.label
  pd <- faith.pd(tr, pa)
  mpd <- Fmpd(d, pa)
  e.pd <- expected.pd(tr, k)
  e.mpd <- expected.mpd(d, k)
  expect_equal(mean(pd), e.pd$expected, tolerance = 0.03)
  expect_equal(sd(pd), e.pd$sd, tolerance = 0.05)
  expect_equal(mean(mpd), e.mpd$expected, tolerance = 0.03)
  expect_equal(sd(mpd), e.mpd$sd, tolerance = 0.05)
  expect_equal(e.mpd$expected, mean(d[upper.tri(d)]))
})

test_that("resolve.poly gives binary ultrametric trees with positive branches", {
  data(example.poly)
  r <- resolve.poly(example.poly, "equidistant")
  expect_true(ape::is.binary(r))
  expect_true(ape::is.ultrametric(r))
  expect_true(all(r$edge.length >= 0))
  expect_equal(max(ape::node.depth.edgelength(r)), max(ape::node.depth.edgelength(example.poly)))
  set.seed(1)
  rr <- resolve.poly(example.poly, "random", n = 20)
  expect_s3_class(rr, "multiPhylo")
  expect_true(all(sapply(rr, ape::is.ultrametric)))
  dd <- sapply(rr, function(x) stats::cophenetic(x)["sp1", "sp2"])
  expect_true(length(unique(round(dd, 6))) > 1)
  rd <- resolve.poly(example.poly, "directed", shape = c(1, 9), n = 20)
  expect_true(all(sapply(rd, function(x) all(x$edge.length >= 0))))
})

test_that("lineage.paint returns one row per node", {
  data(null.pd.example.phylo, null.pd.example.sample)
  cdm <- picante::sample2matrix(null.pd.example.sample)
  cdm <- cdm[, null.pd.example.phylo$tip.label[null.pd.example.phylo$tip.label %in% colnames(cdm)]]
  keep <- ape::keep.tip(null.pd.example.phylo, colnames(cdm))
  res <- lineage.paint(keep, cdm, community = 1, plot = FALSE)
  expect_equal(nrow(res), keep$Nnode)
})
