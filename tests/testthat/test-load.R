test_that("the namespace does not import optional packages", {
  imports <- names(getNamespaceImports("GRACKLE"))
  expect_false(any(c("PLIER", "sgnesR", "NMF", "dynutils") %in% imports))
  expect_true(is.function(asNamespace("GRACKLE")$GRACKLE))
})

test_that("the compiled-code stubs are gone", {
  expect_false(exists("nmf_with_graph_regularization", envir = asNamespace("GRACKLE"), inherits = FALSE))
})
