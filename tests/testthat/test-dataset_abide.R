
testthat::test_that("We can load the ABIDE metadata: demographics and brain stats", {
    md = load_ABIDE_metadata();

    testthat::expect_true(is.list(md));
    testthat::expect_true("brainstats" %in% names(md));
    testthat::expect_true("demographics" %in% names(md));
    testthat::expect_true("subjects_list" %in% names(md));
    testthat::expect_true("merged" %in% names(md));
})
