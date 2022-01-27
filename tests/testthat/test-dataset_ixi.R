test_that("We can load the IXI metadata: demographics and brain stats", {
    md = load_IXI_metadata();

    testthat::expect_true(is.list(md));
    testthat::expect_true("brainstats" %in% names(md));
    testthat::expect_true("demographics" %in% names(md));
    testthat::expect_true("subjects_list" %in% names(md));
    testthat::expect_true("merged" %in% names(md));
})
